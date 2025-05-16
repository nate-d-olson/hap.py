#!/usr/bin/env python3
# coding=utf-8
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
#
# 9/9/2014
#
# Preprocessing for a VCF file
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import sys
import os
import argparse
import logging
import subprocess
import multiprocessing
import tempfile
import time
import pipes
from typing import List, Dict, Tuple, Optional, Union, Any, Set

scriptDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python3')))

import Tools
from Tools import vcfextract
from Tools.bcftools import preprocessVCF
from Tools.fastasize import fastaContigLengths
from Tools.bcftools import runBcftools
import Haplo.partialcredit


def hasChrPrefix(chrlist: List[str]) -> bool:
    """Returns if list of chr names has a chr prefix or not

    Args:
        chrlist: List of chromosome names

    Returns:
        True if the chromosomes have a 'chr' prefix, False otherwise
    """
    noprefix = [str(i) for i in range(23)] + ["X", "Y", "MT"]
    withprefix = ["chr" + x for x in [str(i) for i in range(23)] + ["X", "Y", "M"]]

    nwith = 0
    nwithout = 0

    for x in chrlist:
        if x in noprefix:
            nwithout += 1
        if x in withprefix:
            nwith += 1

    return nwith >= nwithout


def fixChrNames(chrlist: List[str], withchr: bool = True) -> List[str]:
    """Returns a list with fixed chromosome names

    Args:
        chrlist: List of chromosome names
        withchr: Whether to add or remove 'chr' prefix

    Returns:
        List of standardized chromosome names
    """
    result = []
    for c in chrlist:
        if withchr:
            # add chr
            if not c.startswith("chr"):
                c = "chr" + c

            # special case for MT
            if c == "chrMT":
                c = "chrM"
        else:
            # remove chr
            if c.startswith("chr"):
                c = c[3:]

            # special case for M
            if c == "M":
                c = "MT"

        result.append(c)
    return result


def runRScript(script: str) -> Tuple[str, str, int]:
    """Run an R script

    Args:
        script: R script to run

    Returns:
        Tuple of (stdout, stderr, return_code)
    """
    logging.info("Running R script %s " % script)

    temp = tempfile.NamedTemporaryFile(delete=False)
    f = temp.file
    f.write(script.encode('utf-8'))
    f.flush()
    f.close()

    cmd_line = "R --slave --args < %s" % temp.name
    logging.info(cmd_line)

    po = subprocess.Popen(cmd_line,
                          shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          universal_newlines=True)  # text mode output

    stdout, stderr = po.communicate()

    po.wait()

    return_code = po.returncode

    if return_code != 0:
        logging.error("Failed to run R script %s" % temp.name)
        # leave the temp file in case of an error
    else:
        logging.info("Removing temporary R script")
        os.unlink(temp.name)

    return stdout, stderr, return_code


def fixVCF(vcf_file: str, destination: str, locations: Optional[str] = None,
           reference: Optional[str] = None, bgzip: bool = True,
           threads: int = multiprocessing.cpu_count()) -> None:
    """Check and fix a VCF

    Args:
        vcf_file: Input VCF file
        destination: Output destination
        locations: Optional list of regions to include
        reference: Optional reference FASTA file
        bgzip: Whether to compress with bgzip
        threads: Number of threads to use
    """
    # TODO handle gzipped files
    if vcf_file.endswith(".gz"):
        orig = vcf_file
        tf = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False)
        tf.close()
        vcf_file = tf.name
        cmd_line = "gzip -dc %s > %s" % (pipes.quote(orig), pipes.quote(vcf_file))
        logging.info(cmd_line)

        po = subprocess.Popen(cmd_line,
                              shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        stdout, stderr = po.communicate()

        po.wait()

        return_code = po.returncode

        if return_code != 0:
            logging.error("gzip error: %s" % stderr)
            raise Exception("Failed to decompress %s" % orig)

    logging.info("Preprocessing VCF %s -> %s" % (vcf_file, destination))

    # if not reference:
    if not os.path.exists(destination) or os.path.getsize(destination) == 0:
        logging.info("Running preprocessVCF(%s, %s)" % (vcf_file, destination))
        preprocessVCF(vcf_file, destination, locations, reference, bgzip, threads)
    else:
        logging.info("Output file already exists.")

    if bgzip:
        cmd_line = "tabix -p vcf %s" % pipes.quote(destination)
        logging.info(cmd_line)

        po = subprocess.Popen(cmd_line,
                              shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        stdout, stderr = po.communicate()

        po.wait()

        return_code = po.returncode
        if return_code != 0:
            # try bcftools
            try:
                stdout, stderr, ret = runBcftools("index", destination)
            except Exception as e:
                logging.error("Tabix error: %s" % stderr)
                logging.error("BCFTools error: %s" % str(e))
                raise Exception("Failed to tabix / bcftools index output file %s" % destination)

    if vcf_file.endswith(".vcf") and vcf_file.startswith(tempfile.gettempdir()):
        logging.info("Cleaning up temp file %s" % vcf_file)
        try:
            os.unlink(vcf_file)
        except:
            pass


def main():
    parser = argparse.ArgumentParser("VCF file preprocessor. Normalises variants and fixes chromosome names.")

    parser.add_argument("input", help="VCF file to preprocess.")
    parser.add_argument("output", help="Destination file.")

    parser.add_argument("-r", "--reference", dest="ref", help="Reference fasta file.", default=None)
    parser.add_argument("-L", "--locations", dest="locations",
                        help="Locations to select, e.g. chr1:10000-20000.", default=None)
    parser.add_argument("-w", "--window", dest="window", help="Add +/- WindowSize (default=1kb) to locations.",
                        default=1000, type=int)
    parser.add_argument("-t", "--type", dest="types",
                        choices=["preprocess", "fixchr", "roc"],
                        help="Preprocessing type. preprocess (default) = all preprocessing steps. "
                             "fixchr = only correct the chromosome names. "
                             "roc = extract ROC info.", default="preprocess")
    parser.add_argument("--clear", dest="clear", help="Don't build on existing ROC file.", default=False,
                        action="store_true")
    parser.add_argument("--prefix-chr", dest="prefixchr", help="Add chr prefix.", default=None, action="store_true")
    parser.add_argument("--noprefix-chr", dest="prefixchr", help="Remove chr prefix.", action="store_false")
    parser.add_argument("-a", "--all-records", dest="all_records",
                        help="For ROC creation: "
                             "Allow non-PASS records and silently fix the VCF if some QQ scores are "
                             "missing when they should be present (e.g. for Platypus and FreeBayes)",
                        default=False, action="store_true")
    parser.add_argument("-T", "--truth", help="For ROC creation: truth file.", default=None)
    parser.add_argument("--filter", help="For ROC creation: filter to extract from VCF. "
                                         "Will try to recompute score if this is QUAL", default=None)
    parser.add_argument("--regions", help="For ROC creation: confident call regions.", default=None)
    parser.add_argument("--qual", help="For ROC creation: create QUAL field from tag", default=None)
    parser.add_argument("--threads", help="Number of threads to use.", default=multiprocessing.cpu_count(), type=int)

    args = parser.parse_args()

    if not args.ref:
        args.ref = Tools.defaultReference()

    if args.locations:
        # add a window to give more room for manipulations
        locs = []
        for p in args.locations.split(","):
            parts = p.split(":")
            if len(parts) == 1:
                locs.append(p)
                continue
            try:
                pos = parts[1].split("-")
                if len(pos) == 1:
                    locs.append("%s:%i-%i" % (parts[0], max(0, int(pos[0]) - args.window), int(pos[0]) + args.window))
                else:
                    locs.append("%s:%i-%i" % (parts[0], max(0, int(pos[0]) - args.window), int(pos[1]) + args.window))
            except:
                locs.append(p)
        args.locations = ",".join(locs)

    # preprocess file
    if args.types == "preprocess":
        fixVCF(args.input, args.output, args.locations, args.ref, True, args.threads)
    elif args.types == "fixchr":
        # Only adjust chromosome names, making sure we get all chromosomes right.
        # It is important to do this before preprocessing, becasue the BCF index creation will
        # potentially mess with the order in which chromosomes are written to the bcf.
        if args.prefixchr is None:
            # autodetect
            import vcf
            reader = None
            try:
                reader = vcf.Reader(filename=args.input)
                # note if this is never accessed, we end up with no contigs from reader
                contig_prefix = None

                # If bcf, just get the list of sequence names
                header_contigs = []
                try:
                    stdin, stdout, stderr = runBcftools("view", "-h", args.input)
                    has_bcf_header = not stderr and "##contig=<ID=" in stdout
                    if has_bcf_header:
                        for x in stdout.splitlines():
                            if "##contig=<ID=" in x:
                                name = x.split("##contig=<ID=")[1].split(",")[0].strip("\"'>")
                                header_contigs.append(name)
                    if header_contigs:
                        args.prefixchr = hasChrPrefix(header_contigs)
                    else:
                        for c in reader.contigs:
                            header_contigs.append(c)

                        if header_contigs:
                            args.prefixchr = hasChrPrefix(header_contigs)
                        else:
                            # try contig names from query variants
                            contigs = set()
                            for variant in reader:
                                contigs.add(str(variant.CHROM))
                                if len(contigs) >= 10:
                                    break
                            args.prefixchr = hasChrPrefix(contigs)
                except:
                    logging.warn("Failed to extract vcf header using bcftools/pyvcf.")

                # get reference contig naming scheme, make target match reference
                if args.ref:
                    rcontigs = fastaContigLengths(args.ref)
                    header_contigs = rcontigs.keys()
                    if header_contigs:
                        args.prefixchr = hasChrPrefix(header_contigs)

            except:
                pass
            finally:
                try:
                    del reader
                except:
                    pass
                try:
                    del contig_prefix
                except:
                    pass

        if args.prefixchr:
            logging.info("Adding chr prefix to chromosome names.")
        else:
            logging.info("Removing chr prefix from chromosome names.")

        with open(args.output, "w") as output_file:
            headerprinted = False
            stdin, stdout, stderr = runBcftools("view", "-h", args.input)

            for l in stdout.splitlines():
                if "##contig=<ID=" in l:
                    name = l.split("##contig=<ID=")[1].split(",")[0].strip("\"'>")
                    if args.prefixchr:
                        if not name.startswith("chr"):
                            name = "chr" + name
                        if name == "chrMT":
                            name = "chrM"
                    else:
                        if name.startswith("chr"):
                            name = name[3:]
                        if name == "M":
                            name = "MT"
                    l = l.split("##contig=<ID=")[0] + "##contig=<ID=" + name + "," + \
                        l.split("##contig=<ID=")[1].split(",", 1)[1]

                output_file.write(l)
                output_file.write("\n")
                headerprinted = True

            if not headerprinted:
                logging.error("Failed to extract header from VCF.")
                return 1

            stdin, stdout, stderr = runBcftools("view", args.input)

            for l in stdout.splitlines():
                if l.startswith("#"):
                    if not "##contig=" in l or not headerprinted:
                        output_file.write(l)
                        output_file.write("\n")
                else:
                    rec = l.split("\t")
                    if len(rec) >= 1:
                        if args.prefixchr:
                            if not rec[0].startswith("chr"):
                                rec[0] = "chr" + rec[0]
                            if rec[0] == "chrMT":
                                rec[0] = "chrM"
                        else:
                            if rec[0].startswith("chr"):
                                rec[0] = rec[0][3:]
                            if rec[0] == "M":
                                rec[0] = "MT"
                    output_file.write("\t".join(rec))
                    output_file.write("\n")

        if args.output.endswith(".gz"):
            # try to compress output file usign bgzip
            if not os.path.exists(args.output):
                output_bgzip = args.output
                nongz_output = args.output[:-3]  # trim .gz
                if os.path.exists(nongz_output):
                    cmd_line = "bgzip -f %s" % pipes.quote(nongz_output)
                    logging.info(cmd_line)

                    po = subprocess.Popen(cmd_line,
                                          shell=True,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE)

                    stdout, stderr = po.communicate()

                    po.wait()

                    return_code = po.returncode

                    if return_code != 0:
                        logging.error("bgzip error: %s" % stderr)
                        raise Exception("Failed to compress %s" % nongz_output)
        elif args.output.endswith(".bcf"):
            output_bcf = args.output
            nongz_output = args.output[:-4] + ".vcf"  # replace .bcf with .vcf
            if os.path.exists(nongz_output) and not os.path.exists(args.output):
                cmd_line = "bcftools view -O b -o %s %s" % (pipes.quote(output_bcf), pipes.quote(nongz_output))
                logging.info(cmd_line)

                po = subprocess.Popen(cmd_line,
                                      shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)

                stdout, stderr = po.communicate()

                po.wait()

                return_code = po.returncode

                if return_code != 0:
                    logging.error("bcftools view -O b error: %s" % stderr)
                    raise Exception("Failed to compress %s" % nongz_output)
                else:
                    logging.info("Removing temporary vcf file")
                    os.unlink(nongz_output)
    elif args.types == "roc":
        # TODO Find the version of truth and query
        try:
            truth_version = args.truth
            query_version = args.input
            try:
                if args.qual is not None:
                    stdin, stdout, stderr = runBcftools("annotate", "-x", "^INFO/" + args.qual,
                                                        "|", "bcftools", "query", "-f", "[%CHROM\t%POS\t%REF\t%ALT\t%" +
                                                        args.qual + "\\n]")
            except:
                logging.info("Creating QUAL from %s" % args.qual)
                # try to overwrite the qual field
                fd, tempvcf = tempfile.mkstemp(prefix="vcfqual", suffix=".vcf")
                os.close(fd)

                # write a new vcf to the tempfile
                from Tools.bcftools import runBcftoolsView

                try:
                    with open(args.input) as f:
                        with open(tempvcf, "w") as o:
                            stdin, stdout, stderr = runBcftoolsView("-h", args.input)
                            o.write(stdout)
                            stdin, stdout, stderr = runBcftoolsView(args.input)
                            for line in stdout.splitlines():
                                if not line.startswith("#"):
                                    rec = line.split("\t")
                                    if len(rec) >= 8:
                                        info = rec[7].split(";")
                                        for i in range(len(info)):
                                            if info[i] == args.qual:
                                                rec[5] = "."
                                            elif info[i].startswith(args.qual + "="):
                                                rec[5] = info[i].split("=")[1]
                                    o.write("\t".join(rec))
                                    o.write("\n")
                                else:
                                    o.write(line)
                                    o.write("\n")
                    args.input = tempvcf
                    logging.info("Success.")
                except Exception as e:
                    try:
                        os.unlink(tempvcf)
                    except:
                        pass
                    logging.error("Failed to preprocess VCF: %s " % str(e))
                    raise e
        except Exception as e:
            logging.warn("Cannot determine version: %s" % str(e))

        quantifyVCF(args.input, args.truth, args.output, args.filter, args.regions, args.all_records, args.clear)

    return 0


def quantifyVCF(query, truth, output, expr=None, regions=None, allow_non_pass=False,
                clear_existing=False, include_homref=False):
    """Extract quantification info from query VCF for ROC from TP / FP status.

    Args:
        query: Query VCF
        truth: Truth VCF
        output: Output file
        expr: Expression to filter on
        regions: Regions to restrict to
        allow_non_pass: Allow non-PASS records
        clear_existing: Clear existing output file
        include_homref: Include homozygous reference variants
    """
    # special case for expr == QUAL
    q_is_qual = False
    if expr == "QUAL":
        q_is_qual = True
        expr = None

    logging.info("Creating ROC information from %s" % query)
    vroc = vcfextract.VcfRocExtractor(query, not include_homref, True)
    vroc.setRegions(regions)
    vroc.setTruthVcf(truth)
    if expr:
        vroc.setInfoExpr(expr)
    elif not q_is_qual:
        vroc.useQual(True)

    roc_df = vroc.getROC(allow_non_pass, clear_existing)
    if len(roc_df) <= 0:
        logging.warn("VCF record set is empty.")

    if output.endswith(".csv") or output.endswith(".csv.gz"):
        oname = output
    else:
        oname = output + ".csv"
    import pandas
    pandas.set_option('io.hdf.default_format', 'table')
    if output.endswith(".gz"):
        oname = output
        output = output[:-3]
        roc_df.to_csv(output)
        cmd_line = "gzip -f %s" % pipes.quote(output)
        logging.info(cmd_line)

        po = subprocess.Popen(cmd_line,
                              shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
        stdout, stderr = po.communicate()
        po.wait()
        return_code = po.returncode

        if return_code != 0:
            logging.error("gzip error: %s" % stderr)
            raise Exception("Failed to compress %s" % output)
    else:
        if output.endswith(".h5") or output.endswith(".hdf"):
            oname = output
        else:
            oname = output + ".h5"
        roc_df.to_hdf(oname, "roc", complevel=9, complib="zlib")
        roc_df.to_csv(oname + ".csv")


if __name__ == "__main__":
    sys.exit(main())

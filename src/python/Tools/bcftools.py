#!/usr/bin/env python
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

import os
import subprocess
import logging
import pandas
import tempfile
import gzip
import pipes

import Tools


scriptDir = os.path.abspath(os.path.dirname(__file__))


def runShellCommand(*args):
    """ Run a shell command (e.g. bcf tools), and return output
    """
    cmdstr = " ".join(args)
    logging.info("CMD: %s" % cmdstr)

    pipe = subprocess.Popen(cmdstr, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = pipe.communicate()

    # decode bytes to string in Python 3
    if isinstance(stdout, bytes):
        stdout = stdout.decode('utf-8')
    if isinstance(stderr, bytes):
        stderr = stderr.decode('utf-8')

    if pipe.returncode != 0:
        logging.error("Return code: %i from command: %s \nStdErr: %s" %
                     (pipe.returncode, cmdstr, stderr))
        raise Exception("Cannot run command %s" % cmdstr)

    for l in stderr.split("\n"):
        if l and l.strip() and not l.startswith("bcftools::") and not l.startswith("Lines"):
            logging.warning(l)

    return stdout


def runBcftools(*args):
    return runShellCommand("bcftools", *args)


def parseStats(output, colname="count"):
    """ Parse BCFTOOLS Stats Output """

    result = {}
    for x in output.split("\n"):
        if x.startswith("SN"):
            vx = x.split("\t")
            name = vx[2].replace("number of ", "").replace(":", "")
            count = int(vx[3])
            result[name] = count

    result = pandas.DataFrame(list(result.items()), columns=["type", colname])
    return result


def countVCFRows(filename):
    """ Count the number of rows in a VCF
    :param filename: VCF file name
    :return: number of rows
    """
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:  # Use text mode with 'rt'
            count = 0
            for s in f:
                if not s.startswith("#"):
                    count += 1
            return count
    else:
        with open(filename, "r") as f:
            count = 0
            for s in f:
                if not s.startswith("#"):
                    count += 1
            return count


def concatenateParts(output, *args):
    if len(args) == 0:
        raise Exception("No parts to concatenate.")

    if len(args) == 1:
        # single part, just copy
        if args[0] != output:
            runShellCommand("cp", args[0], output)
        return

    to_delete = []

    outputext = ".bcf"
    if output.endswith(".vcf.gz"):
        outputext = ".vcf.gz"

    try:
        if len(args) <= 10:
            vargs = ["bcftools", "concat"]
            vargs += args
            vargs += ["-o", output, "-O", "b" if outputext == ".bcf" else "z"]

            runShellCommand(*vargs)
        else:
            # block in chunks (TODO: make parallel)
            tf1 = tempfile.NamedTemporaryFile(suffix=outputext, delete=False)
            tf2 = tempfile.NamedTemporaryFile(suffix=outputext, delete=False)
            to_delete.append(tf1.name)
            to_delete.append(tf2.name)
            to_delete.append(tf1.name + ".csi")
            to_delete.append(tf2.name + ".csi")
            half1 = [tf1.name] + list(args[:len(args)//2])
            half2 = [tf2.name] + list(args[len(args)//2:])
            concatenateParts(*half1)
            runBcftools("index", tf1.name)
            concatenateParts(*half2)
            runBcftools("index", tf2.name)
            concatenateParts(output, tf1.name, tf2.name)
    finally:
        for f in to_delete:
            try:
                os.unlink(f)
            except:
                pass


# noinspection PyShadowingBuiltins
def preprocessVCF(input_filename, output_filename, location="",
                  pass_only=True,
                  chrprefix=True, norm=False,
                  regions=None, targets=None,
                  reference="fake_reference_path",
                  filters_only=None,
                  somatic_allele_conversion=False,
                  sample="SAMPLE",
                  filter_nonref=True,
                  convert_gvcf=False,
                  num_threads=4):
    """ Preprocess a VCF file with bcftools.

    :param input_filename: name of the input file
    :param output_filename: name of output file
    :param location: target location. Can be 'chr:pos' or 'chr'. If set to an empty string, all records will be
    considered.
    :param pass_only: filter to variant that PASS
    :param chrprefix: add 'chr' prefixes to the VCF records
    :param norm: normalize indels (bcftools norm)
    :param regions: regions bed file name.
    :param targets: targets bed file name.
    :param reference: reference sequence file name.
    :param filters_only: keep variants with a specific filter. Useful for VQSR when we want to see LowQuality records.
    :param somatic_allele_conversion: convert somatic variants. False, "het", "hom", or "hemi"
    :param sample: output sample name
    :param filter_nonref: filter out variants with <NON_REF> in their ALT column
    :param convert_gvcf: convert gvcf to vcf
    :param num_threads: how many threads to use
    :return: name of the resulting VCF file
    """

    tff = tempfile.NamedTemporaryFile(suffix="merged.vcf.gz", delete=False)
    tff.close()

    try:
        vargs = ["bcftools", "view"]

        if pass_only and not filters_only:
            vargs += ["-f", "PASS,.", "-e", "'N_ALT==0'"]
        elif filters_only:
            vargs += ["-f", filters_only, "-e", "'N_ALT==0'"]

        vargs += ["-o", tff.name, "-O", "z", "--threads", str(num_threads)]

        # apply filtering
        # ploidy column is expected to be appended to the sample column in vcfutils
        vargs += [input_filename]

        if filter_nonref:
            vargs += ["|", "grep", "-v", "NON_REF", "|", "bcftools", "view", "-o", tff.name, "-O", "z"]

        # replace chr if required
        if chrprefix and not somatic_allele_conversion:
            runShellCommand(*vargs)
            vargs = ["zcat", tff.name]

            if not norm:
                vargs += ["|", "perl", "-pe", "s/^([0-9XYM])/chr$1/", "|", "perl", "-pe", "s/chrMT/chrM/", "|", "bcftools",
                          "view"]

        if targets:
            vargs += ["-T", targets, "|", "bcftools", "view"]

        if location:
            vargs += ["-t", location, "|", "bcftools", "view"]

        if output_filename.endswith("vcf.gz"):
            int_suffix = "vcf.gz"
        else:
            int_suffix = ".bcf"

        # anything needs tabix? if so do an intermediate stage where we
        # index first
        if regions:
            if int_suffix == "vcf.gz":
                vargs += ["-o", tff.name, "-O", "z"]
                runShellCommand(*vargs)
                runShellCommand('bcftools', "index", "-t", tff.name)
            else:
                vargs += ["-o", tff.name, "-O", "b"]
                runShellCommand(*vargs)
                runShellCommand('bcftools', "index", tff.name)
            vargs = ["bcftools", "view", tff.name, "-R", regions]

        # do normalization ?
        if norm:
            vargs += ["|", "bcftools", "norm", "-f", reference]

        # convert somatic alleles? use custom tool
        if somatic_allele_conversion:
            allele_conv_type = str(somatic_allele_conversion)
            if allele_conv_type == "True":
                allele_conv_type = "half"

            vargs += ["|", os.path.join(scriptDir, "..", "c++", "bin", "somaticAlleleConversion"),
                      "--conversion", allele_conv_type,
                      "--sample-name", sample]

        # convert gVCF?
        if convert_gvcf:
            vargs += ["|", os.path.join(scriptDir, "..", "bin", "gvcfgenotyper")]

        # add output format
        if output_filename.endswith("vcf.gz"):
            vargs += ["-O", "z"]
            istabix = True
        else:
            vargs += ["-O", "b"]
            istabix = False

        runShellCommand(*vargs)

        if istabix:
            runShellCommand('bcftools', "index", "-t", output_filename)
        else:
            runShellCommand('bcftools', "index", output_filename)

    except Exception as ex:
        print("Error running BCFTOOLS; please check your file for compatibility issues issues using vcfcheck")
        raise ex

    finally:
        try:
            os.unlink(tff.name)
        except:
            pass
        try:
            os.unlink(tff.name + ".tbi")
        except:
            pass
        try:
            os.unlink(tff.name + ".csi")
        except:
            pass


def bedOverlapCheck(filename):
    """ Check for overlaps / out of order in a bed file """
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:  # Use text mode with 'rt'
            last = -1
            lines = 1
            thischr = None
            for line in f:
                l = line.split("\t")
                if len(l) < 3:
                    continue
                if thischr is not None and thischr != l[0]:
                    last = -1
                thischr = l[0]
                if (last-1) > int(l[1]):
                    logging.warn("%s has overlapping regions at %s:%i (line %i)" % (filename, l[0], int(l[1]), lines))
                    return 1
                last = int(l[2])
                lines += 1
            return 0
    else:
        with open(filename, "r") as f:
            last = -1
            lines = 1
            thischr = None
            for line in f:
                l = line.split("\t")
                if len(l) < 3:
                    continue
                if thischr is not None and thischr != l[0]:
                    last = -1
                thischr = l[0]
                if (last-1) > int(l[1]):
                    logging.warn("%s has overlapping regions at %s:%i (line %i)" % (filename, l[0], int(l[1]), lines))
                    return 1
                last = int(l[2])
                lines += 1
            return 0

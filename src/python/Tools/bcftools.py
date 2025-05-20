#!/usr/bin/env python3
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt

import contextlib
import gzip
import logging
import os
import pipes
import subprocess
import tempfile
from typing import Tuple

import pandas

scriptDir = os.path.abspath(os.path.dirname(__file__))


def runShellCommand(*args: str) -> Tuple[str, str, int]:
    """Run a shell command (e.g. bcf tools), and return output

    Args:
        *args: Command arguments

    Returns:
        Tuple of (stdout, stderr, return_code)
    """
    qargs = []
    for a in args:
        if a.strip() != "|":
            qargs.append(pipes.quote(a))
        else:
            qargs.append("|")

    cmd_line = " ".join(qargs)
    logging.info(cmd_line)

    po = subprocess.Popen(
        cmd_line,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )  # text mode for Python 3 compatibility

    stdout, stderr = po.communicate()

    po.wait()

    return_code = po.returncode

    if return_code != 0:
        raise Exception(
            "Command line {} got return code {}.\nSTDOUT: {}\nSTDERR: {}".format(
                cmd_line, return_code, stdout, stderr
            )
        )

    return stdout, stderr, return_code


def runBcftools(*args: str) -> str:
    """
    Wraps runShellCommand for compatibility.

    Args:
        *args: Arguments to pass to bcftools

    Returns:
        stdout from the command
    """
    stdout, stderr, return_code = runShellCommand("bcftools", *args)
    return stdout


def parseStats(output, colname="count"):
    """Parse BCFTOOLS Stats Output"""

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
    """Count the number of rows in a VCF
    :param filename: VCF file name
    :return: number of rows
    """
    if filename.endswith(".gz"):
        f = gzip.open(filename, "rt", encoding="utf-8")  # text mode in Python 3
    else:
        f = open(filename, encoding="utf-8")

    count = 0
    for s in f:
        if not s.startswith("#"):
            count += 1

    f.close()
    return count


def concatenateParts(output, *args):
    """Concatenate BCF files


    Trickier than it sounds because when there are many files we might run into
    various limits like the number of open files, or the length of a command line.

    This function will bcftools concat in a tree-like fashion to avoid this.
    """
    to_delete = []
    try:
        if output.endswith(".bcf"):
            outputformat = "b"
            outputext = ".bcf"
        else:
            outputformat = "z"
            outputext = ".vcf.gz"
        if len(args) < 10:
            for x in args:
                if not os.path.exists(x + ".tbi") and not os.path.exists(x + ".csi"):
                    to_delete.append(x + ".csi")
                    runBcftools("index", "-f", x)
            cmdlist = ["concat", "-a", "-O", outputformat, "-o", output, *list(args)]
            runBcftools(*cmdlist)
        else:
            # block in chunks (TODO: make parallel)
            tf1 = tempfile.NamedTemporaryFile(suffix=outputext, delete=False)
            tf2 = tempfile.NamedTemporaryFile(suffix=outputext, delete=False)
            to_delete.append(tf1.name)
            to_delete.append(tf2.name)
            to_delete.append(tf1.name + ".csi")
            to_delete.append(tf2.name + ".csi")
            mid_point = (
                len(args) // 2
            )  # Use integer division for Python 3 compatibility
            half1 = [tf1.name, *list(args[:mid_point])]
            half2 = [tf2.name, *list(args[mid_point:])]
            concatenateParts(*half1)
            runBcftools("index", tf1.name)
            concatenateParts(*half2)
            runBcftools("index", tf2.name)
            concatenateParts(output, tf1.name, tf2.name)
    finally:
        for f in to_delete:
            with contextlib.suppress(Exception):
                os.unlink(f)


# noinspection PyShadowingBuiltins
def preprocessVCF(
    input_filename,
    output_filename,
    location="",
    pass_only=True,
    chrprefix=True,
    norm=False,
    regions=None,
    targets=None,
    reference="fake_reference_path",
    filters_only=None,
    somatic_allele_conversion=False,
    sample="SAMPLE",
    filter_nonref=True,
    convert_gvcf=False,
    num_threads=4,
):
    """Preprocess a VCF + create index

    :param input_filename: the input VCF / BCF / ...
    :param output_filename: the output VCF
    :param location: optional location string -- comma separated
    :param pass_only: only return passing variants
    :param chrprefix: fix chromosome prefix
    :param norm: run through bcftools norm to leftshift indels
    :param regions: specify a subset of regions (traversed using tabix index, which must exist)
    :param targets: specify a subset of target regions (streaming traversal)
    :param reference: reference fasta file to use
    :param filters_only: require a set of filters (overridden by pass_only)
    :param somatic_allele_conversion: assume the input file is a somatic call file and squash
                                      all columns into one, putting all FORMATs into INFO
                                      This is used to treat Strelka somatic files
                                      Possible values for this parameter:
                                      True [="half"] / "hemi" / "het" / "hom" / "half"
                                      to assign one of the following genotypes to the
                                      resulting sample:  1 | 0/1 | 1/1 | ./1
    :param sample: name of the output sample column when using somatic_allele_conversion
    :param filter_nonref: remove any variants genotyped as <NON_REF>
    """
    # ToDo: refactor for simplicity and performance

    vargs = [
        "bcftools",
        "view",
        "--threads",
        str(num_threads),
        "-O",
        "v",
        input_filename,
        "|",
    ]

    if convert_gvcf:
        # prefilter for variant sites (all genome VCF sites have a NON_REF allele, hence use 2 here)
        vargs += ["bcftools", "view", "-I", "-e", "N_ALT < 2", "-O", "u", "|"]
        # strip uninteresting details and arrays which prevent allele trimming
        vargs += [
            "bcftools",
            "annotate",
            "-x",
            "INFO,^FORMAT/GT,FORMAT/DP,FORMAT/GQ",
            "-O",
            "u",
            "|",
        ]
        # trim missing alleles, don't compute the AD/AF fields
        vargs += ["bcftools", "view", "-a", "-I", "-O", "u", "|"]
        # remove variants with NON_REF alleles, don't compute the AD/AF fields
        vargs += [
            "bcftools",
            "view",
            "-I",
            "-e",
            'ALT[*] = "<NON_REF>"',
            "-O",
            "v",
            "|",
        ]

    vargs += ["bcftools", "view", "-O", "v"]

    if filter_nonref:
        vargs += [
            "|",
            "python",
            f"{scriptDir}/remove_nonref_gt_variants.py",
            "|",
            "bcftools",
            "view",
            "-O",
            "v",
        ]

    if type(location) is list:
        location = ",".join(location)

    if pass_only:
        vargs += ["-f", "PASS,."]
    elif filters_only:
        vargs += ["-f", filters_only]

    if chrprefix:
        vargs += [
            "|",
            "perl",
            "-pe",
            "s/^([0-9XYM])/chr$1/",
            "|",
            "perl",
            "-pe",
            "s/chrMT/chrM/",
            "|",
            "bcftools",
            "view",
        ]

    if targets:
        vargs += ["-T", targets, "|", "bcftools", "view"]

    if location:
        vargs += ["-t", location, "|", "bcftools", "view"]

    int_suffix = "vcf.gz" if output_filename.endswith("vcf.gz") else ".bcf"

    tff = tempfile.NamedTemporaryFile(delete=False, suffix=int_suffix)

    try:
        # anything needs tabix? if so do an intermediate stage where we
        # index first
        if regions:
            if int_suffix == "vcf.gz":
                vargs += ["-o", tff.name, "-O", "z"]
                _, _, _ = runShellCommand(*vargs)
                _, _, _ = runShellCommand("bcftools", "index", "-t", tff.name)
            else:
                vargs += ["-o", tff.name, "-O", "b"]
                _, _, _ = runShellCommand(*vargs)
                _, _, _ = runShellCommand("bcftools", "index", tff.name)
            vargs = ["bcftools", "view", tff.name, "-R", regions]

        if somatic_allele_conversion:
            if type(somatic_allele_conversion) is not str:
                somatic_allele_conversion = "half"
            vargs += [
                "|",
                "alleles",
                "-",
                "-o",
                "-.vcf",
                "--gt",
                somatic_allele_conversion,
                "--sample",
                sample,
                "|",
                "bcftools",
                "view",
            ]

        if norm:
            vargs += ["|", "bcftools", "norm", "-f", reference, "-c", "x", "-D"]

        vargs += ["-o", output_filename]
        if int_suffix == "vcf.gz":
            vargs += ["-O", "z"]
            istabix = True
        else:
            vargs += ["-O", "b"]
            istabix = False

        _, _, _ = runShellCommand(*vargs)

        if istabix:
            _, _, _ = runShellCommand("bcftools", "index", "-t", output_filename)
        else:
            _, _, _ = runShellCommand("bcftools", "index", output_filename)

    except Exception as ex:
        print(
            "Error running BCFTOOLS; please check your file for compatibility issues issues using vcfcheck"
        )
        raise ex

    finally:
        with contextlib.suppress(Exception):
            os.unlink(tff.name)

        with contextlib.suppress(Exception):
            os.unlink(tff.name + ".tbi")

        with contextlib.suppress(Exception):
            os.unlink(tff.name + ".csi")


def bedOverlapCheck(filename):
    """Check for overlaps / out of order in a bed file"""
    if filename.endswith(".gz"):
        f = gzip.open(filename, "rt", encoding="utf-8")  # text mode in Python 3
    else:
        f = open(filename, encoding="utf-8")
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
        if (last - 1) > int(l[1]):
            logging.warn(
                "%s has overlapping regions at %s:%i (line %i)"
                % (filename, l[0], int(l[1]), lines)
            )
            return 1
        last = int(l[2])
        lines += 1
    return 0


def run_command(cmd_line, shell=True, fail_message=None):
    """Unified command execution with better error handling

    Args:
        cmd_line: Command line to execute (string or list)
        shell: Whether to use shell execution
        fail_message: Custom error message on failure

    Returns:
        Tuple of (stdout, stderr) as strings

    Raises:
        Exception: If command execution fails
    """
    logging.info(cmd_line if isinstance(cmd_line, str) else " ".join(cmd_line))

    try:
        process = subprocess.Popen(
            cmd_line,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,  # Return strings, not bytes
        )

        stdout, stderr = process.communicate()
        return_code = process.returncode

        if return_code != 0:
            error_msg = fail_message or f"Command failed with return code {return_code}"
            logging.error(f"{error_msg}\nSTDOUT: {stdout}\nSTDERR: {stderr}")
            raise Exception(f"{error_msg}: {stderr}")

        return stdout, stderr
    except Exception as e:
        error_msg = fail_message or "Command execution failed"
        logging.error(f"{error_msg}: {str(e)}")
        raise Exception(f"{error_msg}: {str(e)}")

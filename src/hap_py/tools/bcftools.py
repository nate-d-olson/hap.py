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
import shlex
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd

# Type aliases for better readability
FilePath = Union[str, Path]
CommandOutput = Tuple[str, str, int]  # stdout, stderr, return_code

scriptDir = os.path.abspath(os.path.dirname(__file__))


def runShellCommand(*args: str) -> CommandOutput:
    """Run a shell command (e.g. bcf tools), and return output

    Args:
        *args: Command arguments

    Returns:
        Tuple of (stdout, stderr, return_code)

    Raises:
        Exception: If command execution fails with non-zero return code
    """
    qargs = []
    for a in args:
        if a.strip() != "|":
            qargs.append(shlex.quote(a))
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
            f"Command line {cmd_line} got return code {return_code}.\nSTDOUT: {stdout}\nSTDERR: {stderr}"
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


def parseStats(output: str, colname: str = "count") -> pd.DataFrame:
    """Parse BCFTOOLS Stats Output.

    Args:
        output: String output from bcftools stats
        colname: Name for the count column

    Returns:
        DataFrame with parsed statistics
    """
    result: Dict[str, int] = {}
    for x in output.split("\n"):
        if x.startswith("SN"):
            vx = x.split("\t")
            name = vx[2].replace("number of ", "").replace(":", "")
            count = int(vx[3])
            result[name] = count

    result_df = pd.DataFrame(list(result.items()), columns=["type", colname])
    return result_df


def countVCFRows(filename: str) -> int:
    """Count the number of rows in a VCF.

    Args:
        filename: VCF file name

    Returns:
        Number of rows in the VCF
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


def concatenateParts(output: FilePath, *args: FilePath) -> None:
    """Concatenate BCF files.

    Trickier than it sounds because when there are many files we might run into
    various limits like the number of open files, or the length of a command line.

    This function will bcftools concat in a tree-like fashion to avoid this.

    Args:
        output: Output file path
        *args: Input file paths to concatenate

    Raises:
        Exception: If bcftools concat fails
    """
    to_delete: List[str] = []
    output_str = str(output)  # Convert Path to string if needed

    # Validate input files exist
    for arg in args:
        arg_str = str(arg)
        if not os.path.exists(arg_str):
            raise FileNotFoundError(f"Input file not found: {arg_str}")

    try:
        if output_str.endswith(".bcf"):
            outputformat = "b"
            outputext = ".bcf"
        else:
            outputformat = "z"
            outputext = ".vcf.gz"

        if len(args) < 10:
            # For a small number of files, concatenate directly
            for x in args:
                x_str = str(x)
                if not os.path.exists(x_str + ".tbi") and not os.path.exists(
                    x_str + ".csi"
                ):
                    idx_file = x_str + ".csi"
                    to_delete.append(idx_file)
                    try:
                        logging.info(f"Indexing file {x_str}")
                        runBcftools("index", "-f", x_str)
                    except Exception as e:
                        logging.error(f"Failed to index file {x_str}: {str(e)}")
                        raise

            cmdlist = ["concat", "-a", "-O", outputformat, "-o", output_str]
            cmdlist.extend([str(arg) for arg in args])

            try:
                logging.info(f"Concatenating {len(args)} files into {output_str}")
                runBcftools(*cmdlist)
            except Exception as e:
                logging.error(f"Failed to concatenate files: {str(e)}")
                raise
        else:
            # For many files, use recursive divide-and-conquer approach
            logging.info(
                f"Concatenating {len(args)} files with divide-and-conquer strategy"
            )

            # Create temporary files for intermediate results
            try:
                tf1 = tempfile.NamedTemporaryFile(suffix=outputext, delete=False)
                tf2 = tempfile.NamedTemporaryFile(suffix=outputext, delete=False)
                to_delete.extend(
                    [tf1.name, tf2.name, tf1.name + ".csi", tf2.name + ".csi"]
                )

                # Split input files into two halves
                mid_point = (
                    len(args) // 2
                )  # Integer division for Python 3 compatibility
                half1 = [tf1.name] + [str(arg) for arg in args[:mid_point]]
                half2 = [tf2.name] + [str(arg) for arg in args[mid_point:]]

                # Process each half recursively
                concatenateParts(*half1)
                try:
                    runBcftools("index", tf1.name)
                except Exception as e:
                    logging.error(f"Failed to index first half: {str(e)}")
                    raise

                concatenateParts(*half2)
                try:
                    runBcftools("index", tf2.name)
                except Exception as e:
                    logging.error(f"Failed to index second half: {str(e)}")
                    raise

                # Merge the two halves
                concatenateParts(output_str, tf1.name, tf2.name)
            except Exception as e:
                logging.error(f"Error during recursive concatenation: {str(e)}")
                raise
    except Exception as e:
        logging.error(f"Failed to concatenate VCF/BCF files: {str(e)}")
        # Re-raise to propagate the exception
        raise
    finally:
        # Clean up temporary files
        for f in to_delete:
            if os.path.exists(f):
                try:
                    os.unlink(f)
                except Exception as e:
                    logging.warning(f"Failed to delete temporary file {f}: {str(e)}")


# noinspection PyShadowingBuiltins
def preprocessVCF(
    input_filename: str,
    output_filename: str,
    location: Union[str, List[str]] = "",
    pass_only: bool = True,
    chrprefix: Optional[bool] = True,
    norm: bool = False,
    regions: Optional[str] = None,
    targets: Optional[str] = None,
    reference: str = "fake_reference_path",
    filters_only: Optional[str] = None,
    somatic_allele_conversion: Union[bool, str] = False,
    sample: str = "SAMPLE",
    filter_nonref: bool = True,
    convert_gvcf: bool = False,
    num_threads: int = 4,
) -> None:
    """Preprocess a VCF + create index

    Args:
        input_filename: The input VCF / BCF / ...
        output_filename: The output VCF
        location: Optional location string -- comma separated
        pass_only: Only return passing variants
        chrprefix: Fix chromosome prefix
        norm: Run through bcftools norm to leftshift indels
        regions: Specify a subset of regions (traversed using tabix index, which must exist)
        targets: Specify a subset of target regions (streaming traversal)
        reference: Reference fasta file to use
        filters_only: Require a set of filters (overridden by pass_only)
        somatic_allele_conversion: Assume the input file is a somatic call file and squash
                                   all columns into one, putting all FORMATs into INFO.
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


def bedOverlapCheck(filename: str) -> int:
    """Check for overlaps / out of order in a bed file

    Args:
        filename: Path to the BED file to check

    Returns:
        0 if no overlaps, 1 if overlaps found
    """
    if filename.endswith(".gz"):
        f = gzip.open(filename, "rt", encoding="utf-8")  # text mode in Python 3
    else:
        f = open(filename, encoding="utf-8")
    last = -1
    lines = 1
    thischr = None
    for line in f:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        if thischr is not None and thischr != parts[0]:
            last = -1
        thischr = parts[0]
        if (last - 1) > int(parts[1]):
            logging.warning(
                "%s has overlapping regions at %s:%i (line %i)"
                % (filename, parts[0], int(parts[1]), lines)
            )
            return 1
        last = int(parts[2])
        lines += 1
    return 0


def run_command(
    cmd_line: Union[str, List[str]],
    shell: bool = True,
    fail_message: Optional[str] = None,
) -> Tuple[str, str]:
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
        raise Exception(f"{error_msg}: {str(e)}") from e

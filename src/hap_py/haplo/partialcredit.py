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
# 3/9/2015
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

"""
Module for partial credit variant comparison processing.
Allows for preprocessing and integration of partial credit comparisons.
Modernized to remove blocksplit dependency - processes entire files directly.
"""

import contextlib
import itertools
import logging
import os
import shlex
import subprocess
import tempfile
import time
from typing import Any, Dict, List, Optional, Tuple, Union

from ..tools.bcftools import concatenateParts, runBcftools
from ..tools.parallel import getPool, runParallel
from ..tools.vcfextract import extractHeadersJSON


def preprocessWrapper(
    file_and_location: Tuple[str, str], args: Dict[str, Any]
) -> Optional[str]:
    """Process a VCF file with the preprocess tool.

    Args:
        file_and_location: Tuple of (filename, location_str)
        args: Arguments for preprocessing

    Returns:
        Path to the preprocessed output file or None if processing failed
    """
    starttime = time.time()
    filename, location_str = file_and_location
    int_suffix = "bcf" if args["bcf"] else "vcf.gz"
    temp_file_path = None

    try:
        with tempfile.NamedTemporaryFile(
            delete=False, prefix=f"input.{location_str}", suffix=f".prep.{int_suffix}"
        ) as tf:
            temp_file_path = tf.name  # Store the file path for later use or cleanup

        # Quote filenames for shell safety
        quoted_filename = shlex.quote(filename)
        quoted_location = shlex.quote(location_str) if location_str else ""
        quoted_reference = shlex.quote(args["reference"])

        location_param = f"-l {quoted_location} " if location_str else ""

        to_run = (
            f"preprocess {quoted_filename}:* {location_param}-o {temp_file_path} "
            f"-V {args['decompose']} -L {args['leftshift']} -r {quoted_reference}"
        )

        if args["haploid_x"]:
            to_run += " --haploid-x 1"

        with (
            tempfile.NamedTemporaryFile(
                delete=False, prefix="stderr", suffix=".log"
            ) as tfe,
            tempfile.NamedTemporaryFile(
                delete=False, prefix="stdout", suffix=".log"
            ) as tfo,
        ):
            finished = False
            try:
                logging.info(f"Running '{to_run}'")
                subprocess.check_call(shlex.split(to_run), stdout=tfo, stderr=tfe)
                finished = True
            finally:
                if finished:
                    # Close files and read their contents for logging
                    with open(tfo.name, encoding="utf-8") as file:
                        for line in file:
                            logging.info(line.rstrip())
                    os.unlink(tfo.name)

                    with open(tfe.name, encoding="utf-8") as file:
                        for line in file:
                            logging.warning(line.rstrip())
                    os.unlink(tfe.name)
                else:
                    logging.error(
                        f"Preprocess command {to_run} failed. "
                        f"Outputs are here {tfo.name} / {tfe.name}"
                    )
                    with open(tfo.name, encoding="utf-8") as file:
                        for line in file:
                            logging.error(line.rstrip())
                    with open(tfe.name, encoding="utf-8") as file:
                        for line in file:
                            logging.error(line.rstrip())

                    # Cleanup the temp file if command failed
                    if temp_file_path and os.path.exists(temp_file_path):
                        os.unlink(temp_file_path)
                    return None  # Return None to indicate failure

        elapsed = time.time() - starttime
        logging.info(f"preprocess for {location_str} -- time taken {elapsed:.2f}")

        # Index the output file
        try:
            runBcftools("index", temp_file_path)
        except Exception as e:
            logging.error(f"Failed to index file {temp_file_path}: {str(e)}")
            if temp_file_path and os.path.exists(temp_file_path):
                os.unlink(temp_file_path)
            return None

        return temp_file_path
    except Exception as e:
        # Catch any unexpected exceptions and log them
        logging.error(f"Exception in preprocessWrapper for {location_str}: {str(e)}")
        # Clean up the temporary file if it exists
        if temp_file_path and os.path.exists(temp_file_path):
            try:
                os.unlink(temp_file_path)
            except OSError:
                pass
        return None


def directProcessWrapper(location_str: str, bargs: Dict[str, Any]) -> List[str]:
    """Process VCF directly without splitting into blocks.
    
    Modernized replacement for blocksplit - eliminates the overhead and complexity
    of splitting files into chunks and then merging them back together.
    
    Args:
        location_str: Location string in format chrom:start-end
        bargs: Arguments (not used in direct processing, kept for compatibility)
    
    Returns:
        List containing the original location string (no chunking)
    """
    # Direct processing - no splitting means we just return the original location
    # This eliminates the error-prone blocksplit process as requested
    return [location_str]


def partialCredit(
    vcfname: str,
    outputname: str,
    reference: str,
    locations: Optional[Union[str, List[str]]] = None,
    threads: int = 1,
    window: int = 10000,
    leftshift: bool = True,
    decompose: bool = True,
    haploid_x: bool = False,
) -> None:
    """Partial-credit-process a VCF file according to our args.

    Args:
        vcfname: Input VCF filename
        outputname: Output filename
        reference: Reference FASTA
        locations: List of regions or comma-separated string of regions
        threads: Number of threads to use
        window: Window size for blocksplit
        leftshift: Enable left-shifting of variants
        decompose: Decompose complex variants
        haploid_x: Treat X chromosome as haploid
    """
    pool = getPool(int(threads))
    if threads > 1:
        logging.info(f"Partial credit processing uses {threads} parallel processes.")

        if not locations:
            h = extractHeadersJSON(vcfname)
            if not h["tabix"]["chromosomes"]:
                logging.warning("Empty input or not tabix indexed")
                if outputname.endswith(".bcf"):
                    runBcftools("view", "-O", "b", "-o", outputname, vcfname)
                    runBcftools("index", outputname)
                else:
                    runBcftools("view", "-O", "z", "-o", outputname, vcfname)
                    runBcftools("index", "-t", outputname)
                # just return the same file
                return
            locations = h["tabix"]["chromosomes"]
        elif isinstance(locations, str):
            locations = locations.split(",")

        # use blocksplit to subdivide input
        res = runParallel(
            pool,
            directProcessWrapper,
            locations,
            {"vcf": vcfname, "dist": window, "pieces": min(40, threads * 4)},
        )

        # Direct processing - no splitting means no failures to filter
        # Flatten list of lists
        locations = [
            item for sublist in res if sublist for item in sublist
        ]
        if not locations:
            logging.warning(
                "No locations to process."
            )
            locations = [""]
    else:
        locations = [""]

    res = []
    try:
        # In Python 3, we need to create the zipped iterator explicitly
        file_loc_pairs = list(zip(itertools.repeat(vcfname), locations))

        res = runParallel(
            pool,
            preprocessWrapper,
            file_loc_pairs,
            {
                "reference": reference,
                "decompose": decompose,
                "leftshift": leftshift,
                "haploid_x": haploid_x,
                "bcf": outputname.endswith(".bcf"),
            },
        )

        # Filter out None values from results
        valid_results = [r for r in res if r is not None]

        if len(valid_results) != len(res):
            # Some preprocessing jobs failed
            logging.error(
                f"One or more preprocess jobs failed. Expected {len(res)} results, got {len(valid_results)} valid results."
            )
            raise Exception("One of the preprocess jobs failed")

        if not valid_results:
            raise Exception(
                f"No blocks were processed. List of locations: {list(locations)}"
            )

        # Update res to contain only valid results
        res = valid_results

        concatenateParts(outputname, *res)
        if outputname.endswith(".vcf.gz"):
            runBcftools("index", "-f", "-t", outputname)
        else:  # use bcf
            runBcftools("index", "-f", outputname)
    finally:
        # Clean up temporary files
        for r in res:
            if r is not None:  # Skip None values to prevent TypeError
                for suffix in ["", ".tbi", ".csi"]:
                    with contextlib.suppress(OSError):
                        os.unlink(r + suffix)

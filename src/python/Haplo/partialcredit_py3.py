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
# 3/9/2015
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

"""
Module for partial credit variant comparison processing.
Allows for preprocess, blocksplit and integration of partial credit comparisons.
"""

import itertools
import logging
import os
import shlex
import subprocess
import tempfile
import time
from typing import Any, Dict, List, Optional, Tuple, Union

from Tools.bcftools import concatenateParts, runBcftools
from Tools.parallel import getPool, runParallel
from Tools.vcfextract import extractHeadersJSON


def preprocessWrapper(file_and_location: Tuple[str, str], args: Dict[str, Any]) -> str:
    """Process a VCF file with the preprocess tool.

    Args:
        file_and_location: Tuple of (filename, location_str)
        args: Arguments for preprocessing

    Returns:
        Path to the preprocessed output file
    """
    starttime = time.time()
    filename, location_str = file_and_location
    if args["bcf"]:
        int_suffix = "bcf"
    else:
        int_suffix = "vcf.gz"

    with tempfile.NamedTemporaryFile(
        delete=False, prefix=f"input.{location_str}", suffix=f".prep.{int_suffix}"
    ) as tf:
        pass  # Just creating the file

    # Quote filenames for shell safety
    quoted_filename = shlex.quote(filename)
    quoted_location = shlex.quote(location_str) if location_str else ""
    quoted_reference = shlex.quote(args["reference"])

    location_param = f"-l {quoted_location} " if location_str else ""

    to_run = (
        f"preprocess {quoted_filename}:* {location_param}-o {tf.name} "
        f"-V {args['decompose']} -L {args['leftshift']} -r {quoted_reference}"
    )

    if args["haploid_x"]:
        to_run += " --haploid-x 1"

    with tempfile.NamedTemporaryFile(
        delete=False, prefix="stderr", suffix=".log"
    ) as tfe, tempfile.NamedTemporaryFile(
        delete=False, prefix="stdout", suffix=".log"
    ) as tfo:
        finished = False
        try:
            logging.info(f"Running '{to_run}'")
            subprocess.check_call(to_run, shell=True, stdout=tfo, stderr=tfe)
            finished = True
        finally:
            if finished:
                # Close files and read their contents for logging
                with open(tfo.name, "r", encoding="utf-8") as f:
                    for l in f:
                        logging.info(l.rstrip())
                os.unlink(tfo.name)

                with open(tfe.name, "r", encoding="utf-8") as f:
                    for l in f:
                        logging.warning(l.rstrip())
                os.unlink(tfe.name)
            else:
                logging.error(
                    f"Preprocess command {to_run} failed. "
                    f"Outputs are here {tfo.name} / {tfe.name}"
                )
                with open(tfo.name, "r", encoding="utf-8") as f:
                    for l in f:
                        logging.error(l.rstrip())
                with open(tfe.name, "r", encoding="utf-8") as f:
                    for l in f:
                        logging.error(l.rstrip())

    elapsed = time.time() - starttime
    logging.info(f"preprocess for {location_str} -- time taken {elapsed:.2f}")
    runBcftools("index", tf.name)
    return tf.name


def blocksplitWrapper(location_str: str, bargs: Dict[str, Any]) -> List[str]:
    """Blocksplit for partial credit preprocessing.

    Args:
        location_str: Location string in format chrom:start-end
        bargs: Arguments for blocksplit

    Returns:
        List of location strings for chunks
    """
    starttime = time.time()
    tf = tempfile.NamedTemporaryFile(
        delete=False, prefix=f"result.{location_str}", suffix=".chunks.bed"
    )
    result = None
    try:
        tf.close()

        # Quote strings for shell safety
        quoted_vcf = shlex.quote(bargs["vcf"])
        quoted_location = shlex.quote(location_str)

        to_run = (
            f"blocksplit {quoted_vcf} -l {quoted_location} -o {tf.name} "
            f"--window {bargs['dist']} --nblocks {bargs['pieces']} -f 0"
        )

        with tempfile.NamedTemporaryFile(
            delete=False, prefix="stderr", suffix=".log"
        ) as tfe, tempfile.NamedTemporaryFile(
            delete=False, prefix="stdout", suffix=".log"
        ) as tfo:
            try:
                logging.info(f"Running '{to_run}'")
                subprocess.check_call(to_run, shell=True, stdout=tfo, stderr=tfe)
            finally:
                # Close files and read their contents for logging
                with open(tfo.name, "r", encoding="utf-8") as f:
                    for l in f:
                        logging.info(l.rstrip())
                os.unlink(tfo.name)

                with open(tfe.name, "r", encoding="utf-8") as f:
                    for l in f:
                        logging.warning(l.rstrip())
                os.unlink(tfe.name)

        r = []
        with open(tf.name, "r", encoding="utf-8") as f:
            for l in f:
                ll = l.strip().split("\t", 3)
                if len(ll) < 3:
                    continue
                xchr = ll[0]
                start = int(ll[1]) + 1
                end = int(ll[2])
                r.append(f"{xchr}:{start}-{end}")
        result = r

    finally:
        elapsed = time.time() - starttime
        logging.info(f"blocksplit for {location_str} -- time taken {elapsed:.2f}")
        os.unlink(tf.name)

    return result


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
            blocksplitWrapper,
            locations,
            {"vcf": vcfname, "dist": window, "pieces": min(40, threads * 4)},
        )

        if None in res:
            raise Exception("One of the blocksplit processes failed.")

        # Flatten list of lists
        locations = [item for sublist in res if sublist for item in sublist]
        if not locations:
            logging.warning(
                "Blocksplit returned no blocks. This can happen when "
                "an input contains no valid variants."
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

        if None in res:
            raise Exception("One of the preprocess jobs failed")
        if not res:
            raise Exception(
                f"No blocks were processed. List of locations: {list(locations)}"
            )

        concatenateParts(outputname, *res)
        if outputname.endswith(".vcf.gz"):
            runBcftools("index", "-f", "-t", outputname)
        else:  # use bcf
            runBcftools("index", "-f", outputname)
    finally:
        # Clean up temporary files
        for r in res:
            for suffix in ["", ".tbi", ".csi"]:
                try:
                    os.unlink(r + suffix)
                except OSError:
                    pass

#!/usr/bin/env python33
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
# Diploid ROC Computation
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import contextlib
import logging
import os
import pipes
import subprocess
import tempfile
from typing import Dict, List


def fastaContigLengths(fastafile: str) -> Dict[str, int]:
    """Return contig lengths in a fasta file

    Args:
        fastafile: Path to the FASTA file

    Returns:
        Dictionary mapping contig names to lengths

    Raises:
        Exception: If the FASTA file is not indexed
    """
    if not os.path.exists(fastafile + ".fai"):
        raise Exception(f"Fasta file {fastafile} is not indexed")

    fastacontiglengths = {}

    with open(fastafile + ".fai") as fai:
        for l in fai:
            row = l.strip().split("\t")
            fastacontiglengths[row[0]] = int(row[1])

    return fastacontiglengths


def fastaNonNContigLengths(fastafile: str) -> Dict[str, int]:
    """Return contig lengths in a fasta file excluding
    N bases.

    Args:
        fastafile: Path to the FASTA file

    Returns:
        Dictionary mapping contig names to non-N lengths
    """
    # FIXME -- this could be made more efficient by using Python code
    #          instead of calling out
    #
    # NOTE: This code uses a subprocess call to 'grep' which counts the number of
    # non-N characters in the FASTA file for each contig.
    fd, t = tempfile.mkstemp(prefix="fasta_tmp")
    os.close(fd)
    try:
        cmd_line = "cat {} | grep -v '>' | tr -cd 'ACGTacgt' | wc -c > {}".format(
            pipes.quote(fastafile),
            pipes.quote(t),
        )
        logging.info(cmd_line)

        po = subprocess.Popen(
            cmd_line,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

        stdout, stderr = po.communicate()

        po.wait()

        return_code = po.returncode

        if return_code != 0:
            logging.error("cat | grep | tr | wc error: %s" % stderr)
            raise Exception("Failed to count non-N bases in %s" % fastafile)

        v = int(open(t).read().strip())
        result = {"all": v}

        # also figure contig-by-contig
        cts = fastaContigLengths(fastafile)
        for c in cts:
            cmd_line = f"samtools faidx {pipes.quote(fastafile)} {pipes.quote(c)} | grep -v '>' | tr -cd 'ACGTacgt' | wc -c"
            logging.debug(cmd_line)

            po = subprocess.Popen(
                cmd_line,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )

            stdout, stderr = po.communicate()

            po.wait()

            return_code = po.returncode

            if return_code != 0:
                logging.error("samtools faidx | grep | tr | wc error: %s" % stderr)
                raise Exception(f"Failed to count non-N bases in {fastafile}:{c}")

            result[c] = int(stdout.strip())
    finally:
        with contextlib.suppress(Exception):
            os.unlink(t)

    return result


def fastaSampleRegions(
    fastafile: str, n_regions: int = 10, region_length: int = 10000
) -> List[str]:
    """Sample regions from a fasta file

    Args:
        fastafile: Path to the FASTA file
        n_regions: Number of regions to sample
        region_length: Length of each region

    Returns:
        List of sampled regions in format "chrom:start-end"
    """
    import random

    result = []
    cts = fastaContigLengths(fastafile)

    if n_regions == 0:
        return result
    elif n_regions == 1:
        l = list(cts.keys())
        c = l[random.randint(0, len(l) - 1)]
        start = random.randint(0, max(0, cts[c] - region_length))
        result = ["%s:%i-%i" % (c, start, start + region_length)]
        return result

    total = 0
    for c in cts:
        total += cts[c]

    active_regions = []

    for c in sorted(cts.keys()):
        ar_count = int(
            1.5 + (n_regions - len(active_regions)) * cts[c] / max(1.0, total)
        )
        logging.info("Adding %i regions from %s of %i bp", ar_count, c, cts[c])

        if ar_count > 0 and 11 * region_length < cts[c]:
            for _x in range(ar_count):
                # make sure we leave 10*region_length bp between active regions, and also try
                # reasonable hard to not have +-10*region_length overlap with any existing
                # active region
                position_ok = False
                max_tries = 20
                tries = 0

                start = random.randint(0, max(0, cts[c] - region_length))
                while not position_ok and tries < max_tries:
                    position_ok = True

                    for ar in active_regions:
                        ar_chr = ar.split(":")[0]
                        ar_start = int(ar.split(":")[1].split("-")[0])
                        ar_end = int(ar.split(":")[1].split("-")[1])

                        if ar_chr == c:
                            # require 10*region_length distance to AR boundaries
                            if (
                                start < ar_start
                                and start + 11 * region_length > ar_start
                            ) or (
                                start > ar_start - 11 * region_length
                                and start < ar_end + 11 * region_length
                            ):
                                position_ok = False
                                tries += 1
                                start = random.randint(
                                    0, max(0, cts[c] - region_length)
                                )
                                break

                if position_ok:
                    ar = "%s:%i-%i" % (c, start, start + region_length)
                    active_regions.append(ar)

    return active_regions

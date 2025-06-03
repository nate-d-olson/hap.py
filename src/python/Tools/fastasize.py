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

import ast
import logging
import os
import re
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

    with open(fastafile + ".fai", encoding="utf-8") as fai:
        for l in fai:
            row = l.strip().split("\t")
            fastacontiglengths[row[0]] = int(row[1])

    return fastacontiglengths


def calculateLength(fastacontiglengths, locations) -> int:
    """Calculate total base count for given regions or contigs.

    fastacontiglengths: dict or str; mapping contig name to length
    locations: str or iterable of 'contig' or 'contig:start-end' specs
    """
    # Parse lengths dict if given as string
    if isinstance(fastacontiglengths, str):
        try:
            fastacontiglengths = ast.literal_eval(fastacontiglengths)
        except Exception:
            raise ValueError("Invalid contig lengths specification")
    if not isinstance(fastacontiglengths, dict):
        raise ValueError("fastacontiglengths must be a dict or string repr")
    # Build list of location specs
    if isinstance(locations, str):
        # split on whitespace or commas
        locs = [tok for tok in re.split(r"[\s,]+", locations) if tok]
    else:
        locs = list(locations)
    total = 0
    for spec in locs:
        if ":" in spec and "-" in spec:
            name, rng = spec.split(":", 1)
            parts = rng.split("-", 1)
            try:
                start = int(parts[0])
                end = int(parts[1])
            except ValueError:
                continue
            length = max(0, end - start + 1)
        else:
            name = spec
            length = fastacontiglengths.get(name, 0)
        total += length
    return total


def fastaNonNContigLengths(fastafile: str) -> Dict[str, int]:
    """Return contig lengths in a FASTA file excluding N bases.

    Args:
        fastafile: Path to the FASTA file

    Returns:
        Dictionary mapping 'all' and each contig to non-N base counts
    """
    counts: Dict[str, int] = {}
    current = None
    try:
        with open(fastafile, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].strip().split()[0]
                    counts[name] = 0
                    current = name
                else:
                    if current is None:
                        continue
                    seq = line.strip().upper()
                    counts[current] += sum(1 for c in seq if c in ("A", "C", "G", "T"))
    except Exception as e:
        raise Exception(f"Failed to read FASTA file {fastafile}: {e}")
    total = sum(counts.values())
    result: Dict[str, int] = {"all": total}
    result.update(counts)
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

    if n_regions <= 0:
        return []
    if n_regions == 1:
        chromosomes = list(cts.keys())
        chosen = random.choice(chromosomes)
        max_start = max(0, cts[chosen] - region_length)
        start = random.randint(0, max_start)
        return [f"{chosen}:{start}-{start + region_length}"]

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

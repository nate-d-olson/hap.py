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

import os
import logging
import subprocess
import tempfile
import time
import itertools
import multiprocessing
import pipes

from Tools.parallel import runParallel, getPool
from Tools.bcftools import runBcftools, concatenateParts
from Tools.vcfextract import extractHeadersJSON


# Python 3 compatibility for file handling
def open_file(filename, mode="r"):
    """Helper function to open files in the correct mode for both text and binary."""
    if "b" in mode:
        return open(filename, mode)
    else:
        return open(filename, mode, encoding="utf-8")


def preprocessWrapper(file_and_location, args):
    starttime = time.time()
    filename, location_str = file_and_location
    if args["bcf"]:
        int_suffix = "bcf"
    else:
        int_suffix = "vcf.gz"

    tf = tempfile.NamedTemporaryFile(
        delete=False, prefix=f"input.{location_str}", suffix=f".prep.{int_suffix}"
    )
    tf.close()

    to_run = f"preprocess {pipes.quote(filename)}:* {pipes.quote(location_str)} -o {tf.name} -V {args['decompose']} -L {args['leftshift']} -r {pipes.quote(args['reference'])}"
    if args["haploid_x"]:
        to_run += " --haploid-x 1"

    tfe = tempfile.NamedTemporaryFile(delete=False, prefix="stderr", suffix=".log")
    tfo = tempfile.NamedTemporaryFile(delete=False, prefix="stdout", suffix=".log")
    finished = False
    try:
        logging.info(f"Running '{to_run}'")
        subprocess.check_call(to_run, shell=True, stdout=tfo, stderr=tfe)
        finished = True
    finally:
        if finished:
            tfo.close()
            tfe.close()
            with open_file(tfo.name) as f:
                for l in f:
                    logging.info(l.replace("\n", ""))
            os.unlink(tfo.name)
            with open_file(tfe.name) as f:
                for l in f:
                    logging.warn(l.replace("\n", ""))
            os.unlink(tfe.name)
        else:
            logging.error(
                f"Preprocess command {to_run} failed. Outputs are here {tfo.name} / {tfe.name}"
            )
            with open_file(tfo.name) as f:
                for l in f:
                    logging.error(l.replace("\n", ""))
            with open_file(tfe.name) as f:
                for l in f:
                    logging.error(l.replace("\n", ""))

    elapsed = time.time() - starttime
    logging.info(f"preprocess for {location_str} -- time taken {elapsed:.2f}")
    runBcftools("index", tf.name)
    return tf.name


def blocksplitWrapper(location_str, bargs):
    """Blocksplit for partial credit preprocessing"""
    starttime = time.time()
    tf = tempfile.NamedTemporaryFile(
        delete=False, prefix=f"result.{location_str}", suffix=".chunks.bed"
    )
    result = None
    try:
        tf.close()

        to_run = f"blocksplit {pipes.quote(bargs['vcf'])} -l {pipes.quote(location_str)} -o {tf.name} --window {bargs['dist']} --nblocks {bargs['pieces']} -f 0"

        tfe = tempfile.NamedTemporaryFile(delete=False, prefix="stderr", suffix=".log")
        tfo = tempfile.NamedTemporaryFile(delete=False, prefix="stdout", suffix=".log")
        try:
            logging.info(f"Running '{to_run}'")
            subprocess.check_call(to_run, shell=True, stdout=tfo, stderr=tfe)
        finally:
            tfo.close()
            tfe.close()
            with open_file(tfo.name) as f:
                for l in f:
                    logging.info(l.replace("\n", ""))
            os.unlink(tfo.name)
            with open_file(tfe.name) as f:
                for l in f:
                    logging.warn(l.replace("\n", ""))
            os.unlink(tfe.name)

        r = []
        with open_file(tf.name) as f:
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
    vcfname,
    outputname,
    reference,
    locations,
    threads=1,
    window=10000,
    leftshift=True,
    decompose=True,
    haploid_x=False,
):
    """Partial-credit-process a VCF file according to our args"""

    pool = getPool(int(threads))
    if threads > 1:
        logging.info("Partial credit processing uses %i parallel processes." % threads)

        if not locations:
            h = extractHeadersJSON(vcfname)
            if not h["tabix"]["chromosomes"]:
                logging.warn("Empty input or not tabix indexed")
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

        locations = list(itertools.chain.from_iterable(res))
        if not len(locations):
            logging.warn(
                "Blocksplit returned no blocks. This can happen when "
                "an input contains no valid variants."
            )
            locations = [""]
    else:
        locations = [""]

    res = []
    try:
        res = runParallel(
            pool,
            preprocessWrapper,
            list(zip(itertools.repeat(vcfname), locations)),
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
        for r in res:
            try:
                os.unlink(r)
            except:
                pass
            try:
                os.unlink(r + ".tbi")
            except:
                pass
            try:
                os.unlink(r + ".csi")
            except:
                pass

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
# Diploid VCF File Comparison
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import os
import logging
import tempfile
import time
import subprocess

# Python 3 compatibility for file handling
def open_file(filename, mode='r'):
    """Helper function to open files in the correct mode for both text and binary."""
    if 'b' in mode:
        return open(filename, mode)
    else:
        return open(filename, mode, encoding='utf-8')

def xcmpWrapper(location_str, args):
    """ Haplotype block comparison wrapper function
    """
    starttime = time.time()
    tf = tempfile.NamedTemporaryFile(delete=False,
                                     dir=args.scratch_prefix,
                                     prefix="result.%s" % location_str,
                                     suffix=".bcf")
    tf.close()

    to_run = "xcmp %s %s -l %s -o %s -r %s -f %i -n %i --expand-hapblocks %i " \
             "--window %i --no-hapcmp %i --qq %s" % \
             (args.vcf1.replace(" ", "\\ "),
              args.vcf2.replace(" ", "\\ "),
              location_str,
              tf.name,
              args.ref,
              1 if args.pass_only else 0,  # -f == apply-filtering
              args.max_enum,
              args.hb_expand,
              args.window,
              1 if args.no_hc else 0,
              args.roc if args.roc else "QUAL")

    if args.verbose:
        # this prints information on failed sites
        to_run += " -e -"

    # regions / targets already have been taken care of in blocksplit / preprocessing

    tfe = tempfile.NamedTemporaryFile(delete=False,
                                      dir=args.scratch_prefix,
                                      prefix="stderr",
                                      suffix=".log")
    tfo = tempfile.NamedTemporaryFile(delete=False,
                                      dir=args.scratch_prefix,
                                      prefix="stdout",
                                      suffix=".log")

    try:
        logging.info("Running '%s'" % to_run)
        process = subprocess.Popen(to_run, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        # Handle bytes vs. string in Python 3
        if isinstance(stdout, bytes):
            stdout = stdout.decode('utf-8')
        if isinstance(stderr, bytes):
            stderr = stderr.decode('utf-8')
            
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, to_run, output=stdout, stderr=stderr)
            
        with open_file(tfo.name, 'w') as f:
            f.write(stdout)
        with open_file(tfe.name, 'w') as f:
            f.write(stderr)
    finally:
        with open_file(tfo.name) as f:
            for l in f:
                logging.info(l.replace("\n", ""))
        os.unlink(tfo.name)
        with open_file(tfe.name) as f:
            for l in f:
                logging.warn(l.replace("\n", ""))
        os.unlink(tfe.name)

    elapsed = time.time() - starttime
    logging.info("xcmp for chunk %s -- time taken %.2f" % (location_str, elapsed))

    return tf.name

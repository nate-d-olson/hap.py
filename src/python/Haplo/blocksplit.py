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
import pipes

# Python 3 compatibility for file handling
def open_file(filename, mode='r'):
    """Helper function to open files in the correct mode for both text and binary."""
    if 'b' in mode:
        return open(filename, mode)
    else:
        return open(filename, mode, encoding='utf-8')

def blocksplitWrapper(location_str, args):
    starttime = time.time()
    tf = tempfile.NamedTemporaryFile(delete=False,
                                     dir=args.scratch_prefix,
                                     prefix="result.%s" % location_str,
                                     suffix=".chunks.bed")
    tf.close()

    if location_str:
        loc = " -l %s" % pipes.quote(location_str)
    else:
        loc = ""
    to_run = "blocksplit %s %s%s -o %s --window %i --nblocks %i -f 0" % \
             (pipes.quote(args.vcf1),
              pipes.quote(args.vcf2),
              loc,
              tf.name,
              args.window*2,
              args.pieces)

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
    logging.info("blocksplit for %s -- time taken %.2f" % (location_str, elapsed))
    return tf.name

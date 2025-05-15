#!/usr/bin/env python3
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
import shlex
from typing import List, Dict, Tuple, Optional, Union, Any


def blocksplitWrapper(location_str: Optional[str], args: Any) -> Tuple[str, float]:
    """Run the blocksplit utility to split regions into chunks
    
    Args:
        location_str: Optional region string to restrict analysis to
        args: Arguments with vcf1, vcf2, scratch_prefix, window, pieces fields
    
    Returns:
        Tuple of (output bed file path, execution time)
    """
    starttime = time.time()
    
    # Create a temporary file for results
    with tempfile.NamedTemporaryFile(delete=False,
                                   dir=args.scratch_prefix,
                                   prefix=f"result.{location_str or 'all'}" ,
                                   suffix=".chunks.bed") as tf:
        output_filename = tf.name

    # Build location argument if provided
    loc_arg = f" -l {shlex.quote(location_str)}" if location_str else ""
    
    # Construct the command
    to_run = f"blocksplit {shlex.quote(args.vcf1)} {shlex.quote(args.vcf2)}{loc_arg} " \
             f"-o {output_filename} --window {args.window*2} --nblocks {args.pieces} -f 0"

    # Create temporary files for stdout/stderr
    with tempfile.NamedTemporaryFile(delete=False,
                                    dir=args.scratch_prefix,
                                    prefix="stderr",
                                    suffix=".log") as tfe, \
         tempfile.NamedTemporaryFile(delete=False,
                                    dir=args.scratch_prefix,
                                    prefix="stdout",
                                    suffix=".log") as tfo:
        stderr_file = tfe.name
        stdout_file = tfo.name
    
    logging.info(f"Running: {to_run}")
    
    try:
        with open(stdout_file, "w") as out_fh, open(stderr_file, "w") as err_fh:
            p = subprocess.Popen(to_run, shell=True,
                               stdout=out_fh,
                               stderr=err_fh)
            p.wait()
            
            # Check if the process failed
            if p.returncode != 0:
                with open(stderr_file, "r") as err_fh:
                    error_output = err_fh.read()
                logging.error(f"Command failed with return code {p.returncode}\nError output: {error_output}")
                raise Exception(f"Command failed: {to_run}")
    except Exception as e:
        logging.error(f"Failed to run blocksplit: {e}")
        raise
    finally:
        # Clean up temporary files
        for f in [stderr_file, stdout_file]:
            try:
                if os.path.exists(f):
                    os.unlink(f)
            except Exception as e:
                logging.warning(f"Failed to remove temporary file {f}: {e}")
    
    execution_time = time.time() - starttime
    logging.info(f"blocksplit completed in {execution_time:.2f} seconds")
    
    return output_filename, execution_time

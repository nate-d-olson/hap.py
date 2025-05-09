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
# 10/02/2017
#
# Run gvcf2bed

import os
import tempfile
import subprocess
import json
import logging
import Tools


# Python 3 compatibility for file handling
def open_file(filename, mode="r"):
    """Helper function to open files in the correct mode for both text and binary."""
    if "b" in mode:
        return open(filename, mode)
    else:
        return open(filename, mode, encoding="utf-8")


def gvcf2bed(vcf, ref, regions=None, scratch_prefix=None):
    """Run gvcf2bed and return temporary region bed file in temp folder"""

    tf = tempfile.NamedTemporaryFile(dir=scratch_prefix, suffix=".bed")
    tf.close()
    cmdline = "gvcf2bed %s -r %s -o %s" % (pipes.quote(vcf), pipes.quote(ref), tf.name)
    if regions:
        cmdline += " -T %s" % pipes.quote(regions)
    logging.info("Running gvcf2bed: '%s'" % cmdline)

    # Use Popen instead of check_call for better bytes handling in Python 3
    process = subprocess.Popen(
        cmdline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()

    # Handle bytes vs. string in Python 3
    if isinstance(stdout, bytes):
        stdout = stdout.decode("utf-8")
    if isinstance(stderr, bytes):
        stderr = stderr.decode("utf-8")

    if process.returncode != 0:
        raise subprocess.CalledProcessError(
            process.returncode, cmdline, output=stdout, stderr=stderr
        )

    if stdout.strip() or stderr.strip():
        logging.info("gvcf2bed output: \n%s\n / \n%s\n" % (stdout, stderr))

    return tf.name

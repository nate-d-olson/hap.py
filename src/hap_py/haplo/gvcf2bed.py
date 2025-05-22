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

import logging
import shlex
import subprocess
import tempfile


def gvcf2bed(vcf, ref, regions=None, scratch_prefix=None):
    """Run gvcf2bed and return temporary region bed file in temp folder"""

    tf = tempfile.NamedTemporaryFile(dir=scratch_prefix, suffix=".bed")
    tf.close()
    cmdline = f"gvcf2bed {shlex.quote(vcf)} -r {shlex.quote(ref)} -o {tf.name}"
    if regions:
        cmdline += " -T %s" % shlex.quote(regions)
    logging.info("Running gvcf2bed: '%s'" % cmdline)
    subprocess.check_call(shlex.split(cmdline))
    return tf.name

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
# 3/9/2014
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import logging
import os
import shutil
import subprocess
import tempfile
import time

import Haplo.version  # pylint: disable=E0611,E0401


# Python 3 compatibility for file handling
def open_file(filename, mode="r"):
    """Helper function to open files in the correct mode for both text and binary."""
    if "b" in mode:
        return open(filename, mode)
    else:
        return open(filename, mode, encoding="utf-8")


def findVCFEval():
    """Return default version of rtgtools if hap.py was built with
    rtgtools included.
    """
    if Haplo.version.has_vcfeval:
        base = os.path.join(
            os.path.dirname(__file__),
            "..",  # python27
            "..",  # lib
            "..",  # hap.py-base
            "libexec",
            "rtg-tools-install",
        )
        # prefer wrapper when it's there
        bfile = os.path.join(base, "rtg-wrapper.sh")
        bfile2 = os.path.join(base, "rtg")
        if os.path.isfile(bfile) and os.access(bfile, os.X_OK):
            return bfile
        elif os.path.isfile(bfile2) and os.access(bfile2, os.X_OK):
            return bfile2
        else:
            logging.warn(
                f"Could not find our included version of rtg-tools at {base}. "
                + "To use vcfeval for comparison, you might have to specify "
                "its location on the command line."
            )
            return "rtg"
    else:
        # default: return
        return "rtg"


def runVCFEval(vcf1, vcf2, target, args):
    """Run VCFEval and convert it's output to something quantify
    understands
    """
    starttime = time.time()

    vtf = tempfile.NamedTemporaryFile(
        dir=args.scratch_prefix, prefix="vcfeval.result", suffix=".dir"
    )
    vtf.close()

    del_sdf = False

    try:
        if not args.engine_vcfeval_template and os.path.isdir(args.ref[:-3] + ".sdf"):
            logging.info(f"Using vcfeval template from {args.ref[:-3]}.sdf")
            args.engine_vcfeval_template = args.ref[:-3] + ".sdf"
        if not args.engine_vcfeval_template or not os.path.exists(
            args.engine_vcfeval_template
        ):
            logging.warn(
                f"Creating template for vcfeval. You can speed this up by supplying a SDF template that corresponds to {args.ref}"
            )
            del_sdf = True
            stf = tempfile.NamedTemporaryFile(
                dir=args.scratch_prefix, prefix="vcfeval.sdf", suffix=".dir"
            )
            stf.close()
            args.engine_vcfeval_template = stf.name
            runme = f"{findVCFEval()} format -o {args.ref[:-3]}.sdf {args.ref}"
            logging.info(runme)
            po = subprocess.Popen(
                runme, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = po.communicate()

            po.wait()
            return_code = po.returncode
            if return_code != 0:
                raise Exception(
                    f"Error running rtg tools. Return code was {return_code}, output: {stdout} / {stderr} \n"
                )
            elif stdout.strip() or stderr.strip():
                logging.info(f"RTG output: \n{stdout}\n / \n{stderr}\n")

        # construct vcfeval command
        runme = (
            f"{findVCFEval()} vcfeval "
            f"-b {pipes.quote(vcf1)} -c {pipes.quote(vcf2)} "
            f"-t {pipes.quote(args.engine_vcfeval_template)} "
            f"-o {pipes.quote(vtf.name)} "
            f"-T {args.threads} -m ga4gh --ref-overlap"
        )

        if not args.pass_only:
            runme += " --all-records"

        if args.roc:
            runme += " -f %s" % pipes.quote(args.roc)

        logging.info(runme)
        po = subprocess.Popen(
            runme, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        stdout, stderr = po.communicate()

        po.wait()
        return_code = po.returncode
        if return_code != 0:
            raise Exception(
                f"Error running rtg tools / vcfeval. Return code was {return_code}, output: {stdout} / {stderr} \n"
            )
        elif stdout.strip() or stderr.strip():
            logging.info(f"vcfeval output: \n{stdout}\n / \n{stderr}\n")

        # in GA4GH mode, this is what vcfeval should output
        shutil.copy(os.path.join(vtf.name, "output.vcf.gz"), target)
        shutil.copy(os.path.join(vtf.name, "output.vcf.gz.tbi"), target + ".tbi")
    finally:
        # remove temp path
        try:
            shutil.rmtree(vtf.name)
        except Exception as e:
            pass
        if del_sdf:
            try:
                shutil.rmtree(args.engine_vcfeval_template)
            except Exception as e:
                pass

    elapsed = time.time() - starttime
    logging.info("vcfeval for %s vs. %s -- time taken %.2f" % (vcf1, vcf2, elapsed))

    if os.path.exists(target) and os.path.exists(target + ".tbi"):
        return [target, target + ".tbi"]
    else:
        return None

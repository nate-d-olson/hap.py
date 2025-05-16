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

"""
Module for running RTG's vcfeval for variant comparison.
Provides functionality to compare VCF files using the vcfeval tool.
"""

import contextlib
import logging
import os
import shlex
import shutil
import subprocess
import tempfile
import time
from typing import Any, List, Optional

import Haplo.version  # pylint: disable=E0611,E0401


def findVCFEval() -> str:
    """Return default version of rtgtools if hap.py was built with rtgtools included.

    Returns:
        Path to rtg executable or 'rtg' if not found
    """
    if Haplo.version.has_vcfeval:
        base = os.path.join(
            os.path.dirname(__file__),
            "..",  # python
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
            logging.warning(
                f"Could not find our included version of rtg-tools at {base}. "
                "To use vcfeval for comparison, you might have to specify "
                "its location on the command line."
            )
            return "rtg"
    else:
        # default: return
        return "rtg"


def runVCFEval(vcf1: str, vcf2: str, target: str, args: Any) -> Optional[List[str]]:
    """Run VCFEval and convert its output to something quantify understands.

    Args:
        vcf1: First VCF file (baseline)
        vcf2: Second VCF file (query)
        target: Output file path
        args: Command line arguments

    Returns:
        List of output files or None if failed
    """
    starttime = time.time()

    with tempfile.NamedTemporaryFile(
        dir=args.scratch_prefix, prefix="vcfeval.result", suffix=".dir"
    ) as vtf:
        pass  # Just create the file to get the name

    del_sdf = False

    try:
        # Check for existing SDF template
        if not args.engine_vcfeval_template and os.path.isdir(args.ref[:-3] + ".sdf"):
            logging.info(f"Using vcfeval template from {args.ref[:-3] + '.sdf'}")
            args.engine_vcfeval_template = args.ref[:-3] + ".sdf"

        # Create template if needed
        if not args.engine_vcfeval_template or not os.path.exists(
            args.engine_vcfeval_template
        ):
            logging.warning(
                "Creating template for vcfeval. "
                f"You can speed this up by supplying a SDF template that corresponds to {args.ref}"
            )
            del_sdf = True
            with tempfile.NamedTemporaryFile(
                dir=args.scratch_prefix, prefix="vcfeval.sdf", suffix=".dir"
            ) as stf:
                pass  # Just create the file to get the name

            args.engine_vcfeval_template = stf.name

            # Quote paths for shell safety
            quoted_engine = shlex.quote(args.engine_vcfeval)
            quoted_template = shlex.quote(args.engine_vcfeval_template)
            quoted_ref = shlex.quote(args.ref)

            runme = f"{quoted_engine} format -o {quoted_template} {quoted_ref}"

            logging.info(runme)
            process = subprocess.Popen(
                runme,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            stdout, stderr = process.communicate()
            rc = process.returncode

            if rc != 0:
                raise Exception(
                    f"Error running rtg tools. Return code was {rc}, "
                    f"output: {stdout} / {stderr}"
                )
            elif stdout.strip() or stderr.strip():
                logging.info(f"RTG output: \n{stdout}\n / \n{stderr}\n")

        # Quote paths for shell safety
        quoted_engine = shlex.quote(args.engine_vcfeval)
        quoted_vcf1 = shlex.quote(vcf1)
        quoted_vcf2 = shlex.quote(vcf2)
        quoted_template = shlex.quote(args.engine_vcfeval_template)
        quoted_output = shlex.quote(vtf.name)

        runme = (
            f"{quoted_engine} vcfeval -b {quoted_vcf1} -c {quoted_vcf2} "
            f"-t {quoted_template} -o {quoted_output} -T {args.threads} "
            f"-m ga4gh --ref-overlap"
        )

        if not args.pass_only:
            runme += " --all-records"

        if args.roc:
            runme += f" -f {shlex.quote(args.roc)}"

        if args.engine_scmp_distance:
            runme += f" --Xloose-match-distance={args.engine_scmp_distance}"

        logging.info(runme)
        process = subprocess.Popen(
            runme,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

        stdout, stderr = process.communicate()
        rc = process.returncode

        if rc != 0:
            raise Exception(
                f"Error running rtg tools / vcfeval. Return code was {rc}, "
                f"output: {stdout} / {stderr}"
            )
        elif stdout.strip() or stderr.strip():
            logging.info(f"vcfeval output: \n{stdout}\n / \n{stderr}\n")

        # in GA4GH mode, this is what vcfeval should output
        shutil.copy(os.path.join(vtf.name, "output.vcf.gz"), target)
        shutil.copy(os.path.join(vtf.name, "output.vcf.gz.tbi"), target + ".tbi")
    finally:
        # remove temp paths
        with contextlib.suppress(OSError):
            shutil.rmtree(vtf.name)

        if del_sdf:
            with contextlib.suppress(OSError):
                shutil.rmtree(args.engine_vcfeval_template)

    elapsed = time.time() - starttime
    logging.info(f"vcfeval for {vcf1} vs. {vcf2} -- time taken {elapsed:.2f}")

    if os.path.exists(target) and os.path.exists(target + ".tbi"):
        return [target, target + ".tbi"]
    else:
        return None

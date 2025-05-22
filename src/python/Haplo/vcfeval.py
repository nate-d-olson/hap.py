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
from argparse import Namespace
from typing import Any, Dict, List, Optional, Union

# Set up versioning
try:
    from Haplo import version

    has_vcfeval = getattr(version, "has_vcfeval", False)
except ImportError:
    # Version module not available, assume vcfeval is not included
    has_vcfeval = False


def findVCFEval() -> str:
    """Return default version of rtgtools if hap.py was built with rtgtools included.

    Returns:
        Path to rtg executable or 'rtg' if not found
    """
    if has_vcfeval:
        script_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
        base = os.path.abspath(
            os.path.join(
                script_dir,  # Haplo
                "..",  # python
                "..",  # src
                "..",  # hap.py-base
                "libexec",
                "rtg-tools-install",
            )
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


def runVCFEval(
    vcf1: str, vcf2: str, target: str, args: Union[Namespace, Dict[str, Any]]
) -> Optional[List[str]]:
    """Run VCFEval and convert its output to something quantify understands.

    Args:
        vcf1: First VCF file (baseline)
        vcf2: Second VCF file (query)
        target: Output file path
        args: Command line arguments object with the following attributes:
            - scratch_prefix: Directory for temporary files
            - engine_vcfeval_template: Optional path to an existing SDF template
            - ref: Path to reference FASTA file
            - engine_vcfeval: Path to rtg executable
            - threads: Number of threads to use
            - pass_only: Whether to use only PASS variants
            - roc: ROC file attribute (optional)
            - engine_scmp_distance: Distance for loose matching (optional)

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
        # Validate input files exist
        if not os.path.exists(vcf1):
            raise FileNotFoundError(f"Baseline VCF file not found: {vcf1}")
        if not os.path.exists(vcf2):
            raise FileNotFoundError(f"Query VCF file not found: {vcf2}")
        if not os.path.exists(args.ref):
            raise FileNotFoundError(f"Reference file not found: {args.ref}")

        # Ensure required attributes are present
        if not hasattr(args, "engine_vcfeval") or not args.engine_vcfeval:
            logging.warning("engine_vcfeval not specified, using default findVCFEval()")
            args.engine_vcfeval = findVCFEval()

        if not hasattr(args, "engine_vcfeval_template"):
            args.engine_vcfeval_template = None

        # Set default for scratch_prefix if not provided
        if not hasattr(args, "scratch_prefix") or not args.scratch_prefix:
            args.scratch_prefix = tempfile.gettempdir()

        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(target)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

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
            template_dir = None

            try:
                with tempfile.NamedTemporaryFile(
                    dir=args.scratch_prefix, prefix="vcfeval.sdf", suffix=".dir"
                ) as stf:
                    template_dir = stf.name

                # Ensure template dir exists
                if not os.path.exists(template_dir):
                    os.makedirs(template_dir, exist_ok=True)

                args.engine_vcfeval_template = template_dir

                # Quote paths for shell safety
                quoted_engine = shlex.quote(args.engine_vcfeval)
                quoted_template = shlex.quote(args.engine_vcfeval_template)
                quoted_ref = shlex.quote(args.ref)

                runme = f"{quoted_engine} format -o {quoted_template} {quoted_ref}"

                logging.info(f"Creating SDF template with command: {runme}")
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
                    error_msg = (
                        f"Error running rtg tools. Return code was {rc}, "
                        f"output: {stdout} / {stderr}"
                    )
                    logging.error(error_msg)
                    raise subprocess.SubprocessError(error_msg)
                elif stdout.strip() or stderr.strip():
                    logging.info(f"RTG output: \n{stdout}\n / \n{stderr}\n")
            except Exception as e:
                logging.error(f"Failed to create SDF template: {str(e)}")
                if template_dir and os.path.exists(template_dir):
                    with contextlib.suppress(OSError):
                        shutil.rmtree(template_dir)
                raise

        # Check for threads parameter
        if not hasattr(args, "threads") or not args.threads:
            args.threads = 1

        # Check for pass_only parameter
        if not hasattr(args, "pass_only"):
            args.pass_only = False

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

        if hasattr(args, "roc") and args.roc:
            runme += f" -f {shlex.quote(args.roc)}"

        if hasattr(args, "engine_scmp_distance") and args.engine_scmp_distance:
            runme += f" --Xloose-match-distance={args.engine_scmp_distance}"

        logging.info(f"Running vcfeval with command: {runme}")

        # Run vcfeval command
        try:
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
                error_msg = (
                    f"Error running rtg tools / vcfeval. Return code was {rc}, "
                    f"output: {stdout} / {stderr}"
                )
                logging.error(error_msg)
                raise subprocess.SubprocessError(error_msg)
            elif stdout.strip() or stderr.strip():
                logging.info(f"vcfeval output: \n{stdout}\n / \n{stderr}\n")
        except Exception as e:
            logging.error(f"Failed to run vcfeval: {str(e)}")
            raise

        # Copy output files to target location
        output_vcf = os.path.join(vtf.name, "output.vcf.gz")
        output_tbi = os.path.join(vtf.name, "output.vcf.gz.tbi")

        if not os.path.exists(output_vcf):
            logging.error(f"vcfeval did not produce expected output file: {output_vcf}")
            return None

        if not os.path.exists(output_tbi):
            logging.error(f"vcfeval did not produce expected index file: {output_tbi}")
            return None

        try:
            # in GA4GH mode, this is what vcfeval should output
            shutil.copy(output_vcf, target)
            shutil.copy(output_tbi, target + ".tbi")
        except Exception as e:
            logging.error(f"Failed to copy vcfeval output files: {str(e)}")
            return None
    finally:
        # remove temp paths
        try:
            if vtf.name and os.path.exists(vtf.name):
                logging.debug(f"Cleaning up temporary directory: {vtf.name}")
                shutil.rmtree(vtf.name)
        except OSError as e:
            logging.warning(
                f"Failed to remove temporary directory {vtf.name}: {str(e)}"
            )

        if (
            del_sdf
            and hasattr(args, "engine_vcfeval_template")
            and args.engine_vcfeval_template
        ):
            try:
                if os.path.exists(args.engine_vcfeval_template):
                    logging.debug(
                        f"Cleaning up SDF template directory: {args.engine_vcfeval_template}"
                    )
                    shutil.rmtree(args.engine_vcfeval_template)
            except OSError as e:
                logging.warning(
                    f"Failed to remove SDF template directory {args.engine_vcfeval_template}: {str(e)}"
                )

    elapsed = time.time() - starttime
    logging.info(f"vcfeval for {vcf1} vs. {vcf2} -- time taken {elapsed:.2f}")

    # Final verification that output files exist
    if os.path.exists(target) and os.path.exists(target + ".tbi"):
        logging.info(
            f"vcfeval completed successfully, output files: {target}, {target + '.tbi'}"
        )
        return [target, target + ".tbi"]
    else:
        logging.error("vcfeval failed to produce expected output files")
        return None

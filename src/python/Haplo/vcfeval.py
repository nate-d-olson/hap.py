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
from typing import Any, List, Optional, Tuple

# Persistent template cache directory (e.g. ~/.cache/happy/sdf)
_CACHE_DIR = os.environ.get(
    "HAPPY_CACHE_DIR", os.path.join(os.path.expanduser("~"), ".cache", "happy")
)

# Environment variable override for the vcfeval executable
_VCFEVAL_ENV = "HAPPY_VCFEVAL"


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------


def _get_cached_template(ref_fasta: str) -> Tuple[Optional[str], bool]:
    """Return path to cached SDF template for *ref_fasta*.

    Returns (template_path, exists) where *exists* indicates whether the path
    already existed on disk.  Callers are responsible for creating the
    template when *exists* is ``False``.
    """

    fasta_name = os.path.basename(ref_fasta)
    name_no_ext = os.path.splitext(os.path.splitext(fasta_name)[0])[0]
    tmpl_dir = os.path.join(_CACHE_DIR, "sdf", name_no_ext + ".sdf")
    return tmpl_dir, os.path.isdir(tmpl_dir)


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
        # env-override first
        if os.getenv(_VCFEVAL_ENV):
            return os.getenv(_VCFEVAL_ENV)  # type: ignore[return-value]

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

    # Flag only used in legacy cleanup path – keep for compatibility
    del_sdf = False  # noqa: F841

    try:
        # Resolve SDF template
        if args.engine_vcfeval_template:
            logging.info(
                "Using user-provided vcfeval template at %s",
                args.engine_vcfeval_template,
            )
        else:
            # 1) Env cache ~/.cache/happy/sdf/<ref>.sdf
            tmpl_path, tmpl_exists = _get_cached_template(args.ref)
            if tmpl_exists:
                logging.info("Using cached vcfeval template at %s", tmpl_path)
                args.engine_vcfeval_template = tmpl_path
            # 2) Local <ref>.sdf sibling directory
            elif os.path.isdir(args.ref[:-3] + ".sdf"):
                args.engine_vcfeval_template = args.ref[:-3] + ".sdf"
                logging.info(
                    "Using sibling vcfeval template at %s", args.engine_vcfeval_template
                )

        # Create template if needed
        if not args.engine_vcfeval_template or not os.path.exists(
            args.engine_vcfeval_template
        ):
            logging.warning(
                "Creating template for vcfeval. "
                f"You can speed this up by supplying a SDF template that corresponds to {args.ref}"
            )
            # No template available – build one in cache dir
            # We persist templates in cache now -> no cleanup required
            tmpl_path, _ = _get_cached_template(args.ref)
            # dirname is always defined here
            os.makedirs(os.path.dirname(tmpl_path), exist_ok=True)  # type: ignore[arg-type]
            args.engine_vcfeval_template = tmpl_path

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

        # Add ROC feature selection only when ROC computation is enabled.  The
        # CLI exposes a pair of mutually–exclusive flags ``--roc`` / ``--no-roc``
        # via ``happy.qfy.updateArgs`` which set ``args.roc`` (feature name) and
        # ``args.do_roc`` (boolean).  Historically we *always* passed ``-f`` to
        # vcfeval because ``args.roc`` had a default of "QUAL".  This prevented
        # users from disabling ROC generation and incurred unnecessary work in
        # vcfeval.  We now respect ``--no-roc`` by only appending the flag when
        # ``args.do_roc`` is truthy.

        if getattr(args, "do_roc", True):
            # When ROC computation is requested, ``args.roc`` contains the INFO
            # field (or QUAL, GQX, …) to use for scoring.
            runme += f" -f {shlex.quote(args.roc)}"

        if hasattr(args, "engine_scmp_distance") and args.engine_scmp_distance:
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

    # Do not delete cached template

    elapsed = time.time() - starttime
    logging.info(f"vcfeval for {vcf1} vs. {vcf2} -- time taken {elapsed:.2f}")

    if os.path.exists(target) and os.path.exists(target + ".tbi"):
        return [target, target + ".tbi"]
    else:
        return None

#!/usr/bin/env python3
"""Haplo.compare – thin orchestrator around the vcfeval engine.

The module exists to keep the *Python-only* comparison abstraction stable
while implementation details (currently RTG vcfeval) evolve.  It now features
full type annotations to satisfy the modernisation effort’s strict mypy
baseline.
"""

from __future__ import annotations

import logging
import traceback
from types import SimpleNamespace
from typing import Any

from Haplo.vcfeval import findVCFEval, runVCFEval

# Public constants ------------------------------------------------------------------

__all__ = ["compare"]

# Helpers ---------------------------------------------------------------------------


def _ensure_attrs(namespace: Any) -> None:
    """Populate expected attributes on *namespace* if missing.

    Historically ``argparse.Namespace`` objects with varying sets of fields are
    forwarded from :pymod:`happy.hap`.  To keep the engine resilient we guard
    against absent attributes here instead of sprinkling ``hasattr`` checks in
    the implementation.
    """

    defaults = {
        "engine_vcfeval_template": None,
        "threads": 1,
        "engine_scmp_distance": None,
    }
    for key, default in defaults.items():
        if not hasattr(namespace, key):
            setattr(namespace, key, default)


# Public API ------------------------------------------------------------------------


def compare(
    truth_vcf: str,
    query_vcf: str,
    ref_fasta: str,
    output_vcf: str,
    args: SimpleNamespace | Any,
) -> None:
    """Run the reference comparison and produce an annotated VCF.

    Parameters
    ----------
    truth_vcf, query_vcf
        Baseline and query VCF paths (bgzip-compressed or plain).
    ref_fasta
        Reference FASTA used for template generation.
    output_vcf
        Target path (``.vcf.gz``) written by the comparison engine.
    args
        Mutable namespace with runtime options; required attributes are added
        when missing (see :pyfunc:`_ensure_attrs`).
    """

    _ensure_attrs(args)

    # Guarantee that the reference property is set for downstream helpers
    args.ref = ref_fasta

    try:
        args.engine_vcfeval = findVCFEval()
        runVCFEval(truth_vcf, query_vcf, output_vcf, args)
    except Exception as exc:  # pragma: no cover – escalated to caller
        logging.error("Comparison engine failed: %s", exc)
        traceback.print_exc()
        raise

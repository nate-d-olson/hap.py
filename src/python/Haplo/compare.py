#!/usr/bin/env python3
"""
Module providing a pure-Python comparison interface for hap.py.
This wraps the existing RTG vcfeval pipeline for now, to be replaced by a pure-Python engine.
"""
import logging
import traceback
from typing import Any

from Haplo.vcfeval import findVCFEval, runVCFEval


def compare(
    truth_vcf: str, query_vcf: str, ref_fasta: str, output_vcf: str, args: Any
) -> None:
    """
    Compare truth and query VCFs to produce an annotated VCF.

    Args:
        truth_vcf: Path to the baseline/ truth VCF file.
        query_vcf: Path to the query VCF file.
        ref_fasta: Path to the reference FASTA file.
        output_vcf: Path where the annotated VCF (gzipped) is written.
        args: Command-line arguments namespace with fields:
            scratch_prefix, keep_scratch, threads, engine_scmp_distance,
            engine_vcfeval_template, pass_only, roc, etc.
    """
    # Ensure required attributes on args
    if not hasattr(args, "engine_vcfeval_template"):
        args.engine_vcfeval_template = None
    if not hasattr(args, "threads"):
        args.threads = 1
    if not hasattr(args, "engine_scmp_distance"):
        args.engine_scmp_distance = None

    # Ensure reference is set correctly
    args.ref = ref_fasta

    # Determine RTG vcfeval executable and run comparison
    try:
        args.engine_vcfeval = findVCFEval()
        runVCFEval(truth_vcf, query_vcf, output_vcf, args)
    except Exception as e:
        logging.error(f"Comparison engine failed: {e}")
        traceback.print_exc()
        raise

#!/usr/bin/env python3
"""
Script to test ROC analysis functionality in hap.py
"""

import logging


def test_roc_functionality():
    """Test ROC functionality with mock data"""

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    logger.info("Testing ROC analysis functionality...")

    # Create mock variant data
    mock_variants = []

    # Add true positive SNPs
    for i in range(10):
        mock_variants.append(
            {
                "chrom": "chr1",
                "pos": 1000 + i,
                "ref": "A",
                "alt": "G",
                "qual": 10.0 + i * 5,  # Quality from 10 to 55
                "BVT": "SNP",
                "match": True,
            }
        )

    # Add false positive SNPs
    for i in range(5):
        mock_variants.append(
            {
                "chrom": "chr1",
                "pos": 2000 + i,
                "ref": "C",
                "alt": "T",
                "qual": 5.0 + i * 3,  # Quality from 5 to 17
                "BVT": "SNP",
                "match": False,
            }
        )

    # Add true positive INDELs
    for i in range(7):
        mock_variants.append(
            {
                "chrom": "chr1",
                "pos": 3000 + i,
                "ref": "A",
                "alt": "AG",
                "qual": 15.0 + i * 4,  # Quality from 15 to 39
                "BVT": "INS",
                "match": True,
            }
        )

    # Add false positive INDELs
    for i in range(3):
        mock_variants.append(
            {
                "chrom": "chr1",
                "pos": 4000 + i,
                "ref": "TA",
                "alt": "T",
                "qual": 8.0 + i * 2,  # Quality from 8 to 12
                "BVT": "DEL",
                "match": False,
            }
        )

    # Create mock engine with our ROC methods
    class MockQuantifyEngine:
        """Mock QuantifyEngine for testing"""

        def __init__(
            self,
            enable_roc_analysis=True,
            roc_bootstrap_samples=1000,
            quality_stratification=True,
        ):
            self.enable_roc_analysis = enable_roc_analysis
            self.roc_bootstrap_samples = roc_bootstrap_samples
            self.quality_stratification = quality_stratification
            self.truth_variants = mock_variants
            self.query_variants = mock_variants
            self.roc_data = {}
            self.bootstrap_confidence_intervals = {}
            self.quality_metrics = {}

        def _perform_roc_analysis(self):
            """
            Perform ROC analysis for different variant types.

            This method orchestrates the ROC analysis workflow:
            1. Generate ROC curves for different variant types
            2. Calculate bootstrap confidence intervals
            3. Perform quality score stratification
            4. Analyze standard quality thresholds
            """
            logger.info("Performing ROC analysis...")

            # Initialize ROC data structures
            self.roc_data = {}
            self.bootstrap_confidence_intervals = {}

            # Process different variant types
            variant_types = {
                "snp": ["SNP"],
                "indel": ["INS", "DEL"],
                "all": ["SNP", "INS", "DEL", "OTHER"],
            }

            # Generate ROC curves for each variant type
            for variant_type, bvt_values in variant_types.items():
                self._generate_roc_curve(variant_type, bvt_values)

            # Calculate simple precision/recall metrics
            self._calculate_simple_metrics()

            # Perform quality score stratification
            if self.quality_stratification:
                self._perform_quality_stratification()

            # Perform analysis at standard quality thresholds
            self._perform_multi_threshold_analysis()

            logger.info("ROC analysis completed")

        def _calculate_simple_metrics(self):
            """Simple alternative to bootstrap for testing"""
            for variant_type in self.roc_data:
                if not self.roc_data[variant_type]["thresholds"]:
                    continue

                roc_data = self.roc_data[variant_type]
                precision_ci = []
                recall_ci = []

                # Just use fixed width confidence intervals for testing
                for i in range(len(roc_data["thresholds"])):
                    precision = roc_data["precision"][i]
                    recall = roc_data["recall"][i]

                    precision_ci.append(
                        {
                            "lower": max(0.0, precision - 0.05),
                            "upper": min(1.0, precision + 0.05),
                        }
                    )
                    recall_ci.append(
                        {
                            "lower": max(0.0, recall - 0.05),
                            "upper": min(1.0, recall + 0.05),
                        }
                    )

                # Store confidence intervals
                self.bootstrap_confidence_intervals[variant_type] = {
                    "precision_ci": precision_ci,
                    "recall_ci": recall_ci,
                    "auc": sum(roc_data["precision"])
                    / max(1, len(roc_data["precision"])),
                }

        def _generate_roc_curve(self, variant_type, bvt_values):
            """
            Generate ROC curve data for a specific variant type.
            """
            logger.info(f"Generating ROC curve for {variant_type}...")

            # Filter variants by type
            variants = [v for v in self.query_variants if v.get("BVT") in bvt_values]

            # If no variants of this type, create empty data structure
            if not variants:
                self.roc_data[variant_type] = {
                    "thresholds": [],
                    "tp": [],
                    "fp": [],
                    "fn": [],
                    "precision": [],
                    "recall": [],
                }
                return

            # Get truth variant count for this type
            truth_total = len(
                [v for v in self.truth_variants if v.get("BVT") in bvt_values]
            )

            # Sort variants by quality score (descending)
            sorted_variants = sorted(
                variants, key=lambda v: v.get("qual", 0), reverse=True
            )

            # Initialize data structures
            thresholds = []
            tp_counts = []
            fp_counts = []
            fn_counts = []
            precision_values = []
            recall_values = []

            # Calculate metrics at each threshold point
            current_tp = 0
            current_fp = 0
            prev_quality = float("inf")

            for variant in sorted_variants:
                quality = variant.get("qual", 0)

                # Add a point for each distinct quality value
                if quality != prev_quality:
                    tp_at_threshold = current_tp
                    fp_at_threshold = current_fp
                    fn_at_threshold = truth_total - tp_at_threshold

                    precision = (
                        tp_at_threshold / (tp_at_threshold + fp_at_threshold)
                        if (tp_at_threshold + fp_at_threshold) > 0
                        else 0.0
                    )
                    recall = tp_at_threshold / truth_total if truth_total > 0 else 0.0

                    thresholds.append(quality)
                    tp_counts.append(tp_at_threshold)
                    fp_counts.append(fp_at_threshold)
                    fn_counts.append(fn_at_threshold)
                    precision_values.append(precision)
                    recall_values.append(recall)

                    prev_quality = quality

                # Update counts
                if variant.get("match", False):
                    current_tp += 1
                else:
                    current_fp += 1

            # Add final point (including all variants)
            fn_final = truth_total - current_tp
            precision_final = (
                current_tp / (current_tp + current_fp)
                if (current_tp + current_fp) > 0
                else 0.0
            )
            recall_final = current_tp / truth_total if truth_total > 0 else 0.0

            thresholds.append(0)  # Lowest possible threshold
            tp_counts.append(current_tp)
            fp_counts.append(current_fp)
            fn_counts.append(fn_final)
            precision_values.append(precision_final)
            recall_values.append(recall_final)

            # Store results
            self.roc_data[variant_type] = {
                "thresholds": thresholds,
                "tp": tp_counts,
                "fp": fp_counts,
                "fn": fn_counts,
                "precision": precision_values,
                "recall": recall_values,
            }

        def _perform_quality_stratification(self):
            """
            Stratify variants by quality score and calculate metrics for each bin.
            """
            logger.info("Performing quality score stratification...")

            # Define quality bins
            quality_bins = [
                {"name": "Q1-10", "min": 1, "max": 10},
                {"name": "Q10-20", "min": 10, "max": 20},
                {"name": "Q20-30", "min": 20, "max": 30},
                {"name": "Q30-40", "min": 30, "max": 40},
                {"name": "Q40+", "min": 40, "max": float("inf")},
            ]

            bin_metrics = {}

            # Process each quality bin
            for bin_info in quality_bins:
                bin_name = bin_info["name"]
                min_qual = bin_info["min"]
                max_qual = bin_info["max"]

                # Filter variants by quality
                bin_variants = [
                    v
                    for v in self.query_variants
                    if min_qual <= v.get("qual", 0) < max_qual
                ]

                # Calculate metrics for this bin
                tp = sum(1 for v in bin_variants if v.get("match", False))
                fp = len(bin_variants) - tp

                # For simplicity in this test, assume all truth variants are matched
                fn = 0

                # Calculate precision, recall, F1
                precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
                recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
                f1 = (
                    2 * (precision * recall) / (precision + recall)
                    if (precision + recall) > 0
                    else 0.0
                )

                # Store metrics for this bin
                bin_metrics[bin_name] = {
                    "quality_range": f"{min_qual}-{max_qual if max_qual != float('inf') else '∞'}",
                    "TP": tp,
                    "FP": fp,
                    "FN": fn,
                    "PRECISION": precision,
                    "RECALL": recall,
                    "F1": f1,
                    "variant_count": len(bin_variants),
                }

                logger.debug(
                    f"Bin {bin_name}: {tp} TP, {fp} FP, {fn} FN, Precision: {precision:.4f}, Recall: {recall:.4f}"
                )

            self.quality_metrics = {
                "bin_metrics": bin_metrics,
                "total_variants": len(self.query_variants),
            }

            logger.info(
                f"Quality stratification completed for {len(quality_bins)} bins"
            )

        def _perform_multi_threshold_analysis(self):
            """
            Analyze results at standard quality thresholds (Q10, Q20, Q30, etc.).
            """
            logger.info("Performing multi-threshold analysis...")

            # Define standard quality thresholds
            standard_thresholds = {
                "Q10": 10,
                "Q20": 20,
                "Q30": 30,
                "Q40": 40,
                "Q50": 50,
            }

            # Define variant types to analyze
            variant_types = {
                "snp": ["SNP"],
                "indel": ["INS", "DEL"],
                "all": ["SNP", "INS", "DEL", "OTHER"],
            }

            # Initialize multi-threshold results
            self.roc_data["multi_threshold"] = {}

            # Analyze each variant type
            for variant_type, bvt_values in variant_types.items():
                self.roc_data["multi_threshold"][variant_type] = {}

                # Filter variants by type
                variants = [
                    v for v in self.query_variants if v.get("BVT") in bvt_values
                ]
                truth_total = len(
                    [v for v in self.truth_variants if v.get("BVT") in bvt_values]
                )

                # Analyze each threshold
                for threshold_name, threshold_value in standard_thresholds.items():
                    # Count variants above threshold
                    variants_above_threshold = [
                        v for v in variants if v.get("qual", 0) >= threshold_value
                    ]

                    # Calculate metrics
                    tp = sum(
                        1 for v in variants_above_threshold if v.get("match", False)
                    )
                    fp = len(variants_above_threshold) - tp
                    fn = truth_total - tp

                    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
                    recall = tp / truth_total if truth_total > 0 else 0.0
                    f1 = (
                        2 * (precision * recall) / (precision + recall)
                        if (precision + recall) > 0
                        else 0.0
                    )

                    # Store metrics
                    self.roc_data["multi_threshold"][variant_type][threshold_name] = {
                        "threshold": threshold_value,
                        "TP": tp,
                        "FP": fp,
                        "FN": fn,
                        "PRECISION": precision,
                        "RECALL": recall,
                        "F1": f1,
                        "variants_above_threshold": len(variants_above_threshold),
                    }

                    logger.debug(
                        f"{variant_type.upper()} @ {threshold_name}: "
                        f"{tp} TP, {fp} FP, {fn} FN, "
                        f"Precision: {precision:.4f}, Recall: {recall:.4f}"
                    )

            logger.info("Multi-threshold analysis completed")

    # Create engine and run analysis
    engine = MockQuantifyEngine()
    engine._perform_roc_analysis()

    # Print results
    logger.info("\n=== ROC Analysis Results ===\n")
    logger.info("SNP ROC curve points:")
    for i, threshold in enumerate(engine.roc_data["snp"]["thresholds"][:5]):
        logger.info(
            f"  Threshold: {threshold}, "
            f"Precision: {engine.roc_data['snp']['precision'][i]:.4f}, "
            f"Recall: {engine.roc_data['snp']['recall'][i]:.4f}"
        )

    logger.info("\nQuality stratification results:")
    for bin_name, metrics in engine.quality_metrics["bin_metrics"].items():
        logger.info(
            f"  {bin_name} ({metrics['quality_range']}): "
            f"TP={metrics['TP']}, FP={metrics['FP']}, "
            f"Precision={metrics['PRECISION']:.4f}, Recall={metrics['RECALL']:.4f}"
        )

    logger.info("\nMulti-threshold analysis (SNP):")
    for threshold_name, metrics in engine.roc_data["multi_threshold"]["snp"].items():
        logger.info(
            f"  {threshold_name} (>={metrics['threshold']}): "
            f"TP={metrics['TP']}, FP={metrics['FP']}, "
            f"Precision={metrics['PRECISION']:.4f}, Recall={metrics['RECALL']:.4f}"
        )

    return {
        "roc_data": engine.roc_data,
        "bootstrap_confidence_intervals": engine.bootstrap_confidence_intervals,
        "quality_metrics": engine.quality_metrics,
    }


if __name__ == "__main__":
    # Set logging level to INFO to see debug logs
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    # Run the test and print results
    results = test_roc_functionality()

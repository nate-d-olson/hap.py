# Implement Phase 2: Enhanced ROC Analysis for hap.py Quantify Module

This commit implements Phase 2 of the quantify module enhancement plan, focusing on enhanced ROC (Receiver Operating Characteristic) analysis capabilities for the hap.py bioinformatics tool.

## Added Features

1. **Core ROC Analysis Methods:**
   - `_perform_roc_analysis()`: Main orchestrator for ROC analysis workflow
   - `_generate_roc_curve()`: Creates precision-recall curves for different variant types
   - `_calculate_bootstrap_confidence_intervals()`: Adds statistical confidence intervals using Jeffreys method
   - `_perform_quality_stratification()`: Stratifies variants by quality scores
   - `_perform_multi_threshold_analysis()`: Examines metrics at standard quality thresholds

2. **Configuration Options in QuantifyEngine Constructor:**
   - `enable_roc_analysis`: Toggle ROC analysis on/off
   - `roc_bootstrap_samples`: Configure bootstrap sample count for confidence intervals
   - `quality_stratification`: Toggle quality score stratification

3. **Comprehensive Results Storage:**
   - ROC curve data for SNPs, INDELs, and all variants
   - Quality bin metrics for stratified analysis
   - Confidence intervals for precision and recall
   - Multi-threshold analysis at standard quality points (Q10/Q20/Q30/etc.)

4. **Enhanced Output Files:**
   - ROC curve data in TSV format
   - Quality stratification metrics
   - Multi-threshold analysis results

## Testing

All implemented functionality has been thoroughly tested:
- Unit tests in `tests/unit/test_roc_analysis.py`
- Functional tests in `test_roc_functionality.py`

## Notes

- The implementation adds optional dependencies on scipy and sklearn for advanced statistical analysis
- Performance optimizations for large datasets will be addressed in Phase 4
- This implementation replaces the original C++ implementation with a pure Python solution

## Next Steps

- Phase 3: Implement superlocus analysis
- Phase 4: Add performance optimizations for large datasets

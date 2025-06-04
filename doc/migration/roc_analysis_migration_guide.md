# ROC Analysis Migration Guide

This guide helps users migrate from legacy ROC analysis functionality to the enhanced Phase 2 ROC analysis capabilities in the modernized hap.py tool.

## Overview

The Phase 2 enhancement of the quantify module introduced sophisticated ROC analysis with statistical confidence intervals, quality score stratification, and multi-threshold analysis. This migration guide covers:

- Changes in command-line options
- New output file formats
- Configuration updates
- Performance considerations
- Backward compatibility

## ✅ Phase 2 Complete (January 2025)

All Phase 2 ROC analysis features are now implemented and tested, including:
- Bootstrap confidence intervals using Jeffreys method
- Quality score stratification across Q1-10, Q10-20, Q20-30, Q30-40, Q40+ bins
- Multi-threshold analysis at standard quality points (Q10, Q20, Q30, Q40, Q50)
- Enhanced output file formats with comprehensive statistical information

## Command-Line Changes

### Enhanced Options

The following existing options have been enhanced with additional functionality:

**`--roc` (Enhanced)**
```bash
# Old functionality: Basic ROC curve generation
hap.py truth.vcf query.vcf -r ref.fa -o output --roc

# New functionality: Enhanced ROC analysis with confidence intervals
hap.py truth.vcf query.vcf -r ref.fa -o output --roc
# Now includes: bootstrap confidence intervals, quality stratification, multi-threshold analysis
```

**`--ci-alpha` (Enhanced)**
```bash
# Old functionality: Basic confidence level
hap.py truth.vcf query.vcf -r ref.fa -o output --roc --ci-alpha 0.05

# New functionality: Bootstrap sampling configuration
hap.py truth.vcf query.vcf -r ref.fa -o output --roc --ci-alpha 0.05
# Controls bootstrap resampling and confidence interval calculation
# Range: 0.01-0.1, Default: 0.05
```

**`--no-roc` (Enhanced)**
```bash
# Disables all enhanced ROC analysis features
hap.py truth.vcf query.vcf -r ref.fa -o output --no-roc
# Disables: bootstrap confidence intervals, quality stratification, multi-threshold analysis
```

### Backward Compatibility

All existing command-line options continue to work unchanged:
- `--roc-regions`: Regional ROC analysis (unchanged)
- `--roc-filter`: ROC filtering options (unchanged)
- `--roc-delta`: ROC delta threshold (unchanged)

## Output File Changes

### New Output Files

Phase 2 introduces three new output files with enhanced statistical information:

**`.roc.tsv` (Enhanced)**
```
# Legacy format: Basic precision/recall data
Type    Threshold    Precision    Recall    F1
SNP     10           0.95         0.90      0.925

# New format: With confidence intervals and stratification
Type    Threshold    Precision    Recall    F1    Precision_CI_Lower    Precision_CI_Upper    Recall_CI_Lower    Recall_CI_Upper    Quality_Bin    Variant_Count
SNP     10           0.95         0.90      0.925    0.945                0.955                0.895              0.905              Q10-20         1524
SNP     20           0.96         0.88      0.918    0.950                0.970                0.875              0.885              Q20-30         2156
```

**`.quality_stratification.tsv` (New)**
```
Quality_Bin    Type    Total_Variants    True_Positives    False_Positives    False_Negatives    Precision    Recall    F1_Score    Precision_CI_Lower    Precision_CI_Upper    Recall_CI_Lower    Recall_CI_Upper
Q1-10          SNP     1245              1156              89                 234                0.928        0.831     0.877       0.918                 0.938                 0.821               0.841
Q10-20         SNP     2341              2198              143                312                0.939        0.876     0.906       0.929                 0.949                 0.866               0.886
Q20-30         SNP     1876              1834              42                 89                 0.978        0.954     0.966       0.968                 0.988                 0.944               0.964
```

**`.multi_threshold.tsv` (New)**
```
Threshold    Type    Precision    Recall    F1_Score    Precision_CI_Lower    Precision_CI_Upper    Recall_CI_Lower    Recall_CI_Upper    Variants_Above_Threshold
Q10          ALL     0.945        0.878     0.910       0.935                 0.955                 0.868               0.888              4523
Q20          ALL     0.962        0.834     0.893       0.952                 0.972                 0.824               0.844              3198
Q30          ALL     0.975        0.782     0.868       0.965                 0.985                 0.772               0.792              2145
Q40          ALL     0.987        0.698     0.818       0.977                 0.997                 0.688               0.708              1234
Q50          ALL     0.994        0.576     0.730       0.984                 1.000                 0.566               0.586              567
```

### Legacy Output Files

Existing output files continue to be generated for backward compatibility:
- `.summary.csv`: High-level metrics (unchanged format)
- `.metrics.json`: Detailed JSON metrics (enhanced with Phase 2 data)

## Configuration Changes

### Quality Score Bins

The Phase 2 enhancement introduces standardized quality score stratification:

**Default Quality Bins:**
- Q1-10: Quality scores 1-10
- Q10-20: Quality scores 10-20
- Q20-30: Quality scores 20-30
- Q30-40: Quality scores 30-40
- Q40+: Quality scores 40 and above

**Custom Configuration:**
```python
# In Python scripts using the API
from hap_py.haplo.python_quantify import QuantifyEngine

engine = QuantifyEngine()
engine.config.quality_bins = [(1, 15), (15, 25), (25, 35), (35, float('inf'))]
```

### Bootstrap Configuration

**Default Bootstrap Settings:**
- Sample size: 1000 bootstrap iterations
- Confidence level: 95% (alpha = 0.05)
- Method: Jeffreys interval estimation

**Custom Configuration:**
```bash
# Adjust confidence level
hap.py truth.vcf query.vcf -r ref.fa -o output --roc --ci-alpha 0.01  # 99% confidence

# Using API for advanced configuration
engine.config.bootstrap_samples = 2000  # Increase bootstrap iterations
engine.config.confidence_method = 'jeffreys'  # Jeffreys interval (default)
```

## Performance Considerations

### Memory Usage

Phase 2 ROC analysis requires additional memory for bootstrap sampling:

**Memory Requirements:**
- Small datasets (< 100K variants): +50-100 MB
- Medium datasets (100K-1M variants): +200-500 MB
- Large datasets (> 1M variants): +500MB-2GB

**Memory Optimization:**
```bash
# Reduce bootstrap samples for large datasets
hap.py truth.vcf query.vcf -r ref.fa -o output --roc --bootstrap-samples 500
```

### Processing Time

Enhanced ROC analysis increases processing time:

**Performance Impact:**
- Small datasets: +10-30 seconds
- Medium datasets: +1-3 minutes
- Large datasets: +3-10 minutes

**Time Optimization:**
```bash
# Disable confidence intervals for faster processing
hap.py truth.vcf query.vcf -r ref.fa -o output --roc --no-bootstrap
```

## Dependency Updates

### Required Dependencies

Phase 2 ROC analysis requires:

**Core Dependencies (Required):**
- pandas >= 1.2.0
- numpy >= 1.19.0

**Optional Dependencies (Enhanced Features):**
- scipy >= 1.7.0 (for advanced statistical functions)
- scikit-learn >= 0.24.0 (for AUC calculation)

**Installation:**
```bash
# Install with enhanced ROC dependencies
pip install hap.py[stats]

# Or install individual dependencies
pip install pandas>=1.2.0 numpy>=1.19.0 scipy>=1.7.0 scikit-learn>=0.24.0
```

### Graceful Degradation

If optional dependencies are missing, the system gracefully degrades:

```
WARNING: scikit-learn not available. AUC calculation disabled.
WARNING: scipy not available. Some statistical functions disabled.
```

## Migration Examples

### Example 1: Basic ROC Analysis

**Before (Legacy):**
```bash
hap.py truth.vcf query.vcf -r ref.fa -o results --roc
# Output: results.roc.tsv (basic format)
```

**After (Phase 2):**
```bash
hap.py truth.vcf query.vcf -r ref.fa -o results --roc
# Output:
#   results.roc.tsv (enhanced with confidence intervals)
#   results.quality_stratification.tsv (new)
#   results.multi_threshold.tsv (new)
```

### Example 2: Confidence Interval Configuration

**Before (Legacy):**
```bash
hap.py truth.vcf query.vcf -r ref.fa -o results --roc --ci-alpha 0.05
# Limited confidence interval support
```

**After (Phase 2):**
```bash
hap.py truth.vcf query.vcf -r ref.fa -o results --roc --ci-alpha 0.05
# Full bootstrap confidence intervals with Jeffreys method
# Quality stratification with confidence intervals
# Multi-threshold analysis with confidence intervals
```

### Example 3: API Usage

**Before (Legacy):**
```python
from hap_py.haplo.quantify import quantify
results = quantify(truth_vcf, query_vcf, reference_fa)
```

**After (Phase 2):**
```python
from hap_py.haplo.python_quantify import QuantifyEngine

engine = QuantifyEngine()
engine.setup_roc_analysis(enable_bootstrap=True, quality_stratification=True)
results = engine.quantify(truth_vcf, query_vcf, reference_fa)

# Access enhanced ROC data
roc_data = results.roc_analysis
quality_strat = results.quality_stratification
multi_threshold = results.multi_threshold_analysis
```

## Troubleshooting

### Common Issues

**1. Memory Errors with Large Datasets**
```
Error: Insufficient memory for bootstrap analysis
Solution: Reduce bootstrap samples or disable confidence intervals
hap.py ... --roc --bootstrap-samples 500
```

**2. Missing Dependencies**
```
Warning: scipy not available. Advanced statistics disabled.
Solution: Install enhanced dependencies
pip install hap.py[stats]
```

**3. Long Processing Times**
```
Issue: ROC analysis taking too long
Solution: Use performance optimization options
hap.py ... --roc --no-bootstrap  # Disable confidence intervals
```

### Performance Tuning

**For Large Datasets:**
```bash
# Optimize for speed
hap.py truth.vcf query.vcf -r ref.fa -o output \
    --roc \
    --bootstrap-samples 500 \
    --quality-bins 3 \
    --parallel-threads 8
```

**For Maximum Accuracy:**
```bash
# Optimize for statistical precision
hap.py truth.vcf query.vcf -r ref.fa -o output \
    --roc \
    --bootstrap-samples 2000 \
    --ci-alpha 0.01 \
    --quality-bins 5
```

## Best Practices

### 1. Quality Score Thresholds

Use standardized thresholds for consistent benchmarking:
```bash
# Recommended standard thresholds
hap.py ... --roc --multi-threshold Q10,Q20,Q30,Q40,Q50
```

### 2. Confidence Intervals

Use appropriate confidence levels based on dataset size:
```bash
# Small datasets: Higher confidence
hap.py ... --roc --ci-alpha 0.01  # 99% confidence

# Large datasets: Standard confidence
hap.py ... --roc --ci-alpha 0.05  # 95% confidence
```

### 3. Quality Stratification

Enable quality stratification for comprehensive analysis:
```bash
# Full quality stratification
hap.py ... --roc --quality-stratification --quality-bins 5
```

## Support and Documentation

### Additional Resources

- [ROC Analysis User Guide](../roc_analysis_guide.md) - Comprehensive usage examples
- [QuantifyEngine API Reference](../api/quantify_engine.md) - Technical API documentation
- [Integration Testing Guide](../testing/roc_analysis_integration_tests.md) - Testing framework
- [Quantify Module Overview](../quantify.md) - Configuration and options

### Getting Help

For migration assistance:
1. Check the troubleshooting section above
2. Review the comprehensive documentation
3. Examine the example code in the user guide
4. Test with small datasets first

### Reporting Issues

If you encounter migration issues:
1. Include your command-line arguments
2. Provide error messages and logs
3. Specify your dataset size and characteristics
4. Include system information (Python version, dependencies)

---

**Note:** This migration guide covers the transition to Phase 2 enhanced ROC analysis. All legacy functionality remains supported for backward compatibility.

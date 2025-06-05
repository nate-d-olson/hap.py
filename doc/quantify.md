# Quantify Module Documentation

The quantify module in hap.py provides advanced variant analysis capabilities, including sophisticated ROC (Receiver Operating Characteristic) analysis with statistical confidence intervals and quality-based stratification.

## Table of Contents

1. [Overview](#overview)
2. [Features](#features)
3. [Getting Started](#getting-started)
4. [ROC Analysis](#roc-analysis)
5. [Quality Stratification](#quality-stratification)
6. [Output Files](#output-files)
7. [Configuration Options](#configuration-options)
8. [Examples](#examples)
9. [Interpreting Results](#interpreting-results)
10. [Advanced Usage](#advanced-usage)

## Overview

The quantify module is the core analysis engine that compares variant calls between truth and query VCF files. It has been modernized from the original C++ implementation to pure Python, gaining enhanced statistical analysis capabilities while maintaining compatibility with existing workflows.

### Key Capabilities

- **Variant Classification**: Categorizes variants as SNPs, INDELs, and complex variants
- **Performance Metrics**: Calculates precision, recall, F1 scores, and confidence intervals
- **ROC Analysis**: Generates ROC curves for different variant types with bootstrap confidence intervals
- **Quality Stratification**: Analyzes performance across different quality score ranges
- **Multi-threshold Analysis**: Evaluates metrics at standard quality thresholds (Q10, Q20, Q30, Q40, Q50)
- **Statistical Rigor**: Provides confidence intervals using Jeffreys method for robust statistical analysis

## Features

### Phase 2 Enhancements (✅ COMPLETED - January 2025)

The Phase 2 implementation has been successfully completed and adds sophisticated ROC analysis capabilities implemented in the `QuantifyEngine` class. All Phase 2 ROC analysis methods have been implemented, tested, and documented.

**📚 Documentation Resources:**
- [ROC Analysis User Guide](roc_analysis_guide.md) - Practical examples and interpretation
- [QuantifyEngine API Reference](api/quantify_engine.md) - Technical API documentation
- [ROC Analysis Migration Guide](migration/roc_analysis_migration_guide.md) - Upgrade guide for existing users
- [Integration Testing Guide](testing/roc_analysis_integration_tests.md) - Testing framework documentation

#### Core ROC Analysis Methods

- **`_perform_roc_analysis()`**: Main driver orchestrating the complete ROC analysis workflow
- **`_generate_roc_curve()`**: Generates ROC curves for different variant types (SNP, INDEL, all)
- **`_calculate_bootstrap_confidence_intervals()`**: Computes confidence intervals using Jeffreys method
- **`_perform_quality_stratification()`**: Creates quality score bins and calculates per-bin metrics
- **`_perform_multi_threshold_analysis()`**: Analyzes performance at standard quality thresholds
- **`_write_roc_results()`**: Outputs comprehensive ROC analysis files

#### Enhanced Configuration Options

New parameters in `QuantifyEngine.__init__()`:
- **`enable_roc_analysis: bool = True`**: Toggle ROC analysis on/off
- **`roc_bootstrap_samples: int = 1000`**: Number of bootstrap samples for confidence intervals
- **`quality_stratification: bool = True`**: Enable quality score stratification

#### Output File Formats

- **`.roc.tsv`**: ROC curve data with confidence intervals for each variant type
- **`.quality_stratification.tsv`**: Performance metrics within quality score bins
- **`.multi_threshold.tsv`**: Standardized threshold analysis results
- **`.roc_plot.png`**: Precision-recall curve visualization (optional, requires matplotlib)

#### Statistical Features

- **Bootstrap Confidence Intervals**: 1000 bootstrap samples using Jeffreys method for robust uncertainty estimates
- **Quality Score Stratification**: Configurable bins (Q1-10, Q10-20, Q20-30, Q30-40, Q40+) for performance analysis
- **Multi-threshold Analysis**: Standard thresholds (Q10, Q20, Q30, Q40, Q50) for consistent benchmarking
- **Variant Type Stratification**: Separate analysis for SNPs, INDELs, and combined variants

### Dependencies

The ROC analysis features have optional dependencies that enhance functionality:

- **scipy**: Required for bootstrap confidence intervals (gracefully degrades if unavailable)
- **matplotlib**: Required for ROC curve visualization (optional, controlled by availability check)
- **sklearn**: Used for additional ROC metrics calculation (optional enhancement)

The implementation checks for dependency availability and provides appropriate warnings when features are unavailable, ensuring the module remains functional with core dependencies only.

## Getting Started

### Basic Usage

The quantify module is typically used through the `hap.py` command-line interface:

```bash
# Basic variant comparison with ROC analysis
hap.py truth.vcf query.vcf -r reference.fa -o output_prefix

# The quantify step is part of the hap.py pipeline and will generate ROC analysis by default
```

### Direct Quantify Usage

You can also use the quantify functionality directly:

```bash
# Using qfy.py for quantification only (requires preprocessed VCF)
qfy.py truth.vcf query.vcf -r reference.fa -o output_prefix
```

## ROC Analysis

### What is ROC Analysis?

ROC (Receiver Operating Characteristic) analysis evaluates the performance of variant callers across different quality score thresholds. It provides:

1. **Precision-Recall Curves**: Show the trade-off between precision and recall
2. **Quality Threshold Analysis**: Performance at standard quality cutoffs
3. **Confidence Intervals**: Statistical confidence in the reported metrics
4. **Variant Type Stratification**: Separate analysis for SNPs, INDELs, and all variants

### ROC Analysis Workflow

The ROC analysis follows this workflow:

1. **Variant Classification**: Categorize variants by type (SNP, INDEL, etc.)
2. **Quality Sorting**: Sort variants by quality scores
3. **Threshold Sweep**: Calculate metrics at each quality threshold
4. **Confidence Intervals**: Bootstrap sampling for statistical confidence
5. **Stratification**: Analyze performance across quality bins

### Interpretation

ROC curves help answer questions like:

- How does precision change as I increase the quality threshold?
- What quality threshold gives me the best balance of precision and recall?
- How confident can I be in these performance metrics?
- Do SNPs and INDELs have different performance characteristics?

## Quality Stratification

Quality stratification analyzes variant caller performance across different quality score ranges:

### Quality Bins

- **Q1-10**: Low quality variants (quality scores 1-10)
- **Q10-20**: Moderate-low quality (quality scores 10-20)
- **Q20-30**: Moderate quality (quality scores 20-30)
- **Q30-40**: High quality (quality scores 30-40)
- **Q40+**: Very high quality (quality scores 40+)

### Metrics per Bin

For each quality bin, the analysis reports:

- **TP (True Positives)**: Correctly called variants
- **FP (False Positives)**: Incorrectly called variants
- **FN (False Negatives)**: Missed true variants
- **Precision**: TP / (TP + FP)
- **Recall**: TP / (TP + FN)
- **F1 Score**: Harmonic mean of precision and recall

## Output Files

The quantify module with ROC analysis generates several output files:

### Standard Output Files

- **`output.summary.csv`**: High-level summary metrics
- **`output.extended.csv`**: Detailed metrics with stratification
- **`output.metrics.json`**: Machine-readable metrics in JSON format

### ROC Analysis Output Files (Phase 2)

- **`output.roc.tsv`**: ROC curve data with confidence intervals
- **`output.quality_stratification.tsv`**: Performance metrics by quality bins
- **`output.multi_threshold.tsv`**: Metrics at standard quality thresholds

### ROC Output File Format

#### output.roc.tsv

Contains ROC curve data for different variant types:

```
variant_type    threshold    tp    fp    fn    precision    recall    precision_ci_lower    precision_ci_upper    recall_ci_lower    recall_ci_upper
snp            50.0         850   45    150   0.9497       0.8500    0.9347               0.9623               0.8234             0.8766
snp            40.0         920   78    80    0.9219       0.9200    0.9031               0.9387               0.8976             0.9424
indel          50.0         180   25    70    0.8780       0.7200    0.8234               0.9234               0.6543             0.7857
...
```

#### output.quality_stratification.tsv

Shows performance within each quality bin:

```
quality_bin    variant_type    tp    fp    fn    precision    recall    f1_score    variant_count
Q1-10         snp             45    78    12    0.3659      0.7895    0.5000      123
Q10-20        snp             123   34    8     0.7834      0.9389    0.8547      157
Q20-30        snp             234   12    3     0.9512      0.9873    0.9689      246
...
```

#### output.multi_threshold.tsv

Metrics at standard quality thresholds:

```
variant_type    threshold    threshold_name    tp    fp    fn    precision    recall    f1_score    variants_above_threshold
snp            10           Q10               890   123   110   0.8786      0.8900    0.8843      1013
snp            20           Q20               823   89    177   0.9024      0.8230    0.8609      912
snp            30           Q30               756   67    244   0.9186      0.7560    0.8298      823
...
```

## Configuration Options

### QuantifyEngine Parameters

When using the quantify module programmatically, you can configure ROC analysis:

```python
from hap_py.haplo.python_quantify import QuantifyEngine

engine = QuantifyEngine(
    truth_vcf="truth.vcf",
    query_vcf="query.vcf",
    reference_fasta="reference.fa",
    output_prefix="output",

    # ROC Analysis Configuration (Phase 2)
    enable_roc_analysis=True,           # Enable ROC analysis
    roc_bootstrap_samples=1000,         # Bootstrap samples for confidence intervals
    quality_stratification=True,        # Enable quality stratification
)

# Run the analysis
results = engine.quantify()
```

### Command Line Options

For command-line usage, ROC analysis is enabled by default. You can control it through:

```bash
# ROC analysis is enabled by default in hap.py
hap.py truth.vcf query.vcf -r reference.fa -o output

# For direct quantify usage
qfy.py truth.vcf query.vcf -r reference.fa -o output
```

## Examples

### Example 1: Basic ROC Analysis

```bash
# Run hap.py with ROC analysis (enabled by default)
hap.py truth.vcf query.vcf -r reference.fa -o benchmark_results

# This generates:
# - benchmark_results.summary.csv
# - benchmark_results.roc.tsv
# - benchmark_results.quality_stratification.tsv
# - benchmark_results.multi_threshold.tsv
```

### Example 2: Programmatic Usage

```python
from hap_py.haplo.python_quantify import QuantifyEngine

# Create quantify engine with custom configuration
engine = QuantifyEngine(
    truth_vcf="truth.vcf.gz",
    query_vcf="query.vcf.gz",
    reference_fasta="hg38.fa",
    output_prefix="analysis_output",
    enable_roc_analysis=True,
    roc_bootstrap_samples=2000,  # More bootstrap samples for higher confidence
    quality_stratification=True
)

# Run analysis
results = engine.quantify()

# Access ROC data
roc_data = results.get('roc_data', {})
confidence_intervals = results.get('bootstrap_confidence_intervals', {})
quality_metrics = results.get('quality_metrics', {})

print(f"SNP precision at Q30: {roc_data['snp']['precision'][threshold_index]}")
```

### Example 3: Quality Threshold Selection

Use the multi-threshold analysis to select appropriate quality cutoffs:

```python
# Load multi-threshold results
import pandas as pd

thresholds = pd.read_csv("output.multi_threshold.tsv", sep="\t")

# Find optimal threshold for SNPs (maximize F1 score)
snp_thresholds = thresholds[thresholds['variant_type'] == 'snp']
optimal_threshold = snp_thresholds.loc[snp_thresholds['f1_score'].idxmax()]

print(f"Optimal SNP threshold: Q{optimal_threshold['threshold']}")
print(f"F1 score: {optimal_threshold['f1_score']:.4f}")
print(f"Precision: {optimal_threshold['precision']:.4f}")
print(f"Recall: {optimal_threshold['recall']:.4f}")
```

## Interpreting Results

### Understanding ROC Curves

1. **Precision-Recall Trade-off**: Higher quality thresholds typically increase precision but may decrease recall
2. **AUC (Area Under Curve)**: Higher AUC indicates better overall performance
3. **Confidence Intervals**: Wider intervals indicate less reliable estimates

### Quality Stratification Insights

- **Low Quality Bins (Q1-10)**: Often show lower precision, useful for identifying problematic calls
- **High Quality Bins (Q40+)**: Should show high precision, validates quality score calibration
- **Distribution Analysis**: Understanding how variants distribute across quality bins

### Best Practices

1. **Consider Both Precision and Recall**: Don't optimize for one metric alone
2. **Use Confidence Intervals**: Report uncertainty in your measurements
3. **Stratify by Variant Type**: SNPs and INDELs often have different performance characteristics
4. **Quality Score Calibration**: Use stratification to assess if quality scores are well-calibrated

## Advanced Usage

### Custom Quality Bins

While the default quality bins work for most cases, you may want to customize them:

```python
# The quality bins are currently hardcoded but could be made configurable in future versions
# Current bins: Q1-10, Q10-20, Q20-30, Q30-40, Q40+
```

### Bootstrap Configuration

Adjust bootstrap sampling for different confidence levels:

```python
# More bootstrap samples = more accurate confidence intervals but slower computation
engine = QuantifyEngine(
    # ... other parameters ...
    roc_bootstrap_samples=5000,  # Higher for publication-quality confidence intervals
)
```

### Performance Considerations

For large datasets:

1. **Memory Usage**: ROC analysis stores intermediate results; monitor memory usage
2. **Computation Time**: Bootstrap sampling can be time-intensive for large datasets
3. **Disk Space**: ROC output files can be large for datasets with many quality thresholds

### Integration with Other Tools

The ROC analysis output files are designed to be compatible with:

- **R/ggplot2**: For custom visualization
- **Python/matplotlib**: For programmatic plotting
- **Excel/LibreOffice**: For manual analysis and reporting

## Troubleshooting

### Common Issues

1. **Missing Dependencies**: Install optional dependencies for full functionality:
   ```bash
   pip install scipy matplotlib scikit-learn
   ```

2. **Memory Issues**: For very large datasets, consider:
   - Reducing bootstrap samples
   - Processing subsets of data
   - Using quality filtering

3. **Empty ROC Data**: Check that:
   - Quality scores are present in VCF files
   - Variants are being matched correctly
   - Output directory is writable

### Performance Tips

1. **Bootstrap Samples**: Start with 1000 samples, increase for publication-quality results
2. **Quality Filtering**: Pre-filter very low quality variants if not needed
3. **Variant Type Focus**: Disable analysis for variant types not of interest

## Future Enhancements

The quantify module continues to evolve. Planned enhancements include:

- **Phase 3**: Superlocus analysis for more accurate complex variant handling
- **Phase 4**: Performance optimizations for large-scale genomic datasets
- **Phase 5**: GA4GH compliance and standardized benchmarking protocols

## Implementation Phases

The quantify module has been modernized through a structured phase implementation approach to ensure all functionality from the original C++ implementation is preserved while adding enhanced features.

### Phase 1: Core Variant Matching (✅ COMPLETED)

**Status**: Phase 1 is complete with all core variant matching functionality working and tested.

#### Implementation Highlights
- **Sophisticated Variant Matching**: The `_match_variants` method handles both pandas Series and dictionary inputs
- **Multi-allelic Support**: Enhanced compatibility for complex variant representations
- **Performance Testing**: Realistic performance expectations established for the Python implementation
- **Comprehensive Testing**: Robust tests for variant matching algorithms added

#### Key Methods Implemented
- `_match_variants()`: Core variant matching algorithm with sophisticated allele compatibility
- Benchmarking decision tracking for BD, BVT, QQ field handling
- Variant normalization and comparison logic for consistent matching
- Support for both XCMP and GA4GH quantification methods

### Phase 2: Enhanced ROC Analysis (✅ COMPLETED)

**Status**: Phase 2 is complete with all ROC analysis functionality implemented and tested.

#### Implementation Highlights
- **Advanced Statistical Analysis**: Bootstrap confidence intervals using Jeffreys method
- **Comprehensive ROC Curves**: ROC generation for SNPs, INDELs, and all variants
- **Quality Stratification**: Performance analysis across quality score ranges
- **Multi-threshold Analysis**: Standard quality thresholds (Q10, Q20, Q30, Q40, Q50)

#### Key Methods Implemented
- `_perform_roc_analysis()`: Main orchestrator for complete ROC workflow
- `_generate_roc_curve()`: ROC curve generation with confidence intervals
- `_calculate_bootstrap_confidence_intervals()`: Statistical confidence calculation
- `_perform_quality_stratification()`: Quality score binning and analysis
- `_perform_multi_threshold_analysis()`: Standard threshold evaluation
- `_write_roc_results()`: Comprehensive output generation

### Phase 3: Superlocus Analysis (✅ COMPLETED)

**Status**: Phase 3 is complete with all superlocus and region-based functionality working and tested.

#### Implementation Highlights
- **Superlocus Analysis**: Sophisticated algorithms for overlapping and complex genomic regions
- **Region-based Quantification**: BED file integration for stratified analysis
- **Multi-sample Support**: Comparative analysis across multiple samples
- **Integration Testing**: Complete validation of Phase 3 functionality

#### Key Components Implemented
- `RegionBasedQuantifier`: BED file region filtering and intersection logic
- `MultiSampleQuantifier`: Multi-sample comparative analysis capabilities
- Superlocus identification algorithms for complex genomic regions
- Per-region statistics and stratified analysis support

### Phase 4: Performance Optimization (🔄 DEFERRED)

**Status**: Deferred until after full functionality validation. Will be revisited after Phase 5.

#### Planned Optimizations
- Memory optimization for large-scale datasets
- Processing speed enhancements through vectorization
- Streaming VCF processing for memory efficiency
- Multiprocessing for parallel operations

### Phase 5: GA4GH Compliance (✅ COMPLETED)

**Status**: Phase 5 is complete with full GA4GH compliance implemented and documented.

#### Implementation Highlights
- **Comprehensive GA4GH Support**: Complete implementation of GA4GH benchmarking standards
- **Standard Compliance Classes**: GA4GHFormatter, GA4GHStratification, and GA4GHMetrics
- **QuantifyEngine Integration**: Seamless integration with existing functionality
- **Complete Testing Suite**: Unit and integration tests for all GA4GH functionality

#### Key Components Implemented
- GA4GH VCF formatting with standard decision tracking
- Stratification region support for benchmarking standards
- Metrics calculation with confidence intervals
- Integration module for seamless QuantifyEngine compatibility

For detailed GA4GH implementation information, see the [GA4GH Compliance Guide](ga4gh_compliance.md).

## Current Implementation Status

### ✅ Completed Features
- **Phase 1**: Core variant matching algorithms with sophisticated allele compatibility
- **Phase 2**: Advanced ROC analysis with bootstrap confidence intervals
- **Phase 3**: Superlocus analysis and region-based quantification
- **Phase 5**: Complete GA4GH compliance implementation

### 🔄 Future Work
- **Phase 4**: Performance optimization for large-scale datasets (planned for future releases)
- Enhanced parallelization and memory optimization
- Additional benchmarking and validation against original C++ implementation

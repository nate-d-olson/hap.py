# ROC Analysis User Guide - hap.py Quantify Module

This guide provides practical instructions for using the enhanced ROC analysis features in hap.py's quantify module, implemented in Phase 2 of the modernization project.

## Quick Start

### Basic ROC Analysis

The simplest way to get ROC analysis is through the standard hap.py workflow:

```bash
# ROC analysis is enabled by default
hap.py truth.vcf query.vcf -r reference.fa -o my_analysis

# Check the generated ROC files
ls my_analysis.*
# my_analysis.summary.csv
# my_analysis.roc.tsv                    ← ROC curve data with confidence intervals
# my_analysis.quality_stratification.tsv ← Performance by quality bins
# my_analysis.multi_threshold.tsv        ← Standard threshold analysis
# my_analysis.roc_plot.png              ← Visualization (if matplotlib available)
```

### Direct Quantify Usage

For standalone quantification (requires preprocessed VCF files):

```bash
# Use qfy.py directly
qfy.py truth.vcf query.vcf -r reference.fa -o quantify_results
```

## Understanding the Output Files

### 1. ROC Curve Data (`*.roc.tsv`)

This file contains the precision-recall curve data with bootstrap confidence intervals:

```tsv
Type    Threshold    TP    FP    FN    Precision    Recall    Precision_Lower    Precision_Upper    Recall_Lower    Recall_Upper
SNP     0.00         1500  75    25    0.9524       0.9836    0.9312             0.9693            0.9712          0.9921
SNP     10.00        1420  45    105   0.9693       0.9311    0.9534             0.9834            0.9156          0.9466
SNP     20.00        1280  25    245   0.9807       0.8393    0.9656             0.9921            0.8187          0.8599
```

**Key Columns:**
- **Type**: SNP, INDEL, or ALL variants
- **Threshold**: Quality score cutoff
- **TP/FP/FN**: True positives, false positives, false negatives
- **Precision/Recall**: Performance metrics
- **\*_Lower/\*_Upper**: 95% confidence interval bounds

### 2. Quality Stratification (`*.quality_stratification.tsv`)

Shows performance within quality score bins:

```tsv
Quality_Bin    Quality_Range    TP    FP    FN    Precision    Recall    F1       Variant_Count
Q1_10          1-10             45    25    12    0.6429       0.7895    0.7088   82
Q10_20         10-20            185   15    18    0.9250       0.9113    0.9181   218
Q20_30         20-30            420   8     12    0.9813       0.9722    0.9767   440
Q30_40         30-40            680   3     8     0.9956       0.9884    0.9920   691
Q40_plus       40+              170   1     2     0.9942       0.9884    0.9913   173
```

**Use Cases:**
- **Quality Score Calibration**: Check if higher quality scores actually correspond to better performance
- **Threshold Selection**: Identify quality ranges with acceptable performance trade-offs
- **Dataset Characterization**: Understand the quality distribution of your variant calls

### 3. Multi-threshold Analysis (`*.multi_threshold.tsv`)

Performance at standardized quality thresholds for consistent benchmarking:

```tsv
Type     Threshold_Name    Threshold    TP     FP    FN    Precision    Recall    F1       Variants_Above_Threshold
SNP      Q10               10.0         1420   45    105   0.9693       0.9311    0.9498   1465
SNP      Q20               20.0         1280   25    245   0.9807       0.8393    0.9043   1305
SNP      Q30               30.0         1150   15    375   0.9871       0.7541    0.8564   1165
SNP      Q40               40.0         920    8     605   0.9914       0.6033    0.7487   928
SNP      Q50               50.0         680    3     845   0.9956       0.4459    0.6158   683
```

**Use Cases:**
- **Standardized Benchmarking**: Compare results across different studies using common thresholds
- **Filter Impact Analysis**: See how quality filtering affects your variant call set
- **Publication Reporting**: Standard thresholds provide comparable metrics for papers

## Practical Examples

### Example 1: Finding the Optimal Quality Threshold

```bash
# Run analysis
hap.py truth.vcf query.vcf -r ref.fa -o analysis

# Examine multi-threshold results to find optimal F1 score
awk 'NR==1 || $1=="SNP"' analysis.multi_threshold.tsv | sort -k9 -nr | head -5
```

This will show the top 5 SNP thresholds ranked by F1 score.

### Example 2: Quality Score Calibration Check

```bash
# Look at quality stratification to check calibration
cat analysis.quality_stratification.tsv
```

Well-calibrated quality scores should show:
- Q1-10: Lower precision/recall
- Q40+: Higher precision/recall
- Smooth progression between bins

### Example 3: Comparing SNP vs INDEL Performance

```bash
# Extract SNP ROC data
awk '$1=="SNP"' analysis.roc.tsv > snp_roc.tsv

# Extract INDEL ROC data
awk '$1=="INDEL"' analysis.roc.tsv > indel_roc.tsv

# Compare precision at Q30 threshold
awk '$2==30.0 {print $1, $6}' analysis.roc.tsv
```

### Example 4: Confidence Interval Analysis

Look for cases where confidence intervals are wide (indicating uncertainty):

```bash
# Find ROC points with wide precision confidence intervals (>0.05 width)
awk 'NR>1 && ($9-$8)>0.05 {print $1, $2, $6, $8, $9}' analysis.roc.tsv
```

Wide intervals suggest:
- Small sample sizes at that threshold
- Need for more bootstrap samples
- Inherent uncertainty in the measurement

## Interpreting Results

### Understanding Precision-Recall Trade-offs

**High Precision, Lower Recall**: Conservative calling - fewer false positives but missing some true variants
```
Quality Q40: Precision=0.995, Recall=0.603
```

**Lower Precision, High Recall**: Sensitive calling - catching most true variants but with more false positives
```
Quality Q10: Precision=0.969, Recall=0.931
```

### Quality Score Assessment

**Well-Calibrated Quality Scores**: Performance improves smoothly with increasing quality
```
Q10-20: Precision=0.925, Recall=0.911
Q20-30: Precision=0.981, Recall=0.972
Q30-40: Precision=0.996, Recall=0.988
```

**Poorly-Calibrated Quality Scores**: Little performance difference across quality ranges
```
Q10-20: Precision=0.923, Recall=0.915
Q20-30: Precision=0.924, Recall=0.916  ← No improvement
Q30-40: Precision=0.925, Recall=0.917  ← No improvement
```

### Confidence Interval Interpretation

**Narrow Intervals**: High confidence in measurements
```
Precision=0.952 [0.943, 0.961]  ← Width = 0.018, good confidence
```

**Wide Intervals**: Lower confidence, potentially due to small sample size
```
Precision=0.952 [0.891, 0.987]  ← Width = 0.096, less reliable
```

## Advanced Usage

### Programmatic Access

```python
from hap_py.haplo.python_quantify import QuantifyEngine
import pandas as pd

# Run analysis
engine = QuantifyEngine(
    truth_vcf="truth.vcf",
    query_vcf="query.vcf",
    reference_fasta="ref.fa",
    enable_roc_analysis=True,
    roc_bootstrap_samples=2000  # More samples for tighter intervals
)

results = engine.quantify()

# Access ROC data directly
roc_data = results['roc_data']
confidence_intervals = results['bootstrap_confidence_intervals']

# Find optimal SNP threshold
snp_thresholds = pd.DataFrame({
    'threshold': roc_data['snp']['thresholds'],
    'precision': roc_data['snp']['precision'],
    'recall': roc_data['snp']['recall']
})

# Calculate F1 scores
snp_thresholds['f1'] = 2 * (snp_thresholds['precision'] * snp_thresholds['recall']) / (snp_thresholds['precision'] + snp_thresholds['recall'])

# Find optimal F1
optimal_idx = snp_thresholds['f1'].idxmax()
optimal_threshold = snp_thresholds.loc[optimal_idx, 'threshold']
print(f"Optimal SNP threshold: {optimal_threshold}")
```

### Custom Analysis Scripts

```python
# Load ROC data for custom analysis
import pandas as pd
import matplotlib.pyplot as plt

# Load ROC curve data
roc_data = pd.read_csv("analysis.roc.tsv", sep="\t")

# Plot precision-recall curves
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

for variant_type in ['SNP', 'INDEL']:
    data = roc_data[roc_data['Type'] == variant_type]
    ax.plot(data['Recall'], data['Precision'], label=variant_type, marker='o')

    # Add confidence bands (optional)
    ax.fill_between(data['Recall'],
                   data['Precision_Lower'],
                   data['Precision_Upper'],
                   alpha=0.2)

ax.set_xlabel('Recall')
ax.set_ylabel('Precision')
ax.set_title('Precision-Recall Curves with Confidence Intervals')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('custom_roc_plot.png', dpi=300)
```

## Troubleshooting

### Common Issues

**1. Missing ROC Output Files**
```bash
# Check if ROC analysis is enabled (should be by default)
grep -i "roc analysis" hap.py_log.txt
```

**2. Empty or Sparse ROC Data**
```bash
# Check if quality scores are present in VCF
bcftools query -f '%QUAL\n' query.vcf | head -20

# Check for variant matches
grep -c "TP\|FP" analysis.summary.csv
```

**3. Wide Confidence Intervals**
- Increase bootstrap samples: add `--roc-bootstrap-samples 5000` (if available via CLI)
- Check sample size: wide intervals may indicate small variant counts
- Consider filtering very low-quality variants if they dominate the data

**4. Memory Issues with Large Datasets**
- Reduce bootstrap samples to 500 for initial analysis
- Process chromosome subsets if working with whole genomes
- Monitor memory usage during analysis

### Performance Tips

**For Large Datasets:**
- Start with 1000 bootstrap samples, increase gradually
- Use quality pre-filtering if appropriate for your analysis
- Consider processing subsets for initial exploration

**For Publication-Quality Results:**
- Use 2000-5000 bootstrap samples for tighter confidence intervals
- Include both SNP and INDEL analysis
- Document the bootstrap sample size in methods

## Integration with Other Tools

### R/ggplot2 Analysis

```r
library(ggplot2)
library(dplyr)

# Load ROC data
roc_data <- read.table("analysis.roc.tsv", header=TRUE, sep="\t")

# Create precision-recall plot
ggplot(roc_data, aes(x=Recall, y=Precision, color=Type)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=Precision_Lower, ymax=Precision_Upper, fill=Type),
              alpha=0.2, color=NA) +
  theme_minimal() +
  labs(title="ROC Analysis with Confidence Intervals",
       x="Recall", y="Precision") +
  scale_color_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1")
```

### Excel/LibreOffice Analysis

The TSV files can be directly opened in spreadsheet applications:

1. Open `*.roc.tsv` in Excel
2. Create scatter plots with Recall (x-axis) and Precision (y-axis)
3. Add error bars using the confidence interval columns
4. Use conditional formatting to highlight optimal thresholds

## Best Practices

### Threshold Selection

1. **Consider Use Case**: High precision for clinical applications, balanced F1 for research
2. **Examine Confidence Intervals**: Choose thresholds with reliable estimates
3. **Validate on Independent Data**: Test selected thresholds on holdout datasets
4. **Document Methodology**: Report threshold selection criteria and confidence intervals

### Quality Assessment

1. **Check Calibration**: Quality stratification should show performance progression
2. **Compare Variant Types**: SNPs and INDELs often have different optimal thresholds
3. **Consider Context**: Population frequency, genomic regions, and calling method affect performance

### Reporting Results

1. **Include Confidence Intervals**: Report uncertainty in performance estimates
2. **Specify Bootstrap Samples**: Document the number of bootstrap samples used
3. **Show Multiple Thresholds**: Present results across quality ranges, not just single points
4. **Variant Type Stratification**: Report SNP and INDEL performance separately

This comprehensive guide should help users effectively utilize the enhanced ROC analysis capabilities in hap.py's quantify module.

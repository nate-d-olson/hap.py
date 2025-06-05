# GA4GH Compliance in hap.py

## Overview

The GA4GH (Global Alliance for Genomics and Health) has established benchmarking standards for variant calling evaluation. The hap.py tool provides full support for these standards through comprehensive implementation of:

1. GA4GH stratification standards
2. Benchmarking metrics according to GA4GH specifications  
3. Output formatting that adheres to GA4GH requirements
4. Complete VCF intermediate format support
5. Integration with RTG Tools for GA4GH-compliant comparison

## Implementation Architecture

### Core Components

The GA4GH compliance implementation consists of four main classes that work together to provide complete standards compliance:

#### GA4GHFormatter
Implements GA4GH-compliant VCF formatting and annotation:
- **VCF Header Management**: Adds GA4GH-required FORMAT and INFO fields
- **Record Annotation**: Annotates VCF records with GA4GH-specific values
- **Standard Field Support**: Implements BD, BK, QD, QK fields for decision tracking
- **Subtype Classification**: Provides SNP, INDEL, COMPLEX variant classification

Key methods:
- `format_vcf_header()`: Adds GA4GH fields to VCF header
- `annotate_record()`: Annotates records with decision values and details
- `add_ga4gh_info_fields()`: Adds standard INFO fields for regions and status

#### GA4GHStratification  
Handles GA4GH stratification regions for benchmarking:
- **BED File Integration**: Loads and manages stratification regions from BED files
- **Region Assignment**: Determines which regions contain each variant
- **Standard Stratifications**: Supports high-confidence, difficult, and segmental duplication regions
- **Efficient Querying**: Uses pybedtools for fast region overlap detection

Key methods:
- `add_region()`: Adds a stratification region from BED file
- `get_region_ids_for_variant()`: Returns region IDs for a specific variant
- `stratify_variants()`: Groups variants by stratification regions

#### GA4GHMetrics
Calculates benchmarking metrics according to GA4GH standards:
- **Standard Metrics**: Precision, recall, and F1-score calculation
- **Confidence Intervals**: Bootstrap-based confidence interval calculation
- **Stratified Metrics**: Metrics calculation by region and variant type
- **GA4GH Output**: Produces GA4GH-compliant metrics files

Key methods:
- `calculate_precision()`: Calculates precision with optional confidence intervals
- `calculate_recall()`: Calculates recall with optional confidence intervals  
- `calculate_f1()`: Calculates F1-score with optional confidence intervals
- `calculate_metrics()`: Comprehensive metrics calculation

#### GA4GHIntegration
Provides seamless integration with QuantifyEngine:
- **Engine Integration**: Integrates GA4GH components with QuantifyEngine
- **Workflow Management**: Manages GA4GH-specific analysis workflows
- **Output Coordination**: Coordinates GA4GH output generation
- **Backward Compatibility**: Ensures compatibility with existing functionality

## GA4GH VCF Format Specification

### Required Header Fields

#### FORMAT Fields
```
##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision for the benchmark variant (TP/FP/FN/N/UNK)">
##FORMAT=<ID=BK,Number=1,Type=String,Description="Decision detail for the benchmark variant">  
##FORMAT=<ID=QD,Number=1,Type=String,Description="Decision for the query variant (TP/FP/FN/N/UNK)">
##FORMAT=<ID=QK,Number=1,Type=String,Description="Decision detail for the query variant">
```

#### INFO Fields
```
##INFO=<ID=Regions,Number=.,Type=String,Description="List of region IDs this variant is located in">
##INFO=<ID=TruthStatus,Number=1,Type=String,Description="Status of the variant in truth VCF (TP/FN/FP/N/UNK)">
##INFO=<ID=QueryStatus,Number=1,Type=String,Description="Status of the variant in query VCF (TP/FP/FN/N/UNK)">  
##INFO=<ID=Subtype,Number=1,Type=String,Description="Variant subtype classification (SNP/INDEL/COMPLEX/OTHER)">
```

### Standard Decision Values

The implementation uses standardized enums for consistent decision tracking:

#### GA4GHDecision
- **TP**: True positive
- **FP**: False positive  
- **FN**: False negative
- **N**: Non-assessed (variant in non-confident region)
- **UNK**: Unknown/undetermined

#### GA4GHDecisionDetail  
- **GT_MATCH**: Genotype match
- **GT_MISMATCH**: Genotype mismatch
- **ALLELE_MATCH**: Allele match only
- **ALLELE_MISMATCH**: Allele mismatch
- **OUTSIDE_CONFIDENT**: Outside confident regions
- **REFERENCE_MATCH**: Matches reference
- **COMPLEX_MATCH**: Complex representation match
- **NO_MATCH**: No match found

## Usage Examples

### Basic GA4GH Analysis

```python
from hap_py.haplo.python_quantify import QuantifyEngine

# Create QuantifyEngine with GA4GH compliance
engine = QuantifyEngine(
    truth_vcf="truth.vcf",
    query_vcf="query.vcf", 
    quantify_method="ga4gh",
    enable_ga4gh=True
)

# Run analysis with GA4GH output
results = engine.run()
```

### GA4GH Formatting Example

```python
from hap_py.haplo.ga4gh_compliance import GA4GHFormatter, GA4GHDecision, GA4GHDecisionDetail

# Create formatter
formatter = GA4GHFormatter()

# Format VCF header
formatter.format_vcf_header(vcf_header)

# Annotate a record
formatter.annotate_record(
    record,
    truth_decision=GA4GHDecision.TP,
    query_decision=GA4GHDecision.TP,
    truth_detail=GA4GHDecisionDetail.GT_MATCH,
    query_detail=GA4GHDecisionDetail.GT_MATCH
)
```

### Stratification with GA4GH Regions

```python
from hap_py.haplo.ga4gh_compliance import GA4GHStratification

# Create stratification with standard regions
stratification = GA4GHStratification({
    "highconf": "/path/to/highconf.bed",
    "difficult": "/path/to/difficult.bed", 
    "segdup": "/path/to/segmental_duplications.bed"
})

# Get regions for a variant
regions = stratification.get_region_ids_for_variant("chr1", 1000, 1001)

# Stratify all variants
stratified_results = stratification.stratify_variants(variants_df)
```

### Metrics Calculation with Confidence Intervals

```python
from hap_py.haplo.ga4gh_compliance import GA4GHMetrics

# Create metrics calculator with bootstrap confidence intervals
metrics = GA4GHMetrics(bootstrap_iterations=1000, ci_level=0.95)

# Calculate comprehensive metrics
result = metrics.calculate_metrics(tp=90, fp=10, fn=10, with_ci=True)

# Extract metrics with confidence intervals
precision, prec_lower, prec_upper = result["precision"]
recall, rec_lower, rec_upper = result["recall"] 
f1, f1_lower, f1_upper = result["f1"]

print(f"Precision: {precision:.3f} [{prec_lower:.3f}, {prec_upper:.3f}]")
print(f"Recall: {recall:.3f} [{rec_lower:.3f}, {rec_upper:.3f}]")
print(f"F1-score: {f1:.3f} [{f1_lower:.3f}, {f1_upper:.3f}]")
```

## Command-Line Usage

### Basic GA4GH Analysis
```bash
# Run hap.py with GA4GH compliance
hap.py truth.vcf query.vcf --engine=xcmp --quantify-method=ga4gh --output-format=ga4gh
```

### With Stratification Regions
```bash
# Include stratification regions
hap.py truth.vcf query.vcf \
  --engine=xcmp \
  --quantify-method=ga4gh \
  --ga4gh-stratification=regions.bed \
  --output-format=ga4gh
```

### RTG VCFEval Integration
```bash
# Use RTG vcfeval engine with GA4GH output
hap.py truth.vcf query.vcf \
  --engine=vcfeval \
  --quantify-method=ga4gh \
  --engine-vcfeval-path=/path/to/rtg
```
- Precision and recall with confidence intervals
- F-measure with confidence intervals
- False positive and false negative rates
- Quality metrics and ROC analysis

## Usage

To generate GA4GH-compliant outputs, use the `--ga4gh` flag:

```bash
hap.py truth.vcf query.vcf -o output --ga4gh
```

This will generate additional outputs:
- `output.ga4gh.json`: GA4GH-compliant metrics in JSON format
- `output.ga4gh.tsv`: GA4GH-compliant metrics in TSV format
- `output.ga4gh.extended.csv`: Detailed GA4GH metrics with stratifications

## Interoperability with External Tools

The GA4GH-compliant output from hap.py is designed to work seamlessly with other GA4GH-compliant tools:

### RTG Tools Integration
- **VCFEval Engine**: Native integration with RTG vcfeval for GA4GH-compliant comparison
- **SDF Format**: Support for SDF reference genome format required by RTG
- **Intermediate Files**: Compatible with RTG intermediate file formats

### Other GA4GH Tools
- **GA4GH Benchmarking Tools**: Direct compatibility with reference implementations
- **Truvari**: GA4GH output can be processed by Truvari for additional analysis
- **vcflib**: VCF files with GA4GH annotations work with vcflib utilities

## Testing and Validation

### Unit Tests
Comprehensive unit tests validate all GA4GH functionality:
- `tests/unit/test_ga4gh_compliance.py`: Core GA4GH class testing
- Individual class validation for all GA4GH components
- Decision value and annotation testing

### Integration Tests  
End-to-end testing ensures GA4GH workflows function correctly:
- `tests/integration/test_ga4gh_integration.py`: Complete workflow testing
- RTG integration testing with GA4GH output
- Stratification and metrics validation

### Validation Scripts
- `test_ga4gh_implementation.py`: Comprehensive validation of GA4GH implementation
- Verification against GA4GH reference standards
- Performance and compatibility testing

## Troubleshooting

### Common Issues

#### Missing pybedtools Dependency
```bash
# Install pybedtools for stratification support
conda install -c bioconda pybedtools
# or
pip install pybedtools
```

#### VCF Header Issues
Ensure VCF files have proper headers for GA4GH annotation:
```python
# Check if header has required fields before GA4GH annotation
if not formatter.has_ga4gh_fields(vcf_header):
    formatter.format_vcf_header(vcf_header)
```

#### Region File Format
Stratification BED files must follow standard format:
```
chr1    100000    200000    high_confidence
chr1    300000    400000    difficult_region
```

## Implementation Status

### ✅ Completed Features
- Core GA4GH classes (Formatter, Stratification, Metrics, Integration)
- VCF annotation with GA4GH fields
- Standard decision value tracking
- Confidence interval calculation
- QuantifyEngine integration
- Comprehensive testing suite
- Documentation and examples

### 🔄 Future Enhancements
- Performance optimization for large-scale analysis
- Additional stratification region types
- Enhanced metrics visualization
- Extended RTG integration features

## References

- [GA4GH Benchmarking Standards](https://github.com/ga4gh/benchmarking-tools)
- [GA4GH VCF Specification](https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/README.md)
- [RTG Tools Documentation](https://github.com/RealTimeGenomics/rtg-tools)

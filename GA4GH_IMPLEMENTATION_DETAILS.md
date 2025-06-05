# GA4GH Compliance Implementation Details

## Overview

This document details the implementation of GA4GH compliance in the modernized hap.py codebase. The implementation follows the [GA4GH benchmarking standards](https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/README.md) and ensures compatibility with the existing functionality.

## GA4GH Standards Relevant to hap.py

The GA4GH benchmarking standards define specific formats and requirements for variant comparison and benchmarking:

1. **Intermediate VCF Format**: A standardized format for representing benchmark results in VCF files
2. **Stratification Regions**: Standardized regions for benchmarking defined in BED files
3. **Metrics Calculation**: Standard methods for calculating benchmarking metrics

## Implementation Components

### 1. GA4GH Classes

The implementation consists of four main classes:

#### GA4GHFormatter

This class handles GA4GH-compliant VCF formatting:

- **FORMAT fields**: BD, BK, QD, QK fields for decision values and details
- **INFO fields**: Regions, TruthStatus, QueryStatus, Subtype
- **Methods**:
  - `format_vcf_header()`: Adds GA4GH fields to a VCF header
  - `annotate_record()`: Annotates a VCF record with GA4GH-specific values

#### GA4GHStratification

This class handles stratification regions:

- Manages region BED files for stratification
- **Methods**:
  - `add_region()`: Adds a stratification region
  - `get_region_ids_for_variant()`: Assigns variants to regions
  - `stratify_variants()`: Groups variants by stratification region

#### GA4GHMetrics

This class calculates benchmarking metrics:

- **Metrics**:
  - Precision (TP / (TP + FP))
  - Recall (TP / (TP + FN)) 
  - F1-score (2 * precision * recall / (precision + recall))
- **Confidence Intervals**: Provides methods to calculate confidence intervals for metrics
- **Methods**:
  - `calculate_precision()`: Calculates precision with optional confidence intervals
  - `calculate_recall()`: Calculates recall with optional confidence intervals
  - `calculate_f1()`: Calculates F1-score with optional confidence intervals
  - `calculate_metrics()`: Calculates all metrics with optional confidence intervals

### 2. Supporting Enums

The implementation includes several enums to standardize the values used in GA4GH annotations:

#### GA4GHDecision

Standardized decision values for variant classification:

- **TP**: True positive
- **FP**: False positive
- **FN**: False negative
- **N**: Non-assessed (variant in non-confident region)
- **UNK**: Unknown/undetermined

#### GA4GHDecisionDetail

Detailed decision information:

- **GT_MATCH**: Genotype match
- **GT_MISMATCH**: Genotype mismatch
- **ALLELE_MATCH**: Allele match only
- **ALLELE_MISMATCH**: Allele mismatch
- **OUTSIDE_CONFIDENT**: Outside confident regions
- **REFERENCE_MATCH**: Matches reference
- **COMPLEX_MATCH**: Complex representation match
- **NO_MATCH**: No match found

#### GA4GHVariantType

Variant type classification:

- **SNP**: Single nucleotide polymorphism
- **INDEL**: Insertion/deletion
- **COMPLEX**: Complex variation
- **OTHER**: Other variation type

### 3. Integration with QuantifyEngine

The GA4GH classes are integrated with the existing QuantifyEngine:

- **GA4GH Matching**: Enhanced `_perform_ga4gh_matching()` method in QuantifyEngine
- **GA4GH Output**: Support for GA4GH-compliant output formats
- **Configuration**: QuantifyEngine accepts `quantify_method="ga4gh"` to use GA4GH-compliant comparison

## GA4GH VCF Format Specification

### Header Fields

The GA4GH intermediate format requires specific header fields:

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

### Record Format

An example GA4GH-compliant VCF record:

```
chr1  10000  .  G  A  50  PASS  Regions=highconf,exome;TruthStatus=TP;QueryStatus=TP;Subtype=SNP  GT:BD:BK:QD:QK  0/1:TP:gt-match:TP:gt-match
```

## Usage Examples

### Basic Usage

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

### Stratification Example

```python
from hap_py.haplo.ga4gh_compliance import GA4GHStratification

# Create stratification
stratification = GA4GHStratification({
    "highconf": "/path/to/highconf.bed",
    "exome": "/path/to/exome.bed"
})

# Get regions for a variant
regions = stratification.get_region_ids_for_variant("chr1", 1000, 1001)

# Stratify variants
stratified = stratification.stratify_variants(variants_df)
```

### Metrics Example

```python
from hap_py.haplo.ga4gh_compliance import GA4GHMetrics

# Create metrics calculator
metrics = GA4GHMetrics(bootstrap_iterations=1000, ci_level=0.95)

# Calculate metrics with confidence intervals
result = metrics.calculate_metrics(tp=90, fp=10, fn=10, with_ci=True)

# Access metrics
precision, prec_lower, prec_upper = result["precision"]
recall, rec_lower, rec_upper = result["recall"]
f1, f1_lower, f1_upper = result["f1"]
```

## Command-Line Integration

To use GA4GH-compliant comparison mode:

```bash
hap.py truth.vcf query.vcf --engine=xcmp --quantify-method=ga4gh --output-format=ga4gh
```

## Compatibility with External Tools

The GA4GH-compliant output is compatible with:

1. RTG vcfeval
2. GA4GH benchmarking tools
3. Truvari
4. vcflib

### 1. Variant Representation and Benchmarking Standards

- **GA4GH Intermediate Format**: Standardized format for variant benchmarking output
- **VCF/BCF Extensions**: GA4GH-specific extensions to VCF format for variant comparison results
- **Stratification Standards**: Standardized regions for benchmarking comparisons

### 2. Benchmarking Framework Requirements

- **Methodology Standards**: Follow GA4GH benchmarking best practices
- **Metric Definitions**: Implement standardized precision, recall, and F-score calculations
- **Output Formats**: Generate GA4GH-compliant output files for interoperability

### 3. Integration with RTG Tools

- **VCFEval Integration**: Support for RTG Tools vcfeval as a GA4GH-compliant comparison engine
- **SDF Reference Format**: Support for SDF reference genome format
- **Intermediate File Support**: Compatible with GA4GH intermediate file formats

## Implementation Components

### Phase 5.1: GA4GH Core Classes (Completed)

1. **GA4GHFormatter**
   - Implements VCF header modifications for GA4GH compliance
   - Adds GA4GH-required FORMAT fields (BD, BK, QD, QK)
   - Adds GA4GH-required INFO fields (Regions, TruthStatus, QueryStatus, Subtype)

2. **GA4GHStratification**
   - Handles GA4GH stratification regions
   - Loads BED files for standard stratification
   - Determines region overlaps for variants

3. **GA4GHMetrics**
   - Calculates metrics according to GA4GH standards
   - Supports stratification by region and variant type
   - Produces GA4GH-compliant output formats

### Phase 5.2: Integration with QuantifyEngine (In Progress)

1. **GA4GH Mode in QuantifyEngine**
   - Add GA4GH-specific parameters to QuantifyEngine
   - Integrate GA4GH formatter for output generation
   - Support GA4GH stratification regions

2. **RTG Tools Integration**
   - Ensure compatibility with RTG Tools vcfeval engine
   - Support GA4GH intermediate file format
   - Handle SDF reference requirements

### Phase 5.3: Documentation and Testing (Pending)

1. **GA4GH Compliance Documentation**
   - Document GA4GH-specific functionality
   - Provide examples of GA4GH workflows
   - Update CLI documentation for GA4GH options

2. **Comprehensive Testing**
   - Unit tests for GA4GH classes
   - Integration tests for GA4GH workflows
   - Validation against reference implementations

## Technical Details

### GA4GH Format Requirements

```
# FORMAT fields
BD: Decision for baseline (truth) genotype (TP/FP/FN/N/UNK)
BK: Baseline (truth) genotype filter status
QD: Decision for query (test) genotype (TP/FP/FN/N/UNK)
QK: Query (test) genotype filter status

# INFO fields
Regions: Stratification regions
TruthStatus: Status in truth VCF (TP/FP/FN)
QueryStatus: Status in query VCF (TP/FP/FN)
Subtype: Variant subtype classification
```

### Integration with Existing Code

The GA4GH compliance implementation is designed to be optional and non-disruptive, with the following integration points:

1. **QuantifyEngine Integration**
   - New parameter: `enable_ga4gh=True/False`
   - New parameter: `ga4gh_stratification_regions={...}`

2. **Command-Line Interface**
   - New option: `--ga4gh-output`
   - New option: `--ga4gh-stratification=<BED file>`

3. **Output Files**
   - Annotated VCF file with GA4GH fields
   - GA4GH-compliant metrics report
   - Standard ROC curves with GA4GH format

## Testing Strategy

### 1. Unit Tests

- Tests for GA4GHFormatter class
- Tests for GA4GHStratification class
- Tests for GA4GHMetrics class

### 2. Integration Tests

- End-to-end workflow tests with GA4GH output
- Validation against reference implementations
- Performance benchmarks

### 3. Validation Tests

- Compliance with GA4GH specifications
- Interoperability with other GA4GH tools
- Consistency with existing hap.py functionality

## Timeline and Milestones

### Phase 5.1: GA4GH Core Classes (Completed)

- [x] GA4GHFormatter implementation
- [x] GA4GHStratification implementation
- [x] GA4GHMetrics implementation
- [x] Unit tests for core classes

### Phase 5.2: Integration with QuantifyEngine (In Progress)

- [ ] Add GA4GH mode to QuantifyEngine
- [ ] Integrate with RTG Tools vcfeval
- [ ] Update command-line interface
- [ ] Integration tests

### Phase 5.3: Documentation and Testing (Pending)

- [ ] Document GA4GH functionality
- [ ] Create examples and tutorials
- [ ] Validate against reference implementations
- [ ] Performance optimization

## Conclusion

The GA4GH compliance implementation enhances hap.py's interoperability with the broader genomics ecosystem while maintaining backward compatibility with existing workflows. The modular approach ensures that GA4GH functionality can be enabled or disabled as needed, and the comprehensive testing strategy ensures compliance with GA4GH standards.

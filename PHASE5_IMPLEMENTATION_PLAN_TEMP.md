# Phase 5 Implementation Plan: GA4GH Compliance and Standards Support

## Overview

Phase 5 focuses on ensuring the quantify module in the modernized hap.py codebase complies with Global Alliance for Genomics and Health (GA4GH) standards. For additional information about the standard see the [documentation](https://github.com/Illumina/hap.py/blob/master/doc/happy.md) from the original implementation and the code from the repository to ensure the GA4GH benchmark intermediate vcf specification (as defined [here](https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/intermediate.md)) is implemented correctly. Here is the [repository](https://github.com/ga4gh/benchmarking-tools/tree/master) with the benchmarking standard official definition for reference. This implementation plan outlines the necessary steps to achieve GA4GH compliance while maintaining compatibility with the existing codebase.

## GA4GH Standards Relevant to hap.py

### 1. Variant Representation and Benchmarking Standards

- **GA4GH Intermediate Format**: Standardized format for variant benchmarking output as defined in the [GA4GH benchmarking documentation](https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/README.md)
- **VCF/BCF Extensions**: GA4GH-specific extensions to VCF format for variant comparison results, including specialized annotations for TP/FP/FN/N/UNK classifications
- **Stratification Standards**: Standardized regions for benchmarking comparisons following the GA4GH stratification BED files available at [GA4GH stratification resources](https://github.com/ga4gh/benchmarking-tools/tree/master/resources/stratification-bed-files)

### 2. Benchmarking Framework Requirements

- **Methodology Standards**: Follow GA4GH benchmarking best practices
- **Metric Definitions**: Implement standardized precision, recall, and F-score calculations
- **Output Formats**: Generate GA4GH-compliant output files for interoperability

### 3. Integration with RTG Tools

- **VCFEval Integration**: Support for RTG Tools vcfeval as a GA4GH-compliant comparison engine
- **SDF Reference Format**: Support for SDF reference genome format
- **Intermediate File Support**: Compatible with GA4GH intermediate file formats

## Implementation Components

### Phase 5.1: GA4GH Core Implementation (Weeks 1-2)

1. **GA4GH VCF Format Support**
   - Implement GA4GH-compliant VCF output format with required annotations
   - Add support for TP/FP/FN/N/UNK classifications in VCF records
   - Create converters between internal formats and GA4GH formats

2. **GA4GH Stratification Integration**
   - Import GA4GH standard stratification regions
   - Support region-based performance metrics
   - Implement proper region labeling in output files

3. **RTG Tools VCFEval Integration**
   - Ensure proper use of RTG Tools for GA4GH-compliant comparisons
   - Update VCFEval runners to output GA4GH intermediate format
   - Handle SDF reference requirements

### Phase 5.2: Metrics and Output Formatting (Weeks 3-4)

1. **GA4GH Metric Implementation**
   - Implement GA4GH-defined precision, recall, F-score calculations
   - Add confidence interval calculations per GA4GH specifications
   - Create stratified metrics across variant types and regions

2. **Output Format Standardization**
   - Create GA4GH-compliant summary files
   - Implement GA4GH extended CSV/TSV formats
   - Support GA4GH JSON output options

3. **ROC Analysis Enhancement**
   - Update ROC curve generation to follow GA4GH formats
   - Implement proper quality score binning per GA4GH standards
   - Add visualization options for GA4GH-compliant output

### Phase 5.3: Integration and Interoperability (Weeks 5-6)

1. **QuantifyEngine Integration**
   - Update the QuantifyEngine to support GA4GH modes
   - Add GA4GH-specific parameters to the engine
   - Create a GA4GH-specific workflow option

2. **CLI Updates**
   - Add GA4GH-specific command-line options
   - Provide presets for common GA4GH workflows
   - Update documentation for GA4GH usage

3. **Interoperability Testing**
   - Validate outputs against other GA4GH implementations
   - Test with GA4GH benchmark datasets
   - Create validators for GA4GH compliance

## Technical Details

### GA4GH Format Specifications

```python
# GA4GH VCF Header Requirements
GA4GH_FORMAT_FIELDS = [
    "BD", # Decision for baseline (truth) genotype (TP/FP/FN/N/UNK)
    "BK", # Baseline (truth) genotype filter status
    "QD", # Decision for query (test) genotype (TP/FP/FN/N/UNK)
    "QK"  # Query (test) genotype filter status
]

# GA4GH INFO fields
GA4GH_INFO_FIELDS = [
    "Regions", # Stratification regions
    "TruthStatus", # Status in truth VCF (TP/FP/FN)
    "QueryStatus", # Status in query VCF (TP/FP/FN)
    "Subtype"  # Variant subtype classification
]
```

### Implementation Classes

```python
class GA4GHFormatter:
    """
    Handles GA4GH-compliant output formatting.
    """
    def __init__(self, enable_ga4gh=True):
        self.enable_ga4gh = enable_ga4gh
        self.required_fields = GA4GH_FORMAT_FIELDS
        self.required_info = GA4GH_INFO_FIELDS

    def format_vcf_header(self, header):
        """Add GA4GH-specific header fields to VCF"""
        # Add required FORMAT fields
        for field in self.required_fields:
            if field not in header.formats:
                if field == "BD" or field == "QD":
                    header.add_format(field, "1", "String", "Decision (TP/FP/FN/N/UNK)")
                elif field == "BK" or field == "QK":
                    header.add_format(field, "1", "String", "Filter status")

        # Add required INFO fields
        for field in self.required_info:
            if field not in header.info:
                if field == "Regions":
                    header.add_info(field, ".", "String", "Stratification regions")
                elif field in ["TruthStatus", "QueryStatus"]:
                    header.add_info(field, "1", "String", "Variant status")
                elif field == "Subtype":
                    header.add_info(field, "1", "String", "Variant subtype")

        return header

    def annotate_record(self, record, truth_decision, query_decision, regions=None):
        """Annotate a VCF record with GA4GH-compliant annotations"""
        # Add FORMAT fields
        record.samples["TRUTH"]["BD"] = truth_decision
        record.samples["QUERY"]["QD"] = query_decision

        # Add INFO fields
        if regions:
            record.info["Regions"] = regions
        record.info["TruthStatus"] = truth_decision
        record.info["QueryStatus"] = query_decision

        return record
```

### Integration with QuantifyEngine

```python
class GA4GHQuantifyEngine(QuantifyEngine):
    """
    GA4GH-compliant extension of the QuantifyEngine.
    """
    def __init__(self, truth_vcf, query_vcf, **kwargs):
        # Initialize with GA4GH options
        super().__init__(
            truth_vcf=truth_vcf,
            query_vcf=query_vcf,
            **kwargs
        )
        self.ga4gh_formatter = GA4GHFormatter(
            enable_ga4gh=kwargs.get("enable_ga4gh", True)
        )
        self.ga4gh_regions = kwargs.get("ga4gh_stratification_regions", {})

    def quantify(self):
        """Run quantification with GA4GH compliance"""
        results = super().quantify()

        # Add GA4GH-specific processing
        if self.enable_ga4gh:
            self._format_ga4gh_output(results)

        return results

    def _format_ga4gh_output(self, results):
        """Format results according to GA4GH standards"""
        # Format VCF header
        if hasattr(results, "vcf_header"):
            results.vcf_header = self.ga4gh_formatter.format_vcf_header(results.vcf_header)

        # Annotate records
        if hasattr(results, "records"):
            for record in results.records:
                truth_decision = self._get_truth_decision(record)
                query_decision = self._get_query_decision(record)
                regions = self._get_stratification_regions(record)
                self.ga4gh_formatter.annotate_record(
                    record,
                    truth_decision,
                    query_decision,
                    regions
                )
```

## Testing Strategy

### 1. Unit Tests

- **GA4GH Format Tests**
   - Test proper GA4GH VCF header generation
   - Validate GA4GH annotations on records
   - Test GA4GH metric calculations

- **Integration Tests**
   - Verify GA4GH stratification with standard regions
   - Test end-to-end GA4GH workflow with RTG tools
   - Validate GA4GH output files against specifications

- **Conversion Tests**
   - Test conversion between internal formats and GA4GH formats
   - Verify consistency between GA4GH and non-GA4GH results
   - Test backward compatibility with existing workflows

### 2. Validation Tests

- **Reference Implementation Comparison**
   - Compare results against GA4GH reference implementations
   - Validate using GA4GH standard benchmarking datasets
   - Test interoperability with other GA4GH tools

- **GA4GH Specification Compliance**
   - Validate output files against GA4GH schemas
   - Check for required fields and annotations
   - Verify metrics follow GA4GH definitions

### 3. Edge Cases and Error Handling

- **Incomplete/Invalid Input Handling**
   - Test with missing or invalid GA4GH information
   - Verify graceful degradation when GA4GH requirements not met
   - Test error reporting for GA4GH compliance issues

## Documentation Requirements

### 1. GA4GH-Specific Documentation

- **GA4GH Standards Reference**
   - Document supported GA4GH standards and versions
   - Explain differences from reference implementations
   - Provide links to official GA4GH documentation

- **Usage Guide**
   - Create step-by-step guide for GA4GH-compliant analysis
   - Document GA4GH-specific command-line options
   - Provide examples of GA4GH output interpretation

### 2. Integration Documentation

- **Interoperability Guide**
   - Document how to use outputs with other GA4GH tools
   - Explain integration with RTG Tools
   - Provide conversion utilities documentation

- **Stratification Guide**
   - Document standard GA4GH stratification regions
   - Explain region-based metrics and reports
   - Provide guide for custom region definitions

## Timeline and Milestones

### Phase 5.1: GA4GH Core Implementation (2 weeks)
- Week 1: GA4GH VCF format support and RTG tools integration
- Week 2: GA4GH stratification implementation and testing

### Phase 5.2: Metrics and Output Formatting (2 weeks)
- Week 3: GA4GH metric implementation and validation
- Week 4: Output format standardization and ROC analysis updates

### Phase 5.3: Integration and Interoperability (2 weeks)
- Week 5: QuantifyEngine integration and CLI updates
- Week 6: Interoperability testing and documentation

## Success Criteria

- **Core Functionality**
  - [ ] GA4GH VCF output with required annotations
  - [ ] Standard stratification regions support
  - [ ] GA4GH-compliant metrics calculation

- **Integration**
  - [ ] RTG Tools VCFEval integration complete
  - [ ] QuantifyEngine with GA4GH options
  - [ ] CLI support for GA4GH workflows

- **Validation**
  - [ ] Output validated against GA4GH specifications
  - [ ] Interoperability with at least 2 GA4GH tools
  - [ ] Documentation complete for GA4GH features

## Resources Needed

- **Reference Documentation**
  - [GA4GH Benchmarking Framework](https://github.com/ga4gh/benchmarking-tools)
  - [GA4GH Intermediate VCF Specification](https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/intermediate.md)
  - [hap.py GA4GH Documentation](https://github.com/Illumina/hap.py/blob/master/doc/happy.md)

- **Tools and Libraries**
  - RTG Tools with GA4GH support
  - GA4GH validation tools
  - Standard benchmarking datasets

## Risk Management

### 1. Standard Compliance Challenges

- **Risk**: GA4GH standards may evolve during implementation
- **Impact**: Medium to High
- **Mitigation**: Regular verification against current specifications, modular design to accommodate changes

### 2. Performance Considerations

- **Risk**: GA4GH compliance may add overhead to processing
- **Impact**: Low to Medium
- **Mitigation**: Optional GA4GH features, performance optimization for core functions

### 3. Compatibility Issues

- **Risk**: Changes may affect existing workflows
- **Impact**: Medium
- **Mitigation**: Maintain backward compatibility, thorough testing of non-GA4GH paths

## Conclusion

By implementing GA4GH compliance in the hap.py quantify module, we ensure that the tool remains compatible with the broader genomics ecosystem and adheres to industry standards. The phased implementation approach allows us to maintain backward compatibility while incrementally adding GA4GH features, ensuring a smooth transition for existing users and providing enhanced capabilities for new users.

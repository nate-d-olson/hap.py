# Removed Components in the Modernized hap.py

During the Python 3 modernization of hap.py, some components from the original implementation were not migrated to preserve project scope and focus on the most critical functionality. This document details these removed components and how you can request their restoration.

## Removed Components

### 1. som.py - Somatic Comparison Tool

**Description:** som.py was a simple comparison tool based on bcftools for comparing somatic variants. It performed comparison based on variant location and alleles only, without considering genotype or haplotype information.

**Original Functionality:**
- Comparison of somatic variants by location and alleles using bcftools isec
- Support for specifying false-positive regions
- Option to compute the false-positive rate per megabase
- AF (allele frequency) filtering capabilities
- Stratification of variants using ambiguity regions

**References:**
- [Original som.py Documentation](sompy.md)

### 2. scmp - Specialized Comparison Engine

**Description:** scmp was a comparison engine that provided alternative comparison methods:

**Original Functionality:**
- **scmp-distancebased**: Matched variants by location only within a configurable distance
- **scmp-somatic**: Specialized mode for Tumor/Normal VCF file comparisons

**References:**
- [Original engine documentation in happy.md](happy.md#comparison-engines)

## Requesting Component Restoration

If these removed components are important for your genomic analysis workflow, please consider:

1. **Opening an Issue:** Create a GitHub issue expressing interest in having these components restored
   - Describe your use case and how you used the component
   - Specify any particular features that were most valuable

2. **Contributing:** If you have experience with Python development and bioinformatics, consider contributing to restore these components
   - The project follows modern Python development practices with type hints and pytest testing
   - Contributors can access the original implementation for reference

## Current Alternatives

While these tools are not included in this modernized version, you might consider these alternatives:

1. **For somatic variant comparison:**
   - Using the original hap.py repository's som.py tool if Python 2.7 is available
   - Using direct bcftools commands (som.py was primarily a wrapper around bcftools isec)
   - Using other somatic variant comparison tools like SMuRF or somVarIUS

2. **For distance-based or somatic comparison engines:**
   - Using the other available comparison engines in hap.py
   - Exploring external tools designed specifically for somatic variant comparison

## Future Plans

Based on community interest and contributor availability, these components may be restored in future releases. Please submit your feedback to help prioritize development efforts.

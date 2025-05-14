---
applyTo: "**/*.{py,cpp,h,hpp}"
---
# Bioinformatics Coding Standards for hap.py

## Genomic Data Handling
- Use established libraries (pysam, htslib) for handling genomic data formats
- Implement proper error checking for VCF and BAM/SAM files
- Handle large chromosomes efficiently with streaming where possible
- Validate genomic coordinates before use

## Variant Handling
- Follow VCF 4.2+ specification for variant representations
- Handle complex variants with care (MNPs, indels, structural variants)
- Use proper normalization techniques for variant comparison
- Consider phase information when relevant

## Performance Considerations
- Implement memory-efficient algorithms for large genomic regions
- Use appropriate data structures for genomic interval operations
- Consider parallelization for CPU-intensive tasks
- Use efficient C++ implementations for performance-critical sections
- Profile memory usage with genomic-scale datasets

## Bioinformatics Best Practices
- Include appropriate warnings for reference incompatibilities
- Document assumptions about input data
- Follow accepted variant calling benchmarking methodologies
- Provide clear documentation of metrics and statistical methods
- Ensure reproducibility through fixed random seeds when applicable

# Multimerge Implementation Summary

## Overview
This document summarizes the implementation of the `multimerge` functionality in the modernized hap.py codebase. The original C++ tool has been replaced with a Python implementation that maintains backwards compatibility and provides all required functionality.

## Key Features Implemented

1. **Basic Merging Functionality**
   - Merging multiple VCF files into a single output file
   - Support for sample selection via `filename:sample_name` syntax
   - Support for handling BGZ/TBI compressed VCF files

2. **Advanced Variant Processing**
   - Trimming common prefix/suffix from alleles (`--trimalleles`)
   - Left-shifting variants (`--leftshift`)
   - Handling overlapping variant locations (`--merge-by-location`)
   - Collecting unique alleles when merging variants (`--unique-alleles`)

3. **Header Handling**
   - Combining headers from all input VCF files
   - Preserving FORMAT, INFO, and FILTER entries
   - Adding command line information to the output VCF

## Implementation Details

The implementation is split into several key components:

1. **Command Line Interface**
   - Maintains compatibility with the original multimerge tool CLI
   - Added the `--process-full` option for backward compatibility

2. **Core Merging Logic**
   - `merge_vcfs()`: Main function implementing the merging algorithm
   - Handles different processing modes based on specified options

3. **Helper Functions**
   - `trim_alleles()`: Optimized allele trimming function
   - `get_sample_name()`: Logic for handling sample selection

## Usage

```bash
multimerge [options] input_files... -o output.vcf
```

### Key Options

- `-o, --output`: Output VCF file (required)
- `-r, --reference`: Reference FASTA file
- `--trimalleles`: Trim common bases from alleles (0=off, 1=on)
- `--leftshift`: Left-shift variants (0=off, 1=on)
- `--merge-by-location`: Merge variants at the same location (0=off, 1=on)
- `--unique-alleles`: Keep only unique alleles (0=off, 1=on)
- `--process-full`: Process all variants fully (0=off, 1=on)

## Testing

The multimerge implementation has been tested with:

- Basic input/output functionality tests
- Tests for allele trimming and left-shifting
- Integration tests with real VCF data

## Future Improvements

1. **Performance Optimization**
   - Efficient handling of large VCF files
   - Parallel processing options

2. **Additional Features**
   - More sophisticated sample merging strategies
   - Extended FORMAT field processing
   - Support for additional normalization options

## Related Files

- `src/hap_py/haplo/multimerge.py`: Main implementation
- `build/bin/multimerge`: Binary wrapper script

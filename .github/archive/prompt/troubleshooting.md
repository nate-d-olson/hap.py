# Troubleshooting Guide for hap.py Development

This guide provides solutions for common issues encountered when working with GitHub Copilot on the hap.py project.

## GitHub Copilot Issues

### Copilot Not Providing Relevant Suggestions

**Problem**: GitHub Copilot suggestions aren't relevant to bioinformatics or genomics.

**Solution**:
1. Provide more context in your comments
2. Include bioinformatics terms in your prompts
3. Reference specific file formats (VCF, BAM, etc.)
4. Start with existing code patterns from the codebase

**Example Prompt**:
```
@copilot I'm working on parsing VCF files where each line represents a genomic variant.
Help me write a function that filters variants based on quality score and read depth.
```

### Incorrect Python 3 Suggestions

**Problem**: Copilot suggests code that still uses Python 2 idioms.

**Solution**:
1. Explicitly mention Python 3 in your prompt
2. Use the migration examples from the `.github/prompt/migration-examples.md` file
3. Review suggestions carefully for Python 2 patterns

**Example Prompt**:
```
@copilot Using Python 3, help me update this dictionary iteration that currently uses iteritems()
```

## Build and Installation Issues

### CMake Configuration Failures

**Problem**: CMake fails to configure properly during installation.

**Solution**:
1. Check for missing dependencies
2. Ask Copilot to troubleshoot specific error messages

**Example Prompt**:
```
@copilot I'm getting this CMake error during installation. How do I fix it?
CMake Error at external/htslib.cmake:25 (message):
  Could not find zlib
```

### Python Library Import Errors

**Problem**: Import errors when running Python scripts after installation.

**Solution**:
1. Check that the build directory is in PYTHONPATH
2. Verify that dependencies are properly installed

**Example Prompt**:
```
@copilot I'm getting "ImportError: No module named 'haplotypes'" after installing. How do I fix this?
```

## Testing Issues

### Test Failures During Migration

**Problem**: Tests failing after Python 2 to 3 migration.

**Solution**:
1. Capture the exact error message
2. Identify the specific Python 3 incompatibility
3. Ask Copilot for targeted fixes

**Example Prompt**:
```
@copilot This test is failing with "TypeError: 'dict_keys' object does not support indexing" after migration to Python 3. How do I fix it?
```

### Bioinformatics-Specific Test Failures

**Problem**: Tests failing due to variant representation differences.

**Solution**:
1. Check for normalization issues
2. Verify reference handling
3. Look for coordinate system differences

**Example Prompt**:
```
@copilot The variant comparison test is failing because variants appear equivalent but are being called as different. How should I debug this?
```

## Performance Issues

### Memory Usage Problems

**Problem**: Script runs out of memory with large genomic datasets.

**Solution**:
1. Implement streaming approaches
2. Optimize data structures
3. Add chunking for large chromosomes

**Example Prompt**:
```
@copilot This variant processing code uses too much memory with large VCF files. How can I optimize it to process variants in chunks?
```

### Slow Processing

**Problem**: Script takes too long to process large datasets.

**Solution**:
1. Profile the code to identify bottlenecks
2. Consider moving critical sections to C++
3. Implement parallelization

**Example Prompt**:
```
@copilot How can I parallelize this variant normalization function to process multiple chromosomes simultaneously?
```

## Dependency Issues

### Python Package Compatibility

**Problem**: Packages have different APIs in Python 3.

**Solution**:
1. Check package documentation for migration guides
2. Update function signatures to match new APIs

**Example Prompt**:
```
@copilot How do I update this pysam function call for compatibility with the latest version?
```

### C++ Integration Issues

**Problem**: Cython interfaces breaking after Python 3 migration.

**Solution**:
1. Review Cython type definitions
2. Check for byte/string conversion issues
3. Update C++ function signatures

**Example Prompt**:
```
@copilot I'm getting a "TypeError: bytes/str mismatch" in this Cython module after migrating to Python 3. How do I fix it?
```

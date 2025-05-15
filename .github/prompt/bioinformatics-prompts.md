# Bioinformatics Prompts for GitHub Copilot

This document contains useful prompts for working with GitHub Copilot on genomics and bioinformatics tasks in the hap.py project.

## Variant Processing Prompts

### Variant Normalization
```
@copilot Help me implement a function to normalize variants by left-aligning and decomposing complex variants according to best practices in bioinformatics.
```

### VCF Parsing
```
@copilot Write a memory-efficient function to parse large VCF files that follows Python 3 best practices and handles compressed files.
```

### Variant Comparison
```
@copilot Help me implement a function to compare two variants to determine if they are equivalent, considering normalization, reference context, and representation differences.
```

## Performance Optimization Prompts

### Memory Optimization
```
@copilot This code is using too much memory with large chromosomes. How can I optimize it to use streaming or chunking approaches?
```

### Parallel Processing
```
@copilot Help me implement parallel processing for this variant analysis function to better utilize multiple CPU cores.
```

### C++ Integration
```
@copilot How can I expose this Python function to C++ for better performance in the core comparison algorithm?
```

## File Format Handling

### BAM/CRAM Processing
```
@copilot Write a function to extract read information from BAM files using pysam that works efficiently with Python 3.
```

### BED Interval Operations
```
@copilot Implement an efficient function to find overlaps between variants and regions in a BED file.
```

### FASTA Reference Handling
```
@copilot Create a function to extract reference sequences efficiently for variant normalization.
```

## Statistical Analysis

### Precision-Recall Calculation
```
@copilot Implement a function to calculate precision, recall, and F1 score for variant calls, with optional stratification.
```

### Confidence Intervals
```
@copilot Help me add bootstrap confidence interval calculations to the variant comparison metrics.
```

## Testing

### Test Data Generation
```
@copilot Create a function to generate synthetic variant data for testing, covering different variant types (SNVs, indels, MNPs).
```

### Edge Case Testing
```
@copilot Help me write tests for challenging variant comparison edge cases like complex indels, nested variants, and MNPs.
```

## Documentation

### API Documentation
```
@copilot Convert this function documentation to Google-style docstrings, including appropriate bioinformatics context.
```

### User Guide Examples
```
@copilot Help me write an example usage section for the documentation that demonstrates comparing two VCF files with stratification.
```

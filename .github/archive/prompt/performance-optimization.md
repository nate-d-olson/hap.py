# Performance Optimization Guide for hap.py

This document provides guidance on optimizing hap.py code for performance when working with large genomic datasets.

## Common Performance Bottlenecks in Genomic Analysis

### 1. File I/O

Large genomic files (VCF, BAM) can cause I/O bottlenecks.

**Ask Copilot:**
```
@copilot How can I optimize this VCF parsing function to reduce file I/O overhead?
```

**Optimization Strategies:**
- Use streaming parsers (generators/iterators)
- Implement chunked reading
- Use compressed file formats with index-based random access
- Consider memory-mapped files for large reference sequences

### 2. Variant Comparison Algorithms

Comparing variants across large datasets can be compute-intensive.

**Ask Copilot:**
```
@copilot How can I optimize this variant comparison function that compares thousands of variants?
```

**Optimization Strategies:**
- Index variants by position for O(1) lookup
- Pre-normalize variants to reduce comparison complexity
- Use efficient data structures (interval trees, hash tables)
- Consider region-based partitioning of comparisons

### 3. Memory Usage with Large Genomic Regions

Whole-genome analyses can easily exceed available memory.

**Ask Copilot:**
```
@copilot This code loads all variants into memory. How can I modify it to use streaming instead?
```

**Optimization Strategies:**
- Process by chromosome or region rather than whole genome
- Implement sliding window approaches
- Use generator functions to yield results incrementally
- Consider reference-based compression of redundant data

### 4. String Manipulation Overhead

VCF processing involves significant string operations that can be slow.

**Ask Copilot:**
```
@copilot How can I optimize these string operations in my VCF parsing code?
```

**Optimization Strategies:**
- Minimize string copies and concatenations
- Use bytes objects for performance-critical sections
- Consider Cython for string-heavy operations
- Use regular expressions efficiently (compile once, reuse)

## Optimizing Python Code

### Using Python Profiling

Profile your code to identify bottlenecks:

```python
import cProfile
cProfile.run('my_function(args)')
```

**Ask Copilot:**
```
@copilot Help me interpret this profiling output and identify bottlenecks
```

### Moving Hot Spots to C++

For critical performance sections:

**Ask Copilot:**
```
@copilot How can I move this Python function to C++ using Cython?
```

**Implementation Strategy:**
1. Identify performance-critical functions
2. Create Cython (.pyx) wrapper
3. Implement core algorithm in C++
4. Use proper type definitions in Cython interface

### Parallelization

For operations that can be parallelized:

**Ask Copilot:**
```
@copilot How can I parallelize this variant processing across multiple chromosomes?
```

**Parallelization Strategies:**
- Use `multiprocessing` for CPU-bound tasks
- Use chromosome or region-based parallelism
- Consider thread safety when accessing shared resources
- Implement proper merging of results from parallel processes

## Memory Optimization

### Reducing Memory Footprint

**Ask Copilot:**
```
@copilot How can I reduce the memory footprint of this variant data structure?
```

**Memory Optimization Strategies:**
- Use appropriate data types (e.g., int8 for small values)
- Store only essential fields, compute others on demand
- Use numeric codes instead of string labels where appropriate
- Consider specialized containers like NumPy arrays for homogeneous data

### Efficient Reference Genome Access

**Ask Copilot:**
```
@copilot What's the most memory-efficient way to access the reference genome sequence?
```

**Strategies:**
- Use indexed FASTA access (pyfaidx, pysam)
- Implement caching for commonly accessed regions
- Consider memory mapping for random access

## C++ Optimization

### Efficient Data Structures

**Ask Copilot:**
```
@copilot What's the most efficient C++ data structure for storing millions of genomic intervals?
```

**Recommendations:**
- Use `std::vector` for sequential access to variants
- Use `std::unordered_map` for position-based lookups
- Consider specialized interval trees for overlap queries
- Use move semantics to avoid unnecessary copies

### Memory Management

**Ask Copilot:**
```
@copilot How can I optimize memory management in this C++ genomic interval class?
```

**Best Practices:**
- Avoid unnecessary dynamic allocations
- Use smart pointers for managing resources
- Consider custom allocators for specialized needs
- Use string_view for non-owning string references

## I/O Optimization

### Efficient VCF/BAM Reading

**Ask Copilot:**
```
@copilot What's the fastest way to read specific regions from this VCF file?
```

**Strategies:**
- Use indexed access via tabix/CSI indexes
- Implement multi-region queries in a single pass
- Use appropriate buffer sizes for reading

## Benchmarking

Always benchmark optimizations to ensure they're effective:

**Ask Copilot:**
```
@copilot Help me write a benchmarking function to compare the performance of these two variant comparison approaches
```

**Benchmarking Approach:**
1. Create representative test data
2. Measure execution time and memory usage
3. Run multiple iterations to account for variability
4. Compare against baseline implementation

## Examples

### Example: Optimizing Variant Parsing

**Original Code:**
```python
def parse_variants(filename):
    variants = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            variants.append({
                'chrom': fields[0],
                'pos': int(fields[1]),
                'ref': fields[3],
                'alt': fields[4]
            })
    return variants
```

**Ask Copilot:**
```
@copilot Optimize this variant parsing function for memory efficiency with large VCF files
```

**Optimized Code:**
```python
def parse_variants(filename):
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            yield {
                'chrom': fields[0],
                'pos': int(fields[1]),
                'ref': fields[3],
                'alt': fields[4]
            }
```

### Example: Optimizing Variant Comparison

**Original Code:**
```python
def find_matching_variants(query_variants, truth_variants):
    matches = []
    for q_var in query_variants:
        for t_var in truth_variants:
            if (q_var['chrom'] == t_var['chrom'] and
                abs(q_var['pos'] - t_var['pos']) <= 5 and
                q_var['ref'] == t_var['ref'] and
                q_var['alt'] == t_var['alt']):
                matches.append((q_var, t_var))
    return matches
```

**Ask Copilot:**
```
@copilot Optimize this variant comparison function that's currently O(n²)
```

**Optimized Code:**
```python
def find_matching_variants(query_variants, truth_variants):
    # Index truth variants by position for fast lookup
    truth_by_pos = {}
    for t_var in truth_variants:
        chrom = t_var['chrom']
        pos = t_var['pos']
        if (chrom, pos) not in truth_by_pos:
            truth_by_pos[(chrom, pos)] = []
        truth_by_pos[(chrom, pos)].append(t_var)

    matches = []
    for q_var in query_variants:
        q_chrom = q_var['chrom']
        q_pos = q_var['pos']
        # Check within window of +/- 5bp
        for offset in range(-5, 6):
            key = (q_chrom, q_pos + offset)
            if key in truth_by_pos:
                for t_var in truth_by_pos[key]:
                    if (q_var['ref'] == t_var['ref'] and
                        q_var['alt'] == t_var['alt']):
                        matches.append((q_var, t_var))
    return matches
```

# Python 2 to 3 Migration Examples with GitHub Copilot

This document provides practical examples of using GitHub Copilot to assist with common Python 2 to 3 migration challenges in the hap.py codebase.

## String Handling

### Unicode vs Bytes

**Python 2 Code:**
```python
def process_vcf_line(line):
    fields = line.strip().split('\t')
    return fields[0] + ":" + fields[1]
```

**Copilot Prompt:**
```
@copilot Fix this Python 2 string handling for Python 3, considering that VCF files might contain non-ASCII characters
```

**Expected Python 3 Fix:**
```python
def process_vcf_line(line):
    if isinstance(line, bytes):
        line = line.decode('utf-8')
    fields = line.strip().split('\t')
    return fields[0] + ":" + fields[1]
```

## Dictionary Methods

### Dictionary Iteration

**Python 2 Code:**
```python
def count_variants_by_type(variants_dict):
    counts = {}
    for variant_id, variant in variants_dict.iteritems():
        variant_type = variant.get("TYPE", "UNKNOWN")
        counts[variant_type] = counts.get(variant_type, 0) + 1
    return counts
```

**Copilot Prompt:**
```
@copilot Update this dictionary iteration from Python 2 to Python 3
```

**Expected Python 3 Fix:**
```python
def count_variants_by_type(variants_dict):
    counts = {}
    for variant_id, variant in variants_dict.items():
        variant_type = variant.get("TYPE", "UNKNOWN")
        counts[variant_type] = counts.get(variant_type, 0) + 1
    return counts
```

## Integer Division

**Python 2 Code:**
```python
def calculate_variant_position(start, end):
    return (start + end) / 2
```

**Copilot Prompt:**
```
@copilot Fix this division operation for Python 3 to ensure integer division is preserved
```

**Expected Python 3 Fix:**
```python
def calculate_variant_position(start, end):
    return (start + end) // 2  # Use integer division
```

## Exception Handling

**Python 2 Code:**
```python
def parse_vcf_header(header_line):
    try:
        fields = header_line.split('=')
        return fields[1]
    except IndexError, e:
        logger.error("Invalid header format: %s" % header_line)
        return None
```

**Copilot Prompt:**
```
@copilot Update this exception handling syntax for Python 3
```

**Expected Python 3 Fix:**
```python
def parse_vcf_header(header_line):
    try:
        fields = header_line.split('=')
        return fields[1]
    except IndexError as e:
        logger.error("Invalid header format: %s" % header_line)
        return None
```

## Print Statements

**Python 2 Code:**
```python
def report_variant_counts(counts):
    print "Variant counts by type:"
    for variant_type, count in counts.items():
        print "%s: %d" % (variant_type, count)
```

**Copilot Prompt:**
```
@copilot Convert these print statements to Python 3
```

**Expected Python 3 Fix:**
```python
def report_variant_counts(counts):
    print("Variant counts by type:")
    for variant_type, count in counts.items():
        print("%s: %d" % (variant_type, count))
```

## File Handling

**Python 2 Code:**
```python
def read_vcf_file(filename):
    with open(filename, 'r') as f:
        for line in f:
            yield line
```

**Copilot Prompt:**
```
@copilot Update this file reading code for Python 3, considering potential encoding issues with VCF files
```

**Expected Python 3 Fix:**
```python
def read_vcf_file(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            yield line
```

## Module Imports

**Python 2 Code:**
```python
import ConfigParser
import StringIO
```

**Copilot Prompt:**
```
@copilot Update these Python 2 imports to their Python 3 equivalents
```

**Expected Python 3 Fix:**
```python
import configparser
import io
```

## Using These Examples

1. Find similar patterns in the hap.py codebase
2. Copy the relevant Copilot prompt
3. Adapt it to your specific code context
4. Use GitHub Copilot Chat to get specific solutions

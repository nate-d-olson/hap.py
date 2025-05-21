# String Handling in Python 3 Migration

## Overview

One of the biggest challenges in migrating hap.py from Python 2 to Python 3 is handling the str/bytes distinction. This document outlines the approach and utilities developed for handling string conversion consistently, particularly at the C++/Python boundary.

## Background

In Python 2, strings were represented as `str` (byte strings) and `unicode` (Unicode strings). In Python 3, this model changed:

- `str` is now a Unicode string
- `bytes` is a binary data type
- No automatic conversion between `str` and `bytes` in Python 3

This change significantly impacts code that:

1. Interfaces with C++ libraries expecting byte strings
2. Handles file I/O in text vs binary mode
3. Processes biological data that might contain non-ASCII characters

## String Handling Strategy

The migration implements a consistent string handling strategy using dedicated utility functions:

### Core Utility Functions

The `Haplo.string_handling` module provides these key functions:

```python
def ensure_str(text: Union[str, bytes, None]) -> Optional[str]:
    """Convert from bytes to str if needed, or return as str"""

def ensure_bytes(text: Union[str, bytes, None]) -> Optional[bytes]:
    """Convert from str to bytes if needed, or return as bytes"""

def ensure_text_io(text: Union[str, bytes], file_mode: str) -> Union[str, bytes]:
    """Convert text to the appropriate type for the file mode"""
```

### Usage Guidelines

#### When to use each function

1. **ensure_str()**: Use when receiving data from C++ that might be bytes but you need a string for Python processing

   ```python
   # Example with pysam/C++ function returning potential bytes
   record_id = ensure_str(record.id)
   ```

2. **ensure_bytes()**: Use when passing data to C++ that expects bytes

   ```python
   # Example with C++ function expecting bytes
   c_function(ensure_bytes(chromosome_name))
   ```

3. **ensure_text_io()**: Use when reading/writing files with different modes

   ```python
   # Example with file operations
   mode = "rb" if binary_format else "r"
   content = ensure_text_io(data, mode)
   ```

### Common Patterns

#### File I/O

```python
# Reading files
with open(filename, mode) as f:
    content = f.read()
    processed_content = ensure_str(content)

# Writing files
with open(filename, mode) as f:
    f.write(ensure_text_io(content, mode))
```

#### VCF Processing

```python
# Processing VCF records
for record in vcf_reader:
    chrom = ensure_str(record.chrom)
    ref = ensure_str(record.ref)
```

## Testing String Handling

To ensure correct string handling, test both with:

1. ASCII input (standard DNA sequences)
2. Non-ASCII input (comments or IDs with special characters)
3. Binary data (when appropriate)

All string handling functions have unit tests in `test_string_handling.py`.

## Common Issues and Solutions

### Issue: Format String Type Mismatches

Problem:

```python
# Error: cannot interpolate bytes into string
message = f"Processing chromosome {bytes_chromosome}"
```

Solution:

```python
message = f"Processing chromosome {ensure_str(bytes_chromosome)}"
```

### Issue: Bytes Arguments to String Methods

Problem:

```python
# Error: startswith() requires str, not bytes
if bytes_sequence.startswith("ACGT"):
    # ...
```

Solution:

```python
if ensure_str(bytes_sequence).startswith("ACGT"):
    # ...
```

### Issue: Binary File Mode with String Data

Problem:

```python
# Error: write() argument must be bytes
with open(filename, "wb") as f:
    f.write(string_data)
```

Solution:

```python
with open(filename, "wb") as f:
    f.write(ensure_bytes(string_data))
```

## Further Reading

- [Python 3 Text vs. Data](https://docs.python.org/3/howto/unicode.html)
- [Porting Python 2 Code to Python 3](https://docs.python.org/3/howto/pyporting.html)
- [Python Bytes Objects](https://docs.python.org/3/library/stdtypes.html#bytes)

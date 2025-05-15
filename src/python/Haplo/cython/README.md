# Haplo Cython Module

This directory contains the Cython modules for the hap.py Python/C++ integration.

## Overview

The hap.py tool uses C++ for performance-critical algorithms and Python for higher-level functionality.
The integration between these components is handled through Cython modules in this directory.

## Python 3 Migration Notes

In the Python 3 migration, we've updated the Cython integration to:
1. Use Python 3 string handling (Unicode vs bytes)
2. Support modern Cython language features
3. Provide proper error handling and fallbacks
4. Use Python 3's improved memory management

## Files

- `__init__.py` - Package initialization, handles import failures gracefully
- `_internal.pyx` - Main Cython module interfacing with C++ code
- `mock_internal.py` - Mock implementation for testing without C++ dependencies
- `CMakeLists.txt` - Build configuration for Cython module

## Usage

The Cython module is used automatically when importing from `Haplo`. For example:

```python
from Haplo import quantify

# This will use the C++ implementation automatically if available
result = quantify.analyze_variants(truth_vcf, query_vcf)
```

## Development

When adding new functionality:

1. Define the C++ interface in `_internal.pyx`
2. Update the mock implementation in `mock_internal.py`
3. Test both implementations

## Testing Without C++ Components

During development or for environments where C++ compilation is not available:

```python
import os
os.environ["HAPLO_USE_MOCK"] = "1"
import Haplo
```

This will force the use of the Python mock implementation.

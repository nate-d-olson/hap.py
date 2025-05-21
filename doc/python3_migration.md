# Python 3 Migration Guide for hap.py

This document provides information on the Python 3 migration of the hap.py toolkit, including changes to the API, compatibility notes, and instructions for users migrating from the Python 2 version.

## Background

The original hap.py toolkit was developed with Python 2, which reached end-of-life on January 1, 2020. This migration updates the codebase to use Python 3, ensuring continued support and development.

## Major Changes

### Python 3 Compatibility

- All code has been updated to be compatible with Python 3.7 and later
- String handling has been updated to handle Unicode correctly
- File I/O operations now use proper encoding/decoding
- Exception handling uses the modern Python 3 patterns
- Codebase streamlined to focus on core functionality only

### Cython Integration

- Cython modules have been updated for Python 3 compatibility
- String handling at the C++/Python boundary has been improved
- Mock implementations for testing without C++ components
- Better error handling in Cython wrapper modules

### CLI Tools

- Command-line interfaces have been modernized and streamlined
- Only core CLI tools (hap.py, pre.py, qfy.py) are maintained
- Auxiliary tools (ovc.py, cnx.py) have been removed
- Entry points properly defined in pyproject.toml
- All entry point functions return appropriate exit codes
- Improved error reporting and handling
- Better logging and progress reporting
- All scripts now use proper Python 3 shebangs
- Type hints added to main functions
- Consistent error handling pattern across all tools
- See [CLI Updates](cli_updates.md) for detailed information

### Dependencies

- Updated requirements for Python 3 compatibility:
  - pysam ≥0.16.0
  - numpy ≥1.19.0
  - pandas ≥1.3.0
  - scipy ≥1.7.0
  - psutil ≥5.8.0

## Installation

### Using pip (Recommended)

```bash
# Install from PyPI
pip install hap.py

# Or install from source directory
git clone https://github.com/Illumina/hap.py.git
cd hap.py
pip install .
```

### Using the install script

```bash
# Clone the repository
git clone https://github.com/Illumina/hap.py.git
cd hap.py

# Run the install script with a build directory
python3 install.py /path/to/build/dir
```

## String Handling in Python 3

One of the biggest changes in Python 3 is how strings are handled. Python 3 makes a clear distinction between text strings (`str`, which are Unicode) and binary data (`bytes`). This affects many parts of the code, especially at the interface between Python and C++.

### String Handling Utilities

We've created a dedicated module `Haplo.string_handling` to handle string conversion consistently:

```python
from Haplo.string_handling import ensure_str, ensure_bytes, ensure_text_io

# Convert potential bytes to str
text = ensure_str(bytes_or_str)  # Always returns str

# Convert potential str to bytes
binary = ensure_bytes(str_or_bytes)  # Always returns bytes

# Handle file I/O text conversion based on mode
content = ensure_text_io(content, file_mode)  # Returns appropriate type for the mode
```

### Common String Issues

When migrating Python 2 code to Python 3, watch for these common issues:

1. **String literals in Cython**:
   Python 3 string literals are Unicode by default, while C++ functions often expect bytes.

   ```python
   # Wrong
   c_function("string")  # May fail in Python 3

   # Correct
   c_function(b"string")  # Use bytes literal
   # or
   c_function(ensure_bytes("string"))  # Use conversion function
   ```

2. **File I/O**:

   ```python
   # Python 2 style (problematic in Python 3)
   f = open(filename, "r")
   content = f.read()

   # Python 3 style (preferred)
   with open(filename, "r", encoding="utf-8") as f:
       content = f.read()
   ```

3. **Bytestrings in output**:

   ```python
   # May get b'text' in output
   data = some_c_function()
   print(data)  # Might print b'text' instead of 'text'

   # Fix with ensure_str
   print(ensure_str(data))  # Will print 'text'
   ```

## Testing with Mock Implementations

For testing without requiring C++ components, we've implemented mock versions of the Cython modules:

```python
# Set environment variable to use mocks
os.environ["HAPLO_USE_MOCK"] = "1"

# Import modules as normal - they'll use the mock implementations
from Haplo.cython import _internal
```

The mock implementations provide simplified versions of:

- `complement_sequence` and `reverse_complement`
- `read_fasta_index` and `get_reference_sequence`
- `parse_vcf_file` and `write_vcf_file`
- `MockVariantRecord` and `MockHaploCompare` classes

## Integration Tests

The following integration tests verify Python 3 compatibility:

- `tests/test_cython_integration.py`: Verifies Cython module loading
- `tests/test_happy_integration.py`: Tests the main hap.py CLI
- `tests/test_qfy_integration.py`: Tests the quantification module

Run tests with:

```bash
python -m pytest tests/
```

## Migrating Custom Scripts

If you have custom scripts built on top of hap.py, you'll need to update them for Python 3 compatibility. Key areas to check:

1. **String handling**: Update string operations for Python 3's unicode/bytes distinction
2. **Print statements**: Ensure print statements use parentheses `print("text")`
3. **Division**: Integer division now needs `//` instead of `/`
4. **Imports**: Update any imports of moved standard library modules
5. **Exception handling**: Use `except Exception as e` syntax
6. **File I/O**: Add encoding parameters and use context managers

## API Changes

The Python 3 migration maintains API compatibility where possible, but some changes were unavoidable:

- Functions that returned strings in Python 2 now return Unicode strings in Python 3
- Functions expecting binary data now require bytes objects
- File path handling now uses `pathlib.Path` objects internally

## Troubleshooting

### TypeError: a bytes-like object is required, not 'str'

This common error occurs when passing a string to a function that expects bytes. Solutions:

```python
# Option 1: Use bytes literals
function(b"string")

# Option 2: Use the conversion utility
from Haplo.string_handling import ensure_bytes
function(ensure_bytes("string"))
```

### UnicodeDecodeError: 'utf-8' codec can't decode byte

This error occurs when trying to decode binary data that isn't valid UTF-8. Solutions:

```python
# Option 1: Specify error handling
text = binary.decode('utf-8', errors='replace')

# Option 2: Use the conversion utility
from Haplo.string_handling import ensure_str
text = ensure_str(binary)  # Uses 'replace' error handling by default
```

### ModuleNotFoundError: No module named 'Haplo.cython._internal'

This error occurs when the Cython modules haven't been built. Solutions:

1. Run the install script: `python3 install.py /path/to/build/dir`
2. Use mock implementations for testing:

   ```python
   import os
   os.environ["HAPLO_USE_MOCK"] = "1"
   ```

pip install .

```

This will build all necessary components including the C++ parts and install the Python package with command-line entry points.

To install with optional dependencies for C++/Cython extensions (recommended for performance) or development tools:

```bash
pip install .[cpp]      # For C++/Cython accelerated features
pip install .[dev]      # For development tools (testing, linting)
pip install .[cpp,dev]  # For both
```

### Building from Source

If you need to build from source and `pip install .` does not meet your needs:

1. **Prerequisites**:
   - A C++14 compatible compiler (e.g., GCC, Clang, MSVC)
   - CMake (version 3.10 or newer)
   - Python (version 3.7 or newer, including development headers)
   - Boost libraries (version 1.55.0 or newer)
   - Zlib development libraries

2. **Using the install script**:

   ```bash
   python install.py /path/to/install/dir
   ```

## API Changes

### String Handling

In Python 3, all strings are Unicode by default. When interfacing with external components like VCF files or the C++ code, proper encoding/decoding is needed:

```python
# Python 2 (old)
some_text = str(vcf_record.REF)  # This could be a byte string

# Python 3 (new)
some_text = str(vcf_record.REF)  # This is a Unicode string
```

### Exception Handling

Exception handling has been updated to use modern Python 3 patterns:

```python
# Python 2 (old)
try:
    # some code
except IOError:
    # handle error

# Python 3 (new)
try:
    # some code
except IOError as e:
    # handle error with context
    raise ValueError("Failed to process file") from e
```

### File I/O

File operations now require text mode and encoding to be specified explicitly:

```python
# Python 2 (old)
with open(filename, 'w') as f:
    f.write(data)

# Python 3 (new)
with open(filename, 'w', encoding='utf-8') as f:
    f.write(data)
```

## Compatibility Notes

### Backward Compatibility

Most core command-line interfaces should work the same as in the Python 2 version, but with improved error handling and reporting. However, several non-core components have been removed:

- **Removed CLI tools**: ovc.py, cnx.py
- **Removed modules**: bamstats.py
- **Removed C++ components**: scmp directory, XCmpQuantify
- **Removed functionality**: Somatic variant calling, xcmp comparison engine

The Python 3 version focuses exclusively on the vcfeval comparison engine, which is now the only supported engine option.

### VCF Files

The toolkit is compatible with both VCF and gVCF files. When working with older VCFs that might have encoding issues, use the `--fixchr` option to normalize chromosome names.

### Testing Without C++ Components

For testing or environments where building the C++ components is difficult, a mock implementation is available:

```python
# Set environment variable to use mock implementation
os.environ["HAPLO_USE_MOCK"] = "1"
import Haplo
```

## Known Issues

- Unicode in variant IDs or sample names might cause issues with some VCF parsers
- Running with the mock implementation is significantly slower than with C++ components
- Some advanced features require C++ components and have no mock implementation

## Getting Help

If you encounter issues with the Python 3 version of hap.py, please open an issue on GitHub with details about your environment, commands run, and any error messages.

# Python 3 Migration Tools

This document describes the tools created to assist with the Python 3 migration of hap.py.

## Migration Status

As of May 17, 2025, the migration is **93.9% complete** with 46/49 files updated and only 14 issues remaining across 3 files (down from 28 issues in 7 files at the previous checkpoint).

## Available Tools

### 1. `fix_remaining_py3_issues.py`

This script automatically fixes common Python 3 compatibility issues in the remaining problematic files.

**Features:**
- Fixes exception syntax (`except Exception, e` → `except Exception as e`)
- Adds encoding parameters to file operations
- Updates subprocess calls to use text mode
- Improves string/unicode handling

**Usage:**

```bash
# Check for issues without applying fixes
./fix_remaining_py3_issues.py

# Check files and apply fixes
./fix_remaining_py3_issues.py --apply

# Process specific files
./fix_remaining_py3_issues.py --files file1.py,file2.py --apply

# Show verbose output
./fix_remaining_py3_issues.py --verbose
```

### 2. `verify_py3_build.py`

This script verifies the Python 3 build system compatibility, checking CMake files, installation scripts, and required dependencies.

**Features:**
- Checks CMake configuration files for Python 3 compatibility
- Verifies installation scripts for Python 3 syntax
- Validates that required Python dependencies are installed
- Tests CMake's ability to detect Python 3

**Usage:**

```bash
# Run verification checks
./verify_py3_build.py

# Run checks and fix common issues
./verify_py3_build.py --fix

# Skip the CMake test
./verify_py3_build.py --skip-test

# Show verbose output
./verify_py3_build.py --verbose
```

### 3. Other Helper Scripts

Additionally, several other tools are available to help with the migration:

- `fix_exception_syntax.py`: Specifically targets the Python 2 style exception syntax
- `fix_string_unicode.py`: Fixes string and unicode handling issues in Python files
- `update_shebangs.py`: Updates the shebang lines in Python files to use Python 3
- `check_py3_issues.py`: Identifies Python 3 compatibility issues in the codebase
- `test_py3_core.sh`: Tests the core functionality with Python 3
- `test_cython_module_py3.py`: Tests Cython modules with Python 3

## Recommended Workflow

1. **Identify remaining issues:**
   ```bash
   python3 check_py3_issues.py --verbose
   ```

2. **Fix remaining Python 3 compatibility issues:**
   ```bash
   ./fix_remaining_py3_issues.py --apply
   ```

3. **Verify build system compatibility:**
   ```bash
   ./verify_py3_build.py --fix
   ```

4. **Run the updated Python 3 core tests:**
   ```bash
   ./test_py3_core.sh
   ```

5. **Test Cython integration:**
   ```bash
   python3 test_cython_module_py3.py
   ```

6. **Run comprehensive integration tests:**
   ```bash
   ./test_cython_integration_py3.sh
   ```

## Common Issues and Solutions

### String and Unicode Issues

In Python 3, there's a clear distinction between binary data (bytes) and text strings (str). The following changes are often needed:

1. **File operations:** Add `encoding="utf-8"` parameter:
   ```python
   # Python 2
   with open(filename, 'r') as f:

   # Python 3
   with open(filename, 'r', encoding="utf-8") as f:
   ```

2. **Subprocess output handling:** Use `universal_newlines=True` or `text=True`:
   ```python
   # Python 2
   output = subprocess.check_output(["command"])

   # Python 3
   output = subprocess.check_output(["command"], universal_newlines=True)
   ```

3. **Bytes/string conversion:** Explicitly decode bytes to strings:
   ```python
   # Python 2
   text = str(some_bytes)

   # Python 3
   text = some_bytes.decode('utf-8') if isinstance(some_bytes, bytes) else str(some_bytes)
   ```

### Exception Handling

In Python 3, the syntax for catching exceptions has changed:

```python
# Python 2
try:
    do_something()
except Exception, e:
    handle_error(e)

# Python 3
try:
    do_something()
except Exception as e:
    handle_error(e)
```

### Division Operator

In Python 3, the `/` operator performs true division (returning a float), while the `//` operator performs floor division:

```python
# Python 2
result = 5 / 2  # Result: 2

# Python 3
result = 5 / 2  # Result: 2.5
result = 5 // 2  # Result: 2
```

## More Information

For more details on the Python 3 migration, refer to:
- `PYTHON3_MIGRATION_FINAL.md`: Final status report and test plan
- `PYTHON3_MIGRATION_STATUS.md`: Current migration status
- `PYTHON3_CORE.md`: Documentation on core functionality

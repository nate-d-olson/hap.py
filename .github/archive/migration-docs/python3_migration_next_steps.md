# Python 3 Migration Next Steps

This document outlines the specific next actions to continue the Python 3 migration of the hap.py project.

## 1. Complete Code Conversion

Run the enhanced automated migration script on all remaining Python modules:

```bash
# Run the automated updates on all files
python3 automate_py3_updates.py --paths "src/python" --apply

# Verify fixes with the checking script
python3 check_py3_issues.py
```

## 2. Fix String Encoding Issues

Many files require specific handling for string encoding. Focus on:

- File I/O operations (add encoding='utf-8')
- Binary data handling in VCF processing
- Proper bytes vs str handling in Cython interfaces

## 3. Test Cython Module Integration

Test the mock implementations to ensure they work correctly:

```bash
# Test mock implementations
HAPLO_USE_MOCK=1 python3 test_cython_mock.py

# Test basic functionality
python3 test_cython_module_py3.py
```

## 4. Build System Integration

Test the full build process with Python 3:

```bash
# Test full build with Python 3
python3 install_py3.py /tmp/happy-py3-build

# Check for errors
cd /tmp/happy-py3-build
src/sh/run_tests.sh
```

## 5. Testing Core Functionality

Run comprehensive tests to compare Python 2 and 3 versions:

```bash
# Run test script
./test_happy_py3.sh

# Compare outputs
python3 compare_outputs.py
```

## 6. Documentation Updates

Update documentation to reflect Python 3 compatibility:

- Update main README.md
- Update PYTHON3.md with migration notes
- Update example scripts

## Resource Constraints

When working with large genomic data:

- Monitor memory usage of Python 3 code
- Test with real-world VCF files
- Check performance impact of string handling changes

## Tracking Progress

Update the progress tracking file after completing each module:

```bash
# After each significant change
python3 update_migration_progress.py --module "ModuleName" --status "completed"
```

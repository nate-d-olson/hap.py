# Python 3 Migration Scripts Guide

This document provides an overview of the scripts available for the Python 3 migration of hap.py.

## Primary Migration Scripts

### Code Migration

- **`automate_py3_updates.py`**: Comprehensive tool for automating Python 2 to 3 updates
  ```bash
  python3 automate_py3_updates.py --paths "src/python" --apply
  ```

- **`check_py3_issues.py`**: Tool for identifying Python 2/3 compatibility issues and tracking migration progress
  ```bash
  python3 check_py3_issues.py --paths "src/python" --generate-report
  ```

- **`update_cython_modules_py3.py`**: Tool for updating Cython modules for Python 3 compatibility
  ```bash
  python3 update_cython_modules_py3.py --src-dir "src/python"
  ```

- **`py2to3_migrate.sh`**: Main shell script for orchestrating the migration process
  ```bash
  bash py2to3_migrate.sh
  ```

### Testing

- **`test_python3_compatibility_enhanced.py`**: Comprehensive tool for testing Python 3 module compatibility
  ```bash
  python3 test_python3_compatibility_enhanced.py --build-dir /path/to/build
  ```

- **`test_cython_module_py3.py`**: Tool for testing Cython modules with Python 3
  ```bash
  # Test with real implementation
  python3 test_cython_module_py3.py --build-dir /path/to/build

  # Test with mock implementation
  python3 test_cython_module_py3.py --build-dir /path/to/build --mock
  ```

- **`validate_cpp_integration_enhanced.py`**: Tool for validating C++ integration with Python 3
  ```bash
  # Test with real implementation
  python3 validate_cpp_integration_enhanced.py --build-dir /path/to/build

  # Test with mock implementation
  python3 validate_cpp_integration_enhanced.py --build-dir /path/to/build --use-mock
  ```

- **`test_happy_py3.sh`**: Script for testing hap.py functionality with Python 3
  ```bash
  bash test_happy_py3.sh /path/to/build
  ```

- **`test_cpp_integration.sh`**: Shell script for testing C++ integration with Python 3
  ```bash
  bash test_cpp_integration.sh /path/to/build
  ```

## Migration Workflow

1. **Prepare Environment**
   ```bash
   # Create Python 3 virtual environment
   python3 -m venv /path/to/venv
   source /path/to/venv/bin/activate

   # Install dependencies
   pip install -r happy.requirements.py3.txt
   ```

2. **Automate Python 2 to 3 Conversion**
   ```bash
   # Run the automated migration tool
   python3 automate_py3_updates.py --paths "src/python" --apply

   # Check for remaining issues
   python3 check_py3_issues.py --paths "src/python" --generate-report
   ```

3. **Update Cython Modules**
   ```bash
   # Update Cython modules
   python3 update_cython_modules_py3.py --src-dir "src/python"
   ```

4. **Build and Test**
   ```bash
   # Build with Python 3
   python3 install_py3.py /tmp/happy-py3-build

   # Test Python modules
   python3 test_python3_compatibility_enhanced.py --build-dir /tmp/happy-py3-build

   # Test Cython modules
   python3 test_cython_module_py3.py --build-dir /tmp/happy-py3-build

   # Test C++ integration
   python3 validate_cpp_integration_enhanced.py --build-dir /tmp/happy-py3-build

   # Test main functionality
   bash test_happy_py3.sh /tmp/happy-py3-build
   ```

## Deprecated Scripts

The following scripts have been deprecated and should not be used:

- `convert_py2to3.py`: Use `automate_py3_updates.py` instead
- `test_python3_compatibility.py`: Use `test_python3_compatibility_enhanced.py` instead
- `test_migrated_modules.py`: Use `test_python3_compatibility_enhanced.py` instead
- `validate_cpp_integration.py`: Use `validate_cpp_integration_enhanced.py` instead
- `test_cython_integration.py`: Use `test_cython_module_py3.py` instead
- `test_cython_mock.py`: Use `test_cython_module_py3.py --mock` instead

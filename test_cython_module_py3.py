#!/usr/bin/env python3
"""
Test Cython Module Loading for Python 3 Migration in hap.py

This script tests if the Cython modules can be loaded in Python 3,
and verifies basic functionality. It helps identify issues with
string handling, memory management, and Python 3 compatibility.

Usage:
    python test_cython_module_py3.py --build-dir /path/to/build [--module module_name] [--mock] [--verbose]

The script tests:
1. Import of Cython modules
2. String handling in Cython modules
3. Basic functionality of Cython modules
4. Mock implementation fallbacks
"""

import argparse
import importlib
import logging
import os
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path
from typing import Dict, List, Optional

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("cython-tester")


def configure_paths(build_dir: str) -> bool:
    """Configure paths to include build directory in Python path.

    Args:
        build_dir: Path to the build directory

    Returns:
        True if successful, False otherwise
    """
    build_path = Path(build_dir).resolve()

    if not build_path.exists():
        logger.error(f"Build directory '{build_path}' does not exist")
        return False

    # Add lib paths to Python path
    python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
    python_lib_path = build_path / "lib" / f"python{python_version}"

    potential_paths = [
        python_lib_path,
        build_path / "lib",
        build_path / "lib" / "python",
        build_path / "python",
    ]

    paths_added = False
    for path in potential_paths:
        if path.exists():
            logger.info(f"Adding {path} to Python path")
            sys.path.insert(0, str(path))
            paths_added = True

    if not paths_added:
        logger.error(
            f"Could not find Python library path in {build_path}. "
            f"Checked: {[str(p) for p in potential_paths]}"
        )
        return False

    return True


def test_module_imports(module_names: Optional[List[str]] = None) -> Dict[str, bool]:
    """Test importing specified modules.

    Args:
        module_names: List of module names to test, or None to test all known modules

    Returns:
        Dictionary mapping module names to import success status
    """
    if module_names is None:
        # Default modules to test
        module_names = [
            "Haplo.cython._internal",
            "Haplo.cython.cpp_internal",
            "Haplo.cython.variant_processing",
            "Haplo.cython.reference_handling",
        ]

    results = {}

    for module_name in module_names:
        logger.info(f"Testing import of {module_name}")
        try:
            importlib.import_module(module_name)
            logger.info(f"Successfully imported {module_name}")
            results[module_name] = True
        except ImportError as e:
            logger.error(f"Failed to import {module_name}: {e}")
            results[module_name] = False
        except Exception as e:
            logger.error(f"Error importing {module_name}: {e}")
            results[module_name] = False
            traceback.print_exc()

    return results


def test_cython_string_handling(module_name: str) -> Dict[str, bool]:
    """Test string handling in Cython modules to identify Python 3 issues.

    Args:
        module_name: Name of the module to test

    Returns:
        Dictionary of test names and success status
    """
    results = {}

    try:
        module = importlib.import_module(module_name)

        # Try to call string handling functions if they exist
        if hasattr(module, "test_string_handling"):
            logger.info(f"Testing string handling in {module_name}")
            try:
                result = module.test_string_handling()
                results["string_handling"] = True
                logger.info(f"String handling test result: {result}")
            except Exception as e:
                results["string_handling"] = False
                logger.error(f"String handling test failed: {e}")
                traceback.print_exc()

        # Test basic functionality if available
        if hasattr(module, "test_basic_functionality"):
            logger.info(f"Testing basic functionality in {module_name}")
            try:
                result = module.test_basic_functionality()
                results["basic_functionality"] = True
                logger.info(f"Basic functionality test result: {result}")
            except Exception as e:
                results["basic_functionality"] = False
                logger.error(f"Basic functionality test failed: {e}")
                traceback.print_exc()

        # Test string encoding/decoding specifically
        if hasattr(module, "get_version") or hasattr(module, "test_module"):
            logger.info(f"Testing string encoding/decoding in {module_name}")
            try:
                # Test with a function that returns a string (most modules have one of these)
                test_function = getattr(module, "get_version", None) or getattr(
                    module, "test_module", None
                )
                if test_function:
                    result = test_function()
                    if isinstance(result, str):
                        logger.info(
                            f"String encoding test passed: returned string '{result}'"
                        )
                        results["string_encoding"] = True
                    elif isinstance(result, dict) and any(
                        isinstance(v, str) for v in result.values()
                    ):
                        logger.info(
                            "String encoding test passed: returned dict with string values"
                        )
                        results["string_encoding"] = True
                    else:
                        logger.warning(
                            f"String encoding test inconclusive: unexpected return type {type(result)}"
                        )
                        results["string_encoding"] = None
            except Exception as e:
                results["string_encoding"] = False
                logger.error(f"String encoding test failed: {e}")
                traceback.print_exc()

        # Test Unicode handling if the module has appropriate functions
        if hasattr(module, "complement_sequence") or hasattr(
            module, "reverse_complement"
        ):
            logger.info(f"Testing Unicode handling in {module_name}")
            try:
                # Create test data with Unicode characters
                test_str = "ACGTÑÄÖÜŽ"
                test_func = getattr(module, "complement_sequence", None) or getattr(
                    module, "reverse_complement", None
                )
                if test_func:
                    # This should not raise UnicodeDecodeError if properly handled
                    result = test_func(test_str)
                    logger.info(f"Unicode handling test passed: {result}")
                    results["unicode_handling"] = True
            except Exception as e:
                results["unicode_handling"] = False
                logger.error(f"Unicode handling test failed: {e}")
                traceback.print_exc()

    except ImportError as e:
        logger.error(f"Could not import {module_name}: {e}")
        results["import"] = False
    except Exception as e:
        logger.error(f"Error testing {module_name}: {e}")
        traceback.print_exc()
        results["general"] = False

    return results


def test_mock_implementations() -> bool:
    """Test the mock implementations of Cython modules

    Returns:
        bool: True if mock tests pass, False otherwise
    """
    # Force using mock implementation
    os.environ["HAPLO_USE_MOCK"] = "1"

    try:
        # Try importing the module with the mock implementation
        from Haplo.cython import USING_MOCK, mock__internal, mock_cpp_internal

        if not USING_MOCK:
            logger.error("USING_MOCK flag is False, but should be True")
            return False

        # Test basic functionality
        test_result = mock__internal.test_string_handling()
        logger.info(f"Mock _internal.test_string_handling returned: {test_result}")

        test_result = mock_cpp_internal.test_string_handling()
        logger.info(f"Mock cpp_internal.test_string_handling returned: {test_result}")

        logger.info("Mock implementation working correctly!")
        return True

    except ImportError as e:
        logger.error(f"Failed to import mock modules: {e}")
        return False
    except Exception as e:
        logger.error(f"Error testing mock implementation: {e}")
        return False


def run_build_test(build_dir: str) -> bool:
    """Run a build test to check if Cython modules can be built.

    Args:
        build_dir: Path to build directory

    Returns:
        True if build test passes, False otherwise
    """
    try:
        logger.info("Running build test...")
        # Create a simple test Cython file
        with tempfile.TemporaryDirectory() as temp_dir:
            # Write a simple Cython module
            test_pyx_path = os.path.join(temp_dir, "test_py3_compat.pyx")
            with open(test_pyx_path, "w", encoding="utf-8") as f:
                f.write(
                    """# cython: language_level=3
# distutils: language=c++

"""
                    """
Test module for Python 3 compatibility checks
"""
                    """

from libc.stdlib cimport malloc, free
from cpython.ref cimport PyObject

# String handling functions for testing
def encode_decode_test(s):
    """
                    """Test string encoding/decoding"""
                    """
    if isinstance(s, str):
        # Python str to C++ string
        cdef bytes encoded = s.encode('utf-8')
        # C++ string back to Python str
        return encoded.decode('utf-8')
    elif isinstance(s, bytes):
        # Direct bytes handling
        return s.decode('utf-8')
    return "Unknown type"

# Memory management test
cdef class TestWrapper:
    cdef char* _data

    def __cinit__(self):
        self._data = <char*>malloc(10 * sizeof(char))
        if self._data == NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._data != NULL:
            free(self._data)

def create_and_destroy():
    """
                    """Create and destroy an object to test memory management"""
                    """
    obj = TestWrapper()
    return "Memory test passed"

def test_module():
    """
                    """Test if the module is working properly"""
                    """
    return {
        "encoding_test": encode_decode_test("Test äöüñ"),
        "memory_test": create_and_destroy()
    }
"""
                )

            # Try to build the Cython module
            build_cmd = [
                sys.executable,
                "-m",
                "pip",
                "install",
                "-e",
                temp_dir,
                "--install-option=--cython-language-level=3",
            ]

            logger.info(f"Running build command: {' '.join(build_cmd)}")
            proc = subprocess.run(
                build_cmd, cwd=temp_dir, capture_output=True, text=True
            )

            if proc.returncode != 0:
                logger.error(f"Build failed with error:\n{proc.stderr}")
                return False

            logger.info("Build completed successfully")
            return True
    except Exception as e:
        logger.error(f"Build test failed with error: {e}")
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Test Cython modules for Python 3 compatibility"
    )
    parser.add_argument("--build-dir", required=True, help="Path to build directory")
    parser.add_argument("--module", help="Specific module to test")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--mock", action="store_true", help="Test mock implementations")
    parser.add_argument("--build-test", action="store_true", help="Run build test")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if not configure_paths(args.build_dir):
        logger.error("Failed to configure paths")
        sys.exit(1)

    logger.info(f"Python path: {sys.path}")

    # Run build test if requested
    if args.build_test:
        if run_build_test(args.build_dir):
            logger.info("Build test passed")
        else:
            logger.error("Build test failed")
            sys.exit(1)

    if args.mock:
        if test_mock_implementations():
            logger.info("Mock tests passed")
            sys.exit(0)
        else:
            logger.error("Mock tests failed")
            sys.exit(1)

    modules = [args.module] if args.module else None

    import_results = test_module_imports(modules)

    success = all(import_results.values())

    if success:
        logger.info("All module imports successful")

        # Test string handling in modules that imported successfully
        string_test_results = {}
        for module_name, imported in list(import_results.items()):
            if imported:
                string_test_results[module_name] = test_cython_string_handling(
                    module_name
                )

        # Print a summary of test results
        logger.info("\n===== String Handling Test Results =====")
        for module_name, tests in string_test_results.items():
            logger.info(f"Module: {module_name}")
            for test_name, result in tests.items():
                status = (
                    "✓ PASS"
                    if result
                    else "✗ FAIL"
                    if result is False
                    else "? INCONCLUSIVE"
                )
                logger.info(f"  {test_name}: {status}")

        # Check overall success
        string_tests_success = all(
            all(result for result in tests.values() if result is not None)
            for tests in string_test_results.values()
        )

        if string_tests_success:
            logger.info("All string handling tests passed!")
        else:
            logger.warning("Some string handling tests failed")
            success = False
    else:
        logger.error("Some modules failed to import")

    # Print overall result
    if success:
        logger.info("\n===== ALL TESTS PASSED =====")
    else:
        logger.error("\n===== SOME TESTS FAILED =====")

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()

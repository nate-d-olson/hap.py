#!/usr/bin/env python3
"""
Test Cython Module Loading for Python 3 Migration in hap.py

This script tests if the Cython modules can be loaded in Python 3,
and verifies basic functionality. It helps identify issues with
string handling, memory management, and Python 3 compatibility.

Usage:
    python test_cython_module_py3.py --build-dir /path/to/build [--module module_name] [--mock]
"""

import argparse
import importlib
import logging
import os
import sys
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

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
            module = importlib.import_module(module_name)
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
        import Haplo.cython
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


def main():
    parser = argparse.ArgumentParser(
        description="Test Cython modules for Python 3 compatibility"
    )
    parser.add_argument("--build-dir", required=True, help="Path to build directory")
    parser.add_argument("--module", help="Specific module to test")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--mock", action="store_true", help="Test mock implementations")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if not configure_paths(args.build_dir):
        logger.error("Failed to configure paths")
        sys.exit(1)

    logger.info(f"Python path: {sys.path}")

    if args.mock:
        if test_mock_implementations():
            logger.info("Mock tests passed")
            sys.exit(0)
        else:
            logger.error("Mock tests failed")
            sys.exit(1)

    if args.module:
        modules = [args.module]
    else:
        modules = None

    import_results = test_module_imports(modules)

    success = all(import_results.values())

    if success:
        logger.info("All module imports successful")

        # Test string handling in modules that imported successfully
        string_test_results = {}
        for module_name, imported in import_results.items():
            if imported:
                string_test_results[module_name] = test_cython_string_handling(
                    module_name
                )

        logger.info(f"String handling test results: {string_test_results}")
    else:
        logger.error("Some modules failed to import")

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()

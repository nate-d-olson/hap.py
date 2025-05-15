#!/usr/bin/env python3
"""
Script to check Python 3 migration status and run tests for converted modules
"""

import argparse
import importlib
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("py3-test")


def setup_environment(build_dir: Optional[str] = None) -> bool:
    """Set up the Python environment

    Args:
        build_dir: Path to the build directory

    Returns:
        True if successful
    """
    # Add source directory to path
    script_dir = Path(__file__).parent.resolve()
    src_py_dir = script_dir / "src" / "python"

    if src_py_dir.exists():
        sys.path.insert(0, str(src_py_dir))
        logger.info(f"Added {src_py_dir} to Python path")

    # Add build directory if specified
    if build_dir:
        build_path = Path(build_dir).resolve()
        if not build_path.exists():
            logger.error(f"Build directory {build_path} does not exist")
            return False

        # Add lib paths
        python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
        lib_paths = [
            build_path / "lib" / f"python{python_version}",
            build_path / "lib" / "python",
            build_path / "lib",
        ]

        for lib_path in lib_paths:
            if lib_path.exists():
                sys.path.insert(0, str(lib_path))
                logger.info(f"Added {lib_path} to Python path")

    return True


def find_py3_modules() -> Dict[str, List[Path]]:
    """Find all Python 3 modules in the source tree

    Returns:
        Dictionary mapping package names to lists of Python 3 modules
    """
    script_dir = Path(__file__).parent.resolve()
    src_py_dir = script_dir / "src" / "python"

    # Group modules by package
    modules_by_package = {}

    if src_py_dir.exists():
        # Find Python modules with .py3 extension
        for py3_file in src_py_dir.glob("**/*.py3"):
            rel_path = py3_file.relative_to(src_py_dir)
            package = rel_path.parent.name if len(rel_path.parts) > 1 else ""

            if package not in modules_by_package:
                modules_by_package[package] = []

            modules_by_package[package].append(py3_file)

    return modules_by_package


def test_import_module(module_path: Path) -> Tuple[bool, str]:
    """Test importing a Python 3 module

    Args:
        module_path: Path to the module file

    Returns:
        Tuple of (success, error_message)
    """
    try:
        # Get the module name
        rel_path = module_path.relative_to(
            Path(__file__).parent.resolve() / "src" / "python"
        )
        module_name = (
            str(rel_path.parent / rel_path.stem).replace("/", ".").replace("\\", ".")
        )

        # If module ends with .py3, remove that suffix
        if module_name.endswith(".py3"):
            module_name = module_name[:-4]

        # Create a spec for the module
        spec = importlib.util.spec_from_file_location(module_name, str(module_path))
        if spec is None or spec.loader is None:
            return False, f"Failed to create module spec for {module_path}"

        # Create the module
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module

        # Execute the module
        spec.loader.exec_module(module)

        return True, "Successfully imported module"

    except Exception as e:
        return False, f"Error importing module {module_path}: {str(e)}"


def run_tests(modules_by_package: Dict[str, List[Path]]) -> Dict[str, Dict[str, bool]]:
    """Run tests for all Python 3 modules

    Args:
        modules_by_package: Dictionary mapping package names to lists of Python 3 modules

    Returns:
        Dictionary with test results by package and module
    """
    results = {}

    for package, modules in modules_by_package.items():
        results[package] = {}

        for module_path in modules:
            module_name = module_path.stem
            if module_name.endswith(".py3"):
                module_name = module_name[:-4]

            logger.info(f"Testing module: {package}/{module_name}")
            success, message = test_import_module(module_path)

            results[package][module_name] = success
            if success:
                logger.info(f"  ✓ Module {package}/{module_name} imported successfully")
            else:
                logger.error(
                    f"  ✗ Module {package}/{module_name} failed to import: {message}"
                )

    return results


def run_cpp_tests(build_dir: str) -> bool:
    """Run C++ integration tests

    Args:
        build_dir: Path to the build directory

    Returns:
        True if all tests pass
    """
    try:
        logger.info("Running C++ integration tests...")

        # Use the enhanced validation script
        result = subprocess.run(
            [
                sys.executable,
                "validate_cpp_integration_enhanced.py",
                "--build-dir",
                build_dir,
                "--use-mock",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

        if result.returncode == 0:
            logger.info("C++ integration tests PASSED")
            return True
        else:
            logger.error(f"C++ integration tests FAILED: {result.stderr}")
            return False

    except Exception as e:
        logger.error(f"Error running C++ integration tests: {e}")
        return False


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Test Python 3 modules")
    parser.add_argument("--build-dir", help="Path to the build directory")
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose output"
    )
    parser.add_argument(
        "--test-cpp", action="store_true", help="Run C++ integration tests"
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Set up environment
    if not setup_environment(args.build_dir):
        return 1

    # Find and test Python 3 modules
    logger.info("Finding Python 3 modules...")
    py3_modules = find_py3_modules()

    if not py3_modules:
        logger.error("No Python 3 modules found")
        return 1

    for package, modules in py3_modules.items():
        package_name = package if package else "root"
        logger.info(f"Found {len(modules)} module(s) in package '{package_name}'")

    # Run tests
    logger.info("Running tests for Python 3 modules...")
    results = run_tests(py3_modules)

    # Print summary
    success_count = 0
    fail_count = 0

    for package, module_results in results.items():
        for module, success in module_results.items():
            if success:
                success_count += 1
            else:
                fail_count += 1

    logger.info(f"Test Summary: {success_count} passed, {fail_count} failed")

    # Run C++ integration tests if requested
    if args.test_cpp and args.build_dir:
        cpp_success = run_cpp_tests(args.build_dir)
        if not cpp_success:
            return 1

    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

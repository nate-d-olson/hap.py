"""
Modern dependency handling for hap.py.

This module provides utility functions for checking and using modern Python libraries
that replace custom C++ implementations in the original codebase.
"""

import importlib
import logging
import os
import sys
from typing import Dict, Optional, Tuple

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("happy.dependencies")

# Required package versions
REQUIRED_VERSIONS = {
    "pysam": "0.15.0",
    "pybedtools": "0.8.0",
    "numpy": "1.16.0",
    "pandas": "0.24.0",
}


def version_to_tuple(version_str: str) -> Tuple[int, ...]:
    """Convert version string to tuple for comparison."""
    return tuple(map(int, version_str.split(".")))


def check_package(package_name: str, min_version: str) -> bool:
    """
    Check if a package is installed with minimum required version.

    Args:
        package_name: Name of the package to check
        min_version: Minimum version required

    Returns:
        bool: True if package meets requirements, False otherwise
    """
    try:
        package = importlib.import_module(package_name)
        if not hasattr(package, "__version__"):
            logger.warning(f"Cannot determine version for {package_name}")
            return True  # Assume it's compatible

        current_version = package.__version__
        min_version_tuple = version_to_tuple(min_version)
        current_version_tuple = version_to_tuple(current_version)

        if current_version_tuple >= min_version_tuple:
            logger.debug(
                f"{package_name} {current_version} meets minimum requirement {min_version}"
            )
            return True
        else:
            logger.warning(
                f"{package_name} {current_version} is below minimum requirement {min_version}"
            )
            return False
    except ImportError:
        logger.warning(f"{package_name} is not installed")
        return False


def check_dependencies(required_only: bool = False) -> Dict[str, bool]:
    """
    Check all required dependencies.

    Args:
        required_only: If True, only check packages that are required

    Returns:
        Dict[str, bool]: Dictionary with package names and status
    """
    results = {}
    for package_name, min_version in REQUIRED_VERSIONS.items():
        results[package_name] = check_package(package_name, min_version)

    return results


def using_vcfeval() -> bool:
    """
    Check if vcfeval is available and configured.

    Returns:
        bool: True if vcfeval is available
    """
    # Check environment variable for vcfeval path
    vcfeval_path = os.environ.get("VCFEVAL_PATH")
    if vcfeval_path and os.path.exists(vcfeval_path):
        return True

    # Check for rtg tools in path
    rtg_path = os.environ.get("RTG_PATH")
    if rtg_path and os.path.exists(rtg_path):
        return True

    return False


def get_alternative_implementation(module_name: str) -> Optional[str]:
    """
    Get alternative implementation module name if available.

    Args:
        module_name: The original module name

    Returns:
        Optional[str]: Name of alternative module or None
    """
    # Map of original modules to their modern alternatives
    alternatives = {
        "happy.vcfutils": "happy.modern.vcfutils",
        "happy.bedutils": "happy.modern.bedutils",
        "happy.bamutils": "happy.modern.bamutils",
    }

    return alternatives.get(module_name)


if __name__ == "__main__":
    # When run directly, report on dependency status
    print("Checking hap.py dependencies:")
    results = check_dependencies()

    all_ok = all(results.values())

    if all_ok:
        print("✓ All dependencies satisfied")
        sys.exit(0)
    else:
        print("✗ Some dependencies are missing or outdated")
        for package, status in results.items():
            if not status:
                print(f"  - {package} >= {REQUIRED_VERSIONS[package]} is required")
        sys.exit(1)

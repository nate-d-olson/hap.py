#!/usr/bin/env python3
"""
Environment verification script for hap.py development.
Run this script to verify your development environment is correctly set up.
"""

import shutil
import subprocess
import sys
from pathlib import Path


def check_python_environment():
    """Check Python version and environment."""
    print("=== Python Environment Check ===")
    print(f"Python executable: {sys.executable}")
    print(f"Python version: {sys.version}")

    # Check if we're in the right environment
    if "happy-dev" in sys.executable:
        print("✅ Running in happy-dev environment")
    else:
        print("❌ Not running in happy-dev environment")
        print("Please run: micromamba activate happy-dev")
        return False

    # Check Python version
    if sys.version_info >= (3, 11):
        print("✅ Python version is 3.11+")
    else:
        print("❌ Python version should be 3.11+")
        return False

    return True


def check_package_installation():
    """Check if hap.py is properly installed."""
    print("\n=== Package Installation Check ===")
    try:
        import hap_py

        print(f"✅ hap_py module found: {hap_py.__file__}")

        # Check if it's the development version
        package_path = Path(hap_py.__file__).parent
        if "src/hap_py" in str(package_path):
            print("✅ Development installation detected")
        else:
            print("⚠️  Not a development installation")

        return True
    except ImportError as e:
        print(f"❌ Cannot import hap_py: {e}")
        print("Please run: pip install -e .")
        return False


def check_external_tools():
    """Check external tool availability."""
    print("\n=== External Tools Check ===")

    # Check RTG
    rtg_path = Path("external/rtg-tools-3.12.1/rtg")
    if rtg_path.exists():
        print(f"✅ RTG tools found: {rtg_path}")
    else:
        print(f"❌ RTG tools not found at: {rtg_path}")
        return False

    # Check other bioinformatics tools (optional)
    tools = ["bcftools", "samtools", "tabix"]
    for tool in tools:
        if shutil.which(tool):
            print(f"✅ {tool} found in PATH")
        else:
            print(f"⚠️  {tool} not found in PATH (optional)")

    return True


def check_pytest():
    """Check pytest availability."""
    print("\n=== Testing Framework Check ===")
    try:
        result = subprocess.run(
            ["pytest", "--version"], capture_output=True, text=True, check=True
        )
        print(f"✅ pytest available: {result.stdout.strip()}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"❌ pytest not available: {e}")
        return False


def main():
    """Run all verification checks."""
    print("hap.py Development Environment Verification\n")

    checks = [
        check_python_environment,
        check_package_installation,
        check_external_tools,
        check_pytest,
    ]

    all_passed = True
    for check in checks:
        passed = check()
        if not passed:
            all_passed = False

    print("\n=== Summary ===")
    if all_passed:
        print("✅ All checks passed! Your development environment is ready.")
        print("You can now run tests with: pytest tests/")
    else:
        print("❌ Some checks failed. Please fix the issues above.")
        print("Refer to the development setup instructions in the README.")

    return all_passed


if __name__ == "__main__":
    sys.exit(0 if main() else 1)

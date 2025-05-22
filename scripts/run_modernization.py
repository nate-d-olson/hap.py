#!/usr/bin/env python3
"""
Run all modernization steps for hap.py.

This script runs through all the modernization steps in sequence:
1. Setup project configuration
2. Fix Python 2 artifacts
3. Run pre-commit hooks
4. Generate progress report
5. Run validation tests
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def run_command(cmd, description, check=True):
    """
    Run a command and print its output.

    Args:
        cmd: Command to run
        description: Description of the command
        check: Whether to check the return code

    Returns:
        True if the command succeeded, False otherwise
    """
    print(f"\n=== {description} ===\n")
    try:
        result = subprocess.run(cmd, check=check, text=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        return False
    except Exception as e:
        print(f"Error: {e}")
        return False


def run_modernization(args):
    """
    Run all modernization steps.

    Args:
        args: Command-line arguments

    Returns:
        0 if all steps succeeded, non-zero otherwise
    """
    success = True
    project_root = Path(__file__).parent.parent
    scripts_dir = project_root / "scripts"

    # Step 1: Setup project configuration
    if args.setup or args.all:
        setup_cmd = [
            sys.executable,
            str(scripts_dir / "setup_project_config.py"),
            "--overwrite" if args.force else "",
        ]
        success = run_command(setup_cmd, "Setting up project configuration") and success

    # Step 2: Fix Python 2 artifacts
    if args.fix_py2 or args.all:
        fix_cmd = [
            sys.executable,
            str(scripts_dir / "fix_python2_artifacts.py"),
            "--dry-run" if args.dry_run else "",
        ]
        success = run_command(fix_cmd, "Fixing Python 2 artifacts") and success

    # Step 3: Run pre-commit hooks
    if args.pre_commit or args.all:
        # First install pre-commit if needed
        install_cmd = ["pip", "install", "pre-commit"]
        run_command(install_cmd, "Installing pre-commit")

        # Install the hooks
        hooks_cmd = ["pre-commit", "install"]
        run_command(hooks_cmd, "Installing pre-commit hooks")

        # Run the hooks
        precommit_cmd = ["pre-commit", "run", "--all-files"]
        success = (
            run_command(precommit_cmd, "Running pre-commit hooks", check=False)
            and success
        )

    # Step 4: Generate progress report
    if args.report or args.all:
        report_cmd = [
            sys.executable,
            str(scripts_dir / "generate_report.py"),
            "--detailed",
            f"--output={project_root}/PROGRESS_REPORT.md",
        ]
        success = run_command(report_cmd, "Generating progress report") and success

        if os.path.exists(f"{project_root}/PROGRESS_REPORT.md"):
            print(f"\nProgress report written to {project_root}/PROGRESS_REPORT.md")

    # Step 5: Run validation tests
    if args.validate or args.all:
        validate_cmd = [sys.executable, str(scripts_dir / "validate_replacements.py")]
        success = (
            run_command(validate_cmd, "Validating Python replacements") and success
        )

    # Step 6: Run unit tests
    if args.test or args.all:
        test_cmd = [
            "pytest",
            "tests/unit/test_python_blocksplit.py",
            "tests/unit/test_sequence_utils.py",
            "tests/unit/test_vcfcheck.py",
            "tests/unit/test_quantify.py",
            "tests/unit/test_python_preprocess.py",
            "-v",
        ]
        success = run_command(test_cmd, "Running unit tests") and success

    # Print summary
    print("\n=== Modernization Summary ===\n")
    if success:
        print("✅ All steps completed successfully!")
    else:
        print("⚠️ Some steps failed, please check the output above.")

    # Print next steps
    print("\nNext steps:")
    print("1. Review the progress report")
    print("2. Fix any issues reported by pre-commit or tests")
    print("3. Implement Python replacements for remaining C++ components")
    print("4. Create or update documentation")

    return 0 if success else 1


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Run hap.py modernization steps")
    parser.add_argument(
        "--all", action="store_true", help="Run all modernization steps"
    )
    parser.add_argument(
        "--setup", action="store_true", help="Setup project configuration"
    )
    parser.add_argument("--fix-py2", action="store_true", help="Fix Python 2 artifacts")
    parser.add_argument(
        "--pre-commit", action="store_true", help="Run pre-commit hooks"
    )
    parser.add_argument(
        "--report", action="store_true", help="Generate progress report"
    )
    parser.add_argument(
        "--validate", action="store_true", help="Validate Python replacements"
    )
    parser.add_argument("--test", action="store_true", help="Run unit tests")
    parser.add_argument(
        "--force", action="store_true", help="Force overwrite of existing files"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Don't actually make changes, just report",
    )

    args = parser.parse_args()

    # If no steps specified, run all
    if not (
        args.setup
        or args.fix_py2
        or args.pre_commit
        or args.report
        or args.validate
        or args.test
        or args.all
    ):
        args.all = True

    return run_modernization(args)


if __name__ == "__main__":
    sys.exit(main())

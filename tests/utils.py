"""
Utility functions for tests.
"""

import filecmp
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple


def get_project_root() -> Path:
    """Return the project root directory."""
    return Path(__file__).parent.parent


def get_build_dir() -> Path:
    """Return the build directory."""
    return get_project_root() / "build"


def get_bin_dir() -> Path:
    """Return the binary directory."""
    return get_build_dir() / "bin"


def get_example_dir() -> Path:
    """Return the example directory."""
    return get_project_root() / "example"


def run_command(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run a command and return the result.

    Args:
        cmd: Command to run as a list of strings
        check: If True, raise an exception on non-zero return code

    Returns:
        CompletedProcess object with return code, stdout, and stderr
    """
    return subprocess.run(cmd, check=check, capture_output=True, text=True)


def get_python_executable() -> str:
    """Return the path to the Python executable."""
    return sys.executable


def find_reference_file() -> Optional[str]:
    """Find a reference file for testing.

    First checks for HGREF environment variable, then looks in example directory.

    Returns:
        Path to reference file or None if not found
    """
    if "HGREF" in os.environ:
        ref_path = os.environ["HGREF"]
        if os.path.exists(ref_path):
            return ref_path

    # Try example directory
    example_ref = get_example_dir() / "chr21.fa"
    if example_ref.exists():
        return str(example_ref)

    return None


def check_tool_availability(tool_name: str) -> bool:
    """Check if a tool is available in the bin directory.

    Args:
        tool_name: Name of the tool to check

    Returns:
        True if the tool is available, False otherwise
    """
    tool_path = get_bin_dir() / tool_name
    return tool_path.exists()


def compare_files(file1: Path, file2: Path, ignore_comments: bool = False) -> bool:
    """Compare two files and return True if they are the same.

    Args:
        file1: First file to compare
        file2: Second file to compare
        ignore_comments: If True, ignore lines starting with #

    Returns:
        True if files are the same, False otherwise
    """
    if not file1.exists() or not file2.exists():
        return False

    if not ignore_comments:
        return filecmp.cmp(file1, file2)

    with open(file1, encoding="utf-8") as f1, open(file2, encoding="utf-8") as f2:
        lines1 = [line for line in f1.readlines() if not line.strip().startswith("#")]
        lines2 = [line for line in f2.readlines() if not line.strip().startswith("#")]
        return lines1 == lines2


def run_shell_command(command: str, cwd: Optional[Path] = None) -> Tuple[int, str, str]:
    """Run a shell command and return the exit code, stdout, and stderr.

    Args:
        command: Shell command to run
        cwd: Working directory

    Returns:
        Tuple of (exit_code, stdout, stderr)
    """
    result = subprocess.run(
        command,
        shell=True,
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    return result.returncode, result.stdout, result.stderr


def prepare_temp_vcf(source_vcf: Path, output_dir: Path) -> Path:
    """Create a temporary VCF file from a source VCF.

    This utility function helps when migrating shell tests that create temporary
    VCF files for processing.

    Args:
        source_vcf: Source VCF file
        output_dir: Directory to write the temporary file

    Returns:
        Path to the temporary VCF file
    """
    temp_vcf = output_dir / f"temp_{source_vcf.stem}.vcf"
    shutil.copy2(source_vcf, temp_vcf)
    return temp_vcf


def get_integration_test_name(test_file: Path) -> str:
    """Extract the test name from a test file path.

    Args:
        test_file: Path to the test file

    Returns:
        Test name
    """
    filename = test_file.name
    if filename.startswith("test_") and filename.endswith(".py"):
        return filename[5:-3]  # Remove 'test_' and '.py'
    return filename


def compare_summary_files(file1: Path, file2: Path, tolerance: float = 0.001) -> bool:
    """Compare two summary CSV files with a tolerance.

    This is a Python implementation of the compare_summaries.py script.

    Args:
        file1: First summary CSV file
        file2: Second summary CSV file
        tolerance: Tolerance for floating point comparisons

    Returns:
        True if files are equivalent within tolerance, False otherwise
    """
    import csv

    if not file1.exists() or not file2.exists():
        return False

    # Read the first file
    with open(file1, "r") as f:
        reader = csv.DictReader(f)
        data1 = list(reader)

    # Read the second file
    with open(file2, "r") as f:
        reader = csv.DictReader(f)
        data2 = list(reader)

    # Check if they have the same number of rows
    if len(data1) != len(data2):
        return False

    # Compare each row, allowing for small differences in floating point values
    for row1, row2 in zip(data1, data2):
        for key in row1.keys():
            if key not in row2:
                return False

            # Try to convert to float for numeric comparisons
            try:
                val1 = float(row1[key]) if row1[key] else 0.0
                val2 = float(row2[key]) if row2[key] else 0.0

                # If values are very small, use absolute difference
                if abs(val1) < tolerance and abs(val2) < tolerance:
                    if abs(val1 - val2) > tolerance:
                        return False
                # Otherwise use relative difference
                elif val2 != 0 and abs((val1 - val2) / val2) > tolerance:
                    return False
            except ValueError:
                # For non-numeric values, compare directly
                if row1[key] != row2[key]:
                    return False

    return True


def check_vcfeval_availability() -> bool:
    """Check if vcfeval is available in the Haplo module.

    This is a Python implementation of checking the version.py file for vcfeval.

    Returns:
        bool: True if vcfeval is available, False otherwise
    """
    try:
        import sys
        import os

        project_root = get_project_root()
        sys.path.insert(0, os.path.join(project_root, "src", "python"))

        # Try to import the version module
        try:
            from Haplo import version

            return getattr(version, "has_vcfeval", 0) != 0
        except (ImportError, AttributeError):
            # If we can't import or the attribute doesn't exist, try to parse the file
            version_py = os.path.join(
                project_root, "src", "python", "Haplo", "version.py"
            )
            if os.path.exists(version_py):
                with open(version_py, "r") as f:
                    for line in f:
                        if "has_vcfeval" in line:
                            return "1" in line
            return False
    except Exception:
        # If any error occurs, assume vcfeval is not available
        return False

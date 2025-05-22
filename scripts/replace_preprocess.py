#!/usr/bin/env python3
"""
Script to replace C++ preprocess with Python implementation.

This script updates the necessary files to use the Python implementation
of preprocess instead of the C++ implementation.
"""

import argparse
import os
import re
import sys


def update_partialcredit_py(path):
    """
    Update partialcredit.py to use the Python preprocess.

    Args:
        path: Path to partialcredit.py
    """
    if not os.path.exists(path):
        print(f"Error: {path} does not exist")
        return False

    with open(path, 'r') as f:
        content = f.read()

    # Replace the preprocessWrapper function
    new_preprocess_wrapper = """
def preprocessWrapper(
    file_and_location: Tuple[str, str], args: Dict[str, Any]
) -> Optional[str]:
    \"\"\"Process a VCF file with the Python preprocess tool.

    Args:
        file_and_location: Tuple of (filename, location_str)
        args: Arguments for preprocessing

    Returns:
        Path to the preprocessed output file or None if processing failed
    \"\"\"
    starttime = time.time()
    filename, location_str = file_and_location
    int_suffix = "bcf" if args["bcf"] else "vcf.gz"
    temp_file_path = None

    try:
        with tempfile.NamedTemporaryFile(
            delete=False, prefix=f"input.{location_str}", suffix=f".prep.{int_suffix}"
        ) as tf:
            temp_file_path = tf.name  # Store the file path for later use or cleanup

        # Create the PreprocessEngine
        from Haplo.python_preprocess import PreprocessEngine

        decompose_level = args["decompose"]
        left_shift = args["leftshift"]

        engine = PreprocessEngine(
            input_vcf=filename,
            reference_fasta=args["reference"],
            output_vcf=temp_file_path,
            decompose_level=decompose_level,
            left_shift=left_shift,
            regions=location_str if location_str else None,
            haploid_x=args.get("haploid_x", False),
            output_bcf=args.get("bcf", False)
        )

        # Process the file
        processed_file = engine.process()

        elapsed = time.time() - starttime
        logging.info(f"preprocess for {location_str} -- time taken {elapsed:.2f}")

        return processed_file
    except Exception as e:
        logging.error(f"Exception in preprocessWrapper for {location_str}: {str(e)}")
        # Clean up temp file if it exists
        if temp_file_path and os.path.exists(temp_file_path):
            try:
                os.unlink(temp_file_path)
            except:
                pass
        return None
"""

    # Use regex to replace the function
    pattern = r"def preprocessWrapper\([^)]*\).*?(?=def|$)"
    replacement = new_preprocess_wrapper
    new_content = re.sub(pattern, replacement, content, flags=re.DOTALL)

    # Write back
    with open(path, 'w') as f:
        f.write(new_content)

    print(f"Updated {path}")
    return True


def update_scripts(scripts_dir):
    """
    Update scripts to use the Python preprocess.

    Args:
        scripts_dir: Path to scripts directory
    """
    success = True

    # Add preprocess to scripts/generate_report.py
    report_path = os.path.join(scripts_dir, "generate_report.py")
    if os.path.exists(report_path):
        with open(report_path, 'r') as f:
            content = f.read()

        # Update component status
        if '"preprocess": {"status": "not-started"' in content:
            content = content.replace(
                '"preprocess": {"status": "not-started"',
                '"preprocess": {"status": "completed"'
            )

        with open(report_path, 'w') as f:
            f.write(content)

        print(f"Updated {report_path}")
    else:
        print(f"Warning: {report_path} does not exist")
        success = False

    return success


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Replace C++ preprocess with Python implementation"
    )
    parser.add_argument(
        "--project-root", type=str, default=None,
        help="Path to hap.py project root (default: auto-detect)"
    )

    args = parser.parse_args()

    # Determine project root
    if args.project_root:
        project_root = args.project_root
    else:
        # Try to auto-detect based on script location
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.abspath(os.path.join(script_dir, ".."))

    print(f"Using project root: {project_root}")

    # Paths
    partialcredit_path = os.path.join(
        project_root, "src", "python", "Haplo", "partialcredit.py"
    )
    scripts_dir = os.path.join(project_root, "scripts")

    # Update files
    success = True
    success = update_partialcredit_py(partialcredit_path) and success
    success = update_scripts(scripts_dir) and success

    # Print summary
    print("\n=== Summary ===")
    if success:
        print("✅ All updates completed successfully")
    else:
        print("⚠️ Some updates failed, please check the output above")

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())

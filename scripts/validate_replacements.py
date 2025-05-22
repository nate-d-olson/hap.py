#!/usr/bin/env python3
"""
Validate Python replacements for C++ components.

This script:
1. Identifies C++ components that have Python replacements
2. Runs both implementations on the same input
3. Compares outputs to ensure functional equivalence
4. Reports performance differences
"""

import argparse
import hashlib
import json
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Tuple


class ComponentValidator:
    """Validate Python replacements for C++ components."""

    def __init__(self, project_root: Path):
        """
        Initialize the component validator.

        Args:
            project_root: Root directory of the project
        """
        self.project_root = project_root
        self.example_dir = project_root / "example"
        self.build_dir = project_root / "build"
        self.python_dir = project_root / "src" / "python"

        # Components to validate (name: {cpp_path, py_path, test_args})
        self.components = {
            "blocksplit": {
                "cpp_path": self.build_dir / "bin" / "blocksplit",
                "py_path": self.python_dir / "Haplo" / "python_blocksplit.py",
                "test_args": [
                    "--input",
                    str(self.example_dir / "hc.vcf.gz"),
                    "--output",
                    "{outdir}/blocks.bed",
                    "--block-size",
                    "1000",
                ],
            },
            "vcfcheck": {
                "cpp_path": self.build_dir / "bin" / "vcfcheck",
                "py_path": self.python_dir / "Haplo" / "python_vcfcheck.py",
                "test_args": [
                    "--input",
                    str(self.example_dir / "hc.vcf.gz"),
                    "--reference",
                    str(self.example_dir / "chr21.fa"),
                    "--output",
                    "{outdir}/vcfcheck.txt",
                ],
            },
            # Add more components as they are implemented
        }

        self.results = {}

    def validate_all(self, components: List[str] = None) -> Dict[str, Dict]:
        """
        Validate all components or a specific list of components.

        Args:
            components: List of component names to validate (default: all)

        Returns:
            Dictionary of validation results
        """
        if not components:
            components = list(self.components.keys())

        print(f"🔍 Validating {len(components)} components")

        for component in components:
            if component not in self.components:
                print(f"⚠️ Unknown component: {component}")
                continue

            print(f"\n=== Validating {component} ===")
            try:
                result = self.validate_component(component)
                self.results[component] = result
                self._print_result(component, result)
            except Exception as e:
                print(f"❌ Error validating {component}: {e}")
                self.results[component] = {"error": str(e), "success": False}

        return self.results

    def validate_component(self, component: str) -> Dict:
        """
        Validate a specific component.

        Args:
            component: Name of the component to validate

        Returns:
            Dictionary of validation results
        """
        component_info = self.components[component]
        cpp_path = component_info["cpp_path"]
        py_path = component_info["py_path"]

        # Check if both implementations exist
        if not cpp_path.exists():
            raise FileNotFoundError(f"C++ implementation not found: {cpp_path}")
        if not py_path.exists():
            raise FileNotFoundError(f"Python implementation not found: {py_path}")

        # Create temporary directories for outputs
        with tempfile.TemporaryDirectory() as cpp_outdir, tempfile.TemporaryDirectory() as py_outdir:

            # Prepare arguments
            cpp_args = [
                str(arg).replace("{outdir}", cpp_outdir)
                for arg in component_info["test_args"]
            ]
            py_args = [
                str(arg).replace("{outdir}", py_outdir)
                for arg in component_info["test_args"]
            ]

            # Run C++ implementation
            cpp_start = time.time()
            cpp_result = self._run_cpp(component, cpp_path, cpp_args)
            cpp_duration = time.time() - cpp_start

            # Run Python implementation
            py_start = time.time()
            py_result = self._run_python(component, py_path, py_args)
            py_duration = time.time() - py_start

            # Compare outputs
            outputs_match, diff_details = self._compare_outputs(cpp_outdir, py_outdir)

            # Calculate performance difference
            perf_ratio = (
                py_duration / cpp_duration if cpp_duration > 0 else float("inf")
            )

            return {
                "success": outputs_match
                and cpp_result["return_code"] == py_result["return_code"],
                "outputs_match": outputs_match,
                "diff_details": diff_details,
                "cpp": {
                    "duration": cpp_duration,
                    "return_code": cpp_result["return_code"],
                    "output": (
                        cpp_result["output"][:500] + "..."
                        if len(cpp_result["output"]) > 500
                        else cpp_result["output"]
                    ),
                },
                "python": {
                    "duration": py_duration,
                    "return_code": py_result["return_code"],
                    "output": (
                        py_result["output"][:500] + "..."
                        if len(py_result["output"]) > 500
                        else py_result["output"]
                    ),
                },
                "performance": {
                    "ratio": perf_ratio,
                    "acceptable": perf_ratio < 5.0,  # Python can be up to 5x slower
                },
            }

    def _run_cpp(self, component: str, executable: Path, args: List[str]) -> Dict:
        """
        Run a C++ component.

        Args:
            component: Component name
            executable: Path to the executable
            args: Arguments to pass to the executable

        Returns:
            Dictionary with return code and output
        """
        print(f"Running C++ {component}...")
        cmd = [str(executable)] + args
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            return {
                "return_code": result.returncode,
                "output": result.stdout + result.stderr,
            }
        except Exception as e:
            return {"return_code": -1, "output": str(e)}

    def _run_python(self, component: str, script: Path, args: List[str]) -> Dict:
        """
        Run a Python component.

        Args:
            component: Component name
            script: Path to the Python script
            args: Arguments to pass to the script

        Returns:
            Dictionary with return code and output
        """
        print(f"Running Python {component}...")
        cmd = [sys.executable, str(script)] + args
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            return {
                "return_code": result.returncode,
                "output": result.stdout + result.stderr,
            }
        except Exception as e:
            return {"return_code": -1, "output": str(e)}

    def _compare_outputs(self, cpp_dir: str, py_dir: str) -> Tuple[bool, str]:
        """
        Compare outputs from C++ and Python implementations.

        Args:
            cpp_dir: Directory containing C++ outputs
            py_dir: Directory containing Python outputs

        Returns:
            Tuple of (match, details)
        """
        cpp_files = sorted(Path(cpp_dir).glob("*"))
        py_files = sorted(Path(py_dir).glob("*"))

        if len(cpp_files) != len(py_files):
            return (
                False,
                f"Different number of output files: C++={len(cpp_files)}, Python={len(py_files)}",
            )

        all_match = True
        details = []

        for cpp_file, py_file in zip(cpp_files, py_files):
            if cpp_file.name != py_file.name:
                details.append(f"Filename mismatch: {cpp_file.name} vs {py_file.name}")
                all_match = False
                continue

            # Compare file content
            cpp_hash = self._file_hash(cpp_file)
            py_hash = self._file_hash(py_file)

            if cpp_hash != py_hash:
                details.append(f"Content mismatch in {cpp_file.name}")
                # Add detailed diff for small text files
                if cpp_file.stat().st_size < 1024 * 1024 and cpp_file.suffix in [
                    ".txt",
                    ".bed",
                    ".vcf",
                    ".csv",
                ]:
                    try:
                        diff_cmd = ["diff", "-u", str(cpp_file), str(py_file)]
                        diff_result = subprocess.run(
                            diff_cmd, capture_output=True, text=True
                        )
                        if diff_result.stdout:
                            details.append(
                                f"Diff for {cpp_file.name}:\n{diff_result.stdout[:500]}..."
                            )
                    except Exception:
                        pass
                all_match = False

        if all_match:
            return True, "All output files match exactly"
        else:
            return False, "\n".join(details)

    def _file_hash(self, file_path: Path) -> str:
        """
        Compute a hash of a file.

        Args:
            file_path: Path to the file

        Returns:
            SHA-256 hash of the file
        """
        h = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                h.update(chunk)
        return h.hexdigest()

    def _print_result(self, component: str, result: Dict):
        """
        Print validation result.

        Args:
            component: Component name
            result: Validation result
        """
        if result["success"]:
            print(f"✅ {component} validation successful!")
        else:
            print(f"❌ {component} validation failed!")

        if "outputs_match" in result:
            if result["outputs_match"]:
                print("✅ Outputs match exactly")
            else:
                print("❌ Outputs differ:")
                print(result["diff_details"])

        if "performance" in result:
            perf = result["performance"]
            if perf["acceptable"]:
                print(
                    f"✅ Performance acceptable: Python is {perf['ratio']:.2f}x slower than C++"
                )
            else:
                print(
                    f"⚠️ Performance concern: Python is {perf['ratio']:.2f}x slower than C++"
                )

        print("\nSummary:")
        print(f"C++ execution time: {result['cpp']['duration']:.3f}s")
        print(f"Python execution time: {result['python']['duration']:.3f}s")
        if result["cpp"]["return_code"] != 0 or result["python"]["return_code"] != 0:
            print(
                f"Return codes: C++={result['cpp']['return_code']}, Python={result['python']['return_code']}"
            )


def main():
    """Main function."""

    parser = argparse.ArgumentParser(
        description="Validate Python replacements for C++ components"
    )
    parser.add_argument(
        "components", nargs="*", help="Components to validate (default: all available)"
    )
    parser.add_argument("--json", action="store_true", help="Output results as JSON")
    parser.add_argument("--output", type=str, help="Output file for JSON results")

    args = parser.parse_args()

    # Get project root directory (parent of the script)
    project_root = Path(__file__).parent.parent

    # Create and run the validator
    validator = ComponentValidator(project_root)
    results = validator.validate_all(args.components)

    # Output results
    if args.json:
        # Convert to JSON-serializable format
        json_results = {}
        for component, result in results.items():
            json_result = {k: v for k, v in result.items()}
            # Handle non-serializable objects
            if "diff_details" in json_result and isinstance(
                json_result["diff_details"], str
            ):
                json_result["diff_details"] = json_result["diff_details"][
                    :1000
                ]  # Truncate long diff
            json_results[component] = json_result

        if args.output:
            with open(args.output, "w") as f:
                json.dump(json_results, f, indent=2)
            print(f"Results written to {args.output}")
        else:
            print(json.dumps(json_results, indent=2))

    # Return success if all validations passed
    return 0 if all(r.get("success", False) for r in results.values()) else 1


if __name__ == "__main__":
    import sys

    sys.exit(main())

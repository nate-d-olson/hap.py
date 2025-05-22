#!/usr/bin/env python3
"""
Progress tracking script for hap.py modernization.

This script monitors the progress of Python 3 modernization efforts,
tracks C++ dependency reduction, and generates progress reports.
"""

import json
import subprocess
from pathlib import Path
from typing import Dict


class ModernizationTracker:
    """Track progress of modernization efforts."""

    def __init__(self, config_file: str = "modernization_progress.json"):
        self.config_file = Path(config_file)
        self.project_root = Path(__file__).parent.parent
        self.load_progress()

    def load_progress(self):
        """Load existing progress data."""
        if self.config_file.exists():
            with open(self.config_file) as f:
                self.progress_data = json.load(f)
        else:
            self.progress_data = {
                "phases": {
                    "phase1": {"completed": False, "tasks": []},
                    "phase2": {"completed": False, "tasks": []},
                    "phase3": {"completed": False, "tasks": []},
                    "phase4": {"completed": False, "tasks": []},
                },
                "metrics": {},
                "cpp_components": {},
            }

    def save_progress(self):
        """Save progress data to file."""
        with open(self.config_file, "w") as f:
            json.dump(self.progress_data, f, indent=2)

    def check_code_quality(self) -> Dict[str, any]:
        """Check code quality metrics."""
        results = {}

        # Check if pre-commit is installed and configured
        try:
            result = subprocess.run(
                ["pre-commit", "--version"], capture_output=True, text=True, check=True
            )
            results["pre_commit_installed"] = True
        except (subprocess.CalledProcessError, FileNotFoundError):
            results["pre_commit_installed"] = False

        # Check for Python 2 artifacts
        results["python2_artifacts"] = self.check_python2_artifacts()

        # Check type coverage (simplified)
        results["type_coverage"] = self.check_type_coverage()

        # Check test files
        results["test_files"] = self.count_test_files()

        return results

    def check_python2_artifacts(self) -> Dict[str, int]:
        """Check for remaining Python 2 artifacts."""
        artifacts = {
            "print_statements": 0,
            "unicode_literals": 0,
            "old_super": 0,
            "xrange": 0,
        }

        python_files = list(self.project_root.glob("src/python/**/*.py"))

        for file_path in python_files:
            try:
                with open(file_path, encoding="utf-8") as f:
                    content = f.read()

                # Check for print statements (not function calls)
                lines = content.split("\n")
                for line in lines:
                    stripped = line.strip()
                    if stripped.startswith("print ") and not stripped.startswith(
                        "print("
                    ):
                        artifacts["print_statements"] += 1
                    if "unicode_literals" in line:
                        artifacts["unicode_literals"] += 1
                    if "super(" in line and "super()" not in line:
                        artifacts["old_super"] += 1
                    if "xrange" in line:
                        artifacts["xrange"] += 1

            except Exception as e:
                print(f"Warning: Could not read {file_path}: {e}")

        return artifacts

    def check_type_coverage(self) -> int:
        """Estimate type hint coverage."""
        python_files = list(self.project_root.glob("src/python/**/*.py"))
        total_functions = 0
        typed_functions = 0

        for file_path in python_files:
            try:
                with open(file_path, encoding="utf-8") as f:
                    content = f.read()
                    lines = content.split("\n")

                    for line in lines:
                        stripped = line.strip()
                        if stripped.startswith("def ") and ":" in stripped:
                            total_functions += 1
                            if "->" in stripped or ": " in stripped:
                                typed_functions += 1

            except Exception as e:
                print(f"Warning: Could not read {file_path}: {e}")

        if total_functions == 0:
            return 0
        return int((typed_functions / total_functions) * 100)

    def count_test_files(self) -> Dict[str, int]:
        """Count test files."""
        test_counts = {
            "pytest_files": len(list(self.project_root.glob("**/test_*.py"))),
            "shell_test_files": len(list(self.project_root.glob("src/sh/*test*.sh"))),
            "total_tests": 0,
        }
        test_counts["total_tests"] = (
            test_counts["pytest_files"] + test_counts["shell_test_files"]
        )
        return test_counts

    def check_cpp_dependencies(self) -> Dict[str, Dict[str, any]]:
        """Count and analyze C++ dependencies."""
        cpp_dir = self.project_root / "src" / "c++"
        if not cpp_dir.exists():
            return {}

        components = {}

        # Main C++ executables to track
        main_executables = [
            "blocksplit",
            "quantify",
            "vcfcheck",
            "preprocess",
            "xcmp",
            "scmp",
            "hapcmp",
            "multimerge",
        ]

        for component in main_executables:
            cpp_files = list(cpp_dir.rglob(f"*{component}*.cpp"))
            header_files = list(cpp_dir.rglob(f"*{component}*.h"))

            components[component] = {
                "cpp_files": len(cpp_files),
                "header_files": len(header_files),
                "total_files": len(cpp_files) + len(header_files),
                "has_python_replacement": self.check_python_replacement(component),
                "file_paths": [str(f) for f in cpp_files + header_files],
            }

        return components

    def check_python_replacement(self, component: str) -> bool:
        """Check if a Python replacement exists for a C++ component."""
        python_dir = self.project_root / "src" / "python"

        # Look for Python files that might replace the C++ component
        python_files = list(python_dir.rglob(f"*{component}*.py"))

        # Also check for mentions in Python files
        for py_file in python_dir.rglob("*.py"):
            try:
                with open(py_file, encoding="utf-8") as f:
                    content = f.read()
                    if (
                        f"python_{component}" in content.lower()
                        or f"{component}_python" in content.lower()
                    ):
                        return True
            except Exception:
                continue

        return len(python_files) > 0

    def check_dependency_packages(self) -> Dict[str, bool]:
        """Check if modern Python packages are being used."""
        pyproject_file = self.project_root / "pyproject.toml"
        dependencies = {}

        if pyproject_file.exists():
            try:
                with open(pyproject_file) as f:
                    content = f.read()

                dependencies = {
                    "pysam": "pysam" in content,
                    "biopython": "biopython" in content.lower()
                    or "bio-python" in content.lower(),
                    "numpy": "numpy" in content,
                    "scipy": "scipy" in content,
                    "pandas": "pandas" in content,
                    "pytest": "pytest" in content,
                    "black": "black" in content,
                    "ruff": "ruff" in content,
                }
            except Exception as e:
                print(f"Warning: Could not read pyproject.toml: {e}")

        return dependencies

    def generate_report(self, detailed: bool = False) -> str:
        """Generate modernization progress report."""
        quality = self.check_code_quality()
        cpp_deps = self.check_cpp_dependencies()
        packages = self.check_dependency_packages()

        report = "# hap.py Modernization Progress Report\n\n"

        # Overall status
        report += "## Overall Status\n"
        total_cpp_components = len(cpp_deps)
        migrated_components = sum(
            1 for comp in cpp_deps.values() if comp["has_python_replacement"]
        )
        migration_percentage = (
            (migrated_components / total_cpp_components * 100)
            if total_cpp_components > 0
            else 0
        )

        report += f"- **Migration Progress**: {migration_percentage:.1f}% ({migrated_components}/{total_cpp_components} components)\n"
        report += f"- **Pre-commit Setup**: {'✅ Configured' if quality['pre_commit_installed'] else '❌ Not configured'}\n"
        report += f"- **Type Coverage**: {quality['type_coverage']}%\n"
        report += f"- **Test Files**: {quality['test_files']['pytest_files']} pytest, {quality['test_files']['shell_test_files']} shell scripts\n\n"

        # Code quality status
        report += "## Code Quality Status\n"
        artifacts = quality["python2_artifacts"]
        total_artifacts = sum(artifacts.values())

        if total_artifacts == 0:
            report += "- **Python 2 Artifacts**: ✅ None found\n"
        else:
            report += f"- **Python 2 Artifacts**: ❌ {total_artifacts} issues found\n"
            for artifact, count in artifacts.items():
                if count > 0:
                    report += f"  - {artifact}: {count}\n"

        report += "\n"

        # Package dependencies
        report += "## Modern Package Adoption\n"
        for package, adopted in packages.items():
            status = "✅ Adopted" if adopted else "❌ Missing"
            report += f"- **{package}**: {status}\n"
        report += "\n"

        # C++ component status
        report += "## C++ Component Migration Status\n"
        for component, details in cpp_deps.items():
            if details["has_python_replacement"]:
                status = "✅ Python replacement exists"
            elif details["total_files"] == 0:
                status = "ℹ️ No files found"
            else:
                status = f"🔄 {details['total_files']} C++ files remaining"

            report += f"- **{component}**: {status}\n"

            if detailed and details["total_files"] > 0:
                report += f"  - C++ files: {details['cpp_files']}\n"
                report += f"  - Header files: {details['header_files']}\n"

        report += "\n"

        # Priority recommendations
        report += "## Next Priority Actions\n"

        if total_artifacts > 0:
            report += "1. **Immediate**: Fix remaining Python 2 artifacts\n"

        if not quality["pre_commit_installed"]:
            report += "2. **Setup**: Install and configure pre-commit hooks\n"

        # Find components that should be prioritized
        high_priority_components = []
        for component, details in cpp_deps.items():
            if not details["has_python_replacement"] and details["total_files"] > 0:
                # Prioritize by component importance
                if component in ["blocksplit", "quantify", "vcfcheck", "preprocess"]:
                    high_priority_components.append(component)

        if high_priority_components:
            report += f"3. **Migration**: Start with high-priority components: {', '.join(high_priority_components)}\n"

        if not packages.get("biopython", False):
            report += "4. **Dependencies**: Add BioPython for sequence processing\n"

        return report

    def update_metrics(self):
        """Update and save current metrics."""
        self.progress_data["metrics"] = {
            "code_quality": self.check_code_quality(),
            "cpp_dependencies": self.check_cpp_dependencies(),
            "packages": self.check_dependency_packages(),
            "timestamp": subprocess.check_output(["date", "+%Y-%m-%d %H:%M:%S"])
            .decode()
            .strip(),
        }
        self.save_progress()


def main():
    """Main function to run the tracker."""
    import argparse

    parser = argparse.ArgumentParser(description="Track hap.py modernization progress")
    parser.add_argument(
        "--detailed", action="store_true", help="Generate detailed report"
    )
    parser.add_argument(
        "--update", action="store_true", help="Update metrics and save progress"
    )
    parser.add_argument("--json", action="store_true", help="Output raw JSON data")

    args = parser.parse_args()

    tracker = ModernizationTracker()

    if args.update:
        tracker.update_metrics()
        print("✅ Progress metrics updated")

    if args.json:
        # Update metrics first
        tracker.update_metrics()
        print(json.dumps(tracker.progress_data, indent=2))
    else:
        report = tracker.generate_report(detailed=args.detailed)
        print(report)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Generate a modernization progress report for the hap.py project.

This script analyzes the codebase, tracks migration progress from C++ to Python,
and generates a detailed report in Markdown format.
"""

import argparse
import datetime
from pathlib import Path
from typing import Dict


class ModernizationReporter:
    """Generate a modernization progress report."""

    def __init__(self, project_root: Path):
        """
        Initialize the modernization reporter.

        Args:
            project_root: Root directory of the project
        """
        self.project_root = project_root
        self.cpp_dir = project_root / "src" / "c++"
        self.python_dir = project_root / "src" / "python"
        self.scripts_dir = project_root / "scripts"
        self.test_dir = project_root / "tests"

        # Components to track
        self.cpp_components = {
            "blocksplit": {"category": "core", "priority": "high"},
            "quantify": {"category": "core", "priority": "high"},
            "vcfcheck": {"category": "validation", "priority": "high"},
            "preprocess": {"category": "core", "priority": "high"},
            "xcmp": {"category": "comparison", "priority": "medium"},
            "scmp": {"category": "comparison", "priority": "medium"},
            "hapcmp": {"category": "comparison", "priority": "medium"},
            "multimerge": {"category": "utility", "priority": "low"},
        }

    def generate_report(self, detailed: bool = False) -> str:
        """
        Generate a modernization progress report.

        Args:
            detailed: Whether to include detailed information

        Returns:
            Markdown-formatted report
        """
        # Get metrics
        cpp_stats = self._analyze_cpp_components()
        python_stats = self._analyze_python_replacement_status()
        test_stats = self._analyze_test_coverage()
        quality_stats = self._analyze_code_quality()
        dependency_stats = self._analyze_dependencies()

        # Calculate overall progress
        total_components = len(self.cpp_components)
        migrated_components = sum(
            1 for c, s in python_stats.items() if s["has_replacement"]
        )
        progress_percentage = (
            (migrated_components / total_components) * 100
            if total_components > 0
            else 0
        )

        # Generate report
        report = f"""# hap.py Modernization Progress Report

*Generated on: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}*

## Overall Progress

**Migration Status**: {progress_percentage:.1f}% complete ({migrated_components}/{total_components} components)

![Progress](https://progress-bar.dev/{int(progress_percentage)}/)

## Migration Status by Component

| Component | Category | Priority | Status | Python Replacement |
|-----------|----------|----------|--------|-------------------|
"""

        # Add component rows
        for component, info in sorted(self.cpp_components.items()):
            category = info["category"]
            priority = info["priority"]

            cpp_info = cpp_stats.get(component, {"files": 0, "loc": 0})
            py_info = python_stats.get(
                component, {"has_replacement": False, "files": 0, "loc": 0}
            )

            if py_info["has_replacement"]:
                status = "✅ Migrated"
                py_details = f"`{py_info['files']} files, {py_info['loc']} lines`"
            else:
                status = "❌ Pending"
                cpp_remaining = (
                    f"`{cpp_info['files']} C++ files, {cpp_info['loc']} lines`"
                )
                py_details = cpp_remaining

            report += (
                f"| {component} | {category} | {priority} | {status} | {py_details} |\n"
            )

        # Code quality
        report += """
## Code Quality Status

| Metric | Value |
|--------|-------|
"""

        for metric, value in quality_stats.items():
            report += f"| {metric} | {value} |\n"

        # Dependency status
        report += """
## Dependency Status

| Package | Status | Utilization |
|---------|--------|-------------|
"""

        for pkg, info in dependency_stats.items():
            status = "✅ Adopted" if info["installed"] else "❌ Missing"
            util = (
                f"{info['utilization']}%" if info["utilization"] is not None else "N/A"
            )
            report += f"| {pkg} | {status} | {util} |\n"

        # Test coverage
        report += """
## Test Coverage

| Type | Count | Migration Status |
|------|-------|------------------|
"""

        for test_type, info in test_stats.items():
            status = f"{info['migrated']}/{info['total']} migrated"
            report += f"| {test_type} | {info['total']} | {status} |\n"

        # Priority recommendations
        report += """
## Next Steps

"""

        # Identify next components to migrate
        unmigrated = [
            c
            for c, s in python_stats.items()
            if not s["has_replacement"] and self.cpp_components[c]["priority"] == "high"
        ]

        if unmigrated:
            report += (
                f"1. **Priority Components to Migrate**: {', '.join(unmigrated)}\n"
            )

        # Identify quality issues
        quality_issues = []
        if quality_stats.get("Python 2 Artifacts", 0) > 0:
            quality_issues.append("Fix remaining Python 2 artifacts")
        if quality_stats.get("Type Coverage", "0%").rstrip("%") < "50":
            quality_issues.append("Improve type hint coverage")

        if quality_issues:
            report += f"2. **Code Quality**: {', '.join(quality_issues)}\n"

        # Missing dependencies
        missing_deps = [
            p
            for p, i in dependency_stats.items()
            if not i["installed"] and p in ["pysam", "biopython", "pytest"]
        ]
        if missing_deps:
            report += f"3. **Dependencies**: Add {', '.join(missing_deps)}\n"

        # Add detailed sections if requested
        if detailed:
            report += self._generate_detailed_sections(
                cpp_stats, python_stats, test_stats
            )

        return report

    def _analyze_cpp_components(self) -> Dict[str, Dict]:
        """
        Analyze C++ components in the codebase.

        Returns:
            Dictionary mapping component names to statistics
        """
        stats = {}

        for component in self.cpp_components:
            cpp_files = list(self.cpp_dir.glob(f"**/*{component}*.cpp")) + list(
                self.cpp_dir.glob(f"**/*{component}*.h")
            )

            loc = 0
            for file_path in cpp_files:
                try:
                    with open(file_path, encoding="utf-8") as f:
                        loc += sum(
                            1
                            for line in f
                            if line.strip() and not line.strip().startswith("//")
                        )
                except Exception:
                    pass

            stats[component] = {"files": len(cpp_files), "loc": loc}

        return stats

    def _analyze_python_replacement_status(self) -> Dict[str, Dict]:
        """
        Analyze Python replacement status for C++ components.

        Returns:
            Dictionary mapping component names to replacement status
        """
        stats = {}

        for component in self.cpp_components:
            # Look for direct replacements (python_component.py)
            py_files = list(self.python_dir.glob(f"**/*{component}*.py")) + list(
                self.python_dir.glob(f"**/python_{component}*.py")
            )

            # Check file contents for imports of the component
            if not py_files:
                for py_file in self.python_dir.glob("**/*.py"):
                    try:
                        with open(py_file, encoding="utf-8") as f:
                            content = f.read()
                            if (
                                f"python_{component}" in content.lower()
                                or f"{component}_python" in content.lower()
                            ):
                                py_files.append(py_file)
                    except Exception:
                        pass

            loc = 0
            for file_path in py_files:
                try:
                    with open(file_path, encoding="utf-8") as f:
                        loc += sum(
                            1
                            for line in f
                            if line.strip() and not line.strip().startswith("#")
                        )
                except Exception:
                    pass

            stats[component] = {
                "has_replacement": len(py_files) > 0,
                "files": len(py_files),
                "loc": loc,
            }

        return stats

    def _analyze_test_coverage(self) -> Dict[str, Dict]:
        """
        Analyze test coverage for the codebase.

        Returns:
            Dictionary with test coverage statistics
        """
        stats = {
            "unit": {"total": 0, "migrated": 0},
            "integration": {"total": 0, "migrated": 0},
            "shell": {"total": 0, "migrated": 0},
        }

        # Count pytest files
        pytest_files = list(self.test_dir.glob("**/test_*.py"))
        stats["unit"]["total"] = len(
            [f for f in pytest_files if "integration" not in str(f)]
        )
        stats["integration"]["total"] = len(
            [f for f in pytest_files if "integration" in str(f)]
        )

        # Count shell test scripts
        shell_files = list(self.project_root.glob("src/sh/*test*.sh"))
        stats["shell"]["total"] = len(shell_files)

        # Check for pytest migration
        for py_file in pytest_files:
            try:
                with open(py_file, encoding="utf-8") as f:
                    content = f.read()
                    if "pytest" in content:
                        if "integration" in str(py_file):
                            stats["integration"]["migrated"] += 1
                        else:
                            stats["unit"]["migrated"] += 1
            except Exception:
                pass

        return stats

    def _analyze_code_quality(self) -> Dict[str, any]:
        """
        Analyze code quality metrics.

        Returns:
            Dictionary with code quality metrics
        """
        metrics = {}

        # Check for Python 2 artifacts
        py2_artifacts = self._count_python2_artifacts()
        metrics["Python 2 Artifacts"] = sum(py2_artifacts.values())

        # Check type hint coverage
        type_coverage = self._check_type_coverage()
        metrics["Type Coverage"] = f"{type_coverage}%"

        # Check pre-commit configuration
        pre_commit_file = self.project_root / ".pre-commit-config.yaml"
        metrics["Pre-commit Config"] = (
            "✅ Present" if pre_commit_file.exists() else "❌ Missing"
        )

        # Check pytest configuration
        pytest_conf = self.project_root / "pytest.ini"
        metrics["Pytest Config"] = (
            "✅ Present" if pytest_conf.exists() else "❌ Missing"
        )

        return metrics

    def _count_python2_artifacts(self) -> Dict[str, int]:
        """
        Count Python 2 artifacts in the codebase.

        Returns:
            Dictionary with artifact counts
        """
        artifacts = {
            "print_statements": 0,
            "unicode_literals": 0,
            "old_super": 0,
            "xrange": 0,
            "iteritems": 0,
        }

        for py_file in self.python_dir.glob("**/*.py"):
            try:
                with open(py_file, encoding="utf-8") as f:
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
                        if ".iteritems()" in line:
                            artifacts["iteritems"] += 1
            except Exception:
                pass

        return artifacts

    def _check_type_coverage(self) -> int:
        """
        Check type hint coverage in the codebase.

        Returns:
            Percentage of functions with type hints
        """
        total_functions = 0
        typed_functions = 0

        for py_file in self.python_dir.glob("**/*.py"):
            try:
                with open(py_file, encoding="utf-8") as f:
                    content = f.read()
                    lines = content.split("\n")

                    for line in lines:
                        stripped = line.strip()
                        if stripped.startswith("def ") and ":" in stripped:
                            total_functions += 1
                            if "->" in stripped or (
                                ": " in stripped
                                and "def " not in stripped.split(": ")[1]
                            ):
                                typed_functions += 1
            except Exception:
                pass

        if total_functions == 0:
            return 0
        return int((typed_functions / total_functions) * 100)

    def _analyze_dependencies(self) -> Dict[str, Dict]:
        """
        Analyze dependencies and their utilization.

        Returns:
            Dictionary with dependency status
        """
        dependencies = {
            "pysam": {"installed": False, "utilization": None},
            "biopython": {"installed": False, "utilization": None},
            "numpy": {"installed": False, "utilization": None},
            "pandas": {"installed": False, "utilization": None},
            "pytest": {"installed": False, "utilization": None},
        }

        # Check pyproject.toml for installed packages
        pyproject_file = self.project_root / "pyproject.toml"
        if pyproject_file.exists():
            try:
                with open(pyproject_file) as f:
                    content = f.read()

                    for package in dependencies:
                        dependencies[package]["installed"] = (
                            package.lower() in content.lower()
                        )
            except Exception:
                pass

        # Check utilization in Python files
        package_imports = {pkg: 0 for pkg in dependencies}
        total_files = 0

        for py_file in self.python_dir.glob("**/*.py"):
            total_files += 1
            try:
                with open(py_file, encoding="utf-8") as f:
                    content = f.read()

                    for package in dependencies:
                        if (
                            f"import {package}" in content
                            or f"from {package}" in content
                        ):
                            package_imports[package] += 1
            except Exception:
                pass

        # Calculate utilization percentage
        if total_files > 0:
            for package, count in package_imports.items():
                dependencies[package]["utilization"] = int((count / total_files) * 100)

        return dependencies

    def _generate_detailed_sections(self, cpp_stats, python_stats, test_stats) -> str:
        """
        Generate detailed sections for the report.

        Args:
            cpp_stats: C++ component statistics
            python_stats: Python replacement statistics
            test_stats: Test coverage statistics

        Returns:
            Markdown-formatted detailed sections
        """
        report = """
## Detailed Analysis

### C++ Components by Priority

| Priority | Components | Migration Status |
|----------|------------|------------------|
"""

        # Group by priority
        priorities = {"high": [], "medium": [], "low": []}
        for component, info in self.cpp_components.items():
            priorities[info["priority"]].append(component)

        for priority, components in priorities.items():
            migrated = sum(
                1
                for c in components
                if python_stats.get(c, {}).get("has_replacement", False)
            )
            status = f"{migrated}/{len(components)} migrated"
            report += f"| {priority} | {', '.join(components)} | {status} |\n"

        # Add Python 2 artifacts section if there are any
        py2_artifacts = self._count_python2_artifacts()
        if sum(py2_artifacts.values()) > 0:
            report += """
### Python 2 Artifacts

| Artifact Type | Count |
|---------------|-------|
"""
            for artifact, count in py2_artifacts.items():
                if count > 0:
                    report += f"| {artifact} | {count} |\n"

        return report


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Generate a modernization progress report"
    )
    parser.add_argument(
        "--detailed",
        action="store_true",
        help="Include detailed sections in the report",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file for the report (default: print to stdout)",
    )

    args = parser.parse_args()

    # Get project root directory (parent of the script)
    project_root = Path(__file__).parent.parent

    # Create reporter and generate report
    reporter = ModernizationReporter(project_root)
    report = reporter.generate_report(detailed=args.detailed)

    # Output report
    if args.output:
        with open(args.output, "w") as f:
            f.write(report)
        print(f"Report written to {args.output}")
    else:
        print(report)


if __name__ == "__main__":
    main()

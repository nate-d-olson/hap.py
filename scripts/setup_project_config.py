#!/usr/bin/env python3
"""
Track progress on the pyproject.toml file for the hap.py modernization.

This script creates or updates a pyproject.toml file with modern
Python packaging configuration, dependencies, and development tools.
"""

import argparse
from pathlib import Path


def create_pyproject_toml(project_root: Path, overwrite: bool = False) -> bool:
    """
    Create a pyproject.toml file with modern Python packaging configuration.
    
    Args:
        project_root: Root directory of the project
        overwrite: Whether to overwrite an existing file
    
    Returns:
        True if the file was created or updated, False otherwise
    """
    pyproject_path = project_root / "pyproject.toml"
    
    if pyproject_path.exists() and not overwrite:
        print(f"⚠️ {pyproject_path} already exists. Use --overwrite to replace it.")
        return False
    
    pyproject_content = """[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "happy"
version = "0.4.0"
description = "Haplotype VCF comparison tools"
readme = "README.md"
requires-python = ">=3.8"
license = {text = "MIT"}
authors = [
    {name = "Illumina, Inc.", email = "support@illumina.com"}
]
maintainers = [
    {name = "Bioinformatics Tools Team", email = "support@illumina.com"}
]
keywords = ["bioinformatics", "genomics", "vcf", "variant calling", "benchmarking"]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
dependencies = [
    "pysam>=0.21.0",
    "biopython>=1.80",
    "numpy>=1.20.0",
    "pandas>=1.3.0",
    "scipy>=1.7.0",
    "pybedtools>=0.9.0",
    "pyfaidx>=0.7.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.0.0",
    "pytest-cov>=4.0.0",
    "black>=23.0.0",
    "ruff>=0.0.265",
    "isort>=5.12.0",
    "mypy>=1.0.0",
    "pre-commit>=3.0.0",
    "pyupgrade>=3.0.0",
]
docs = [
    "sphinx>=5.0.0",
    "sphinx-rtd-theme>=1.0.0",
    "myst-parser>=1.0.0",
]

[project.urls]
"Homepage" = "https://github.com/Illumina/hap.py"
"Bug Tracker" = "https://github.com/Illumina/hap.py/issues"
"Documentation" = "https://illumina.github.io/hap.py/"

[project.scripts]
hap.py = "Haplo.happyscript:main"
vcfcheck = "Haplo.python_vcfcheck:main"
blocksplit = "Haplo.python_blocksplit:main"

[tool.setuptools]
package-dir = {"" = "src/python"}
include-package-data = true

[tool.setuptools.packages.find]
where = ["src/python"]

[tool.black]
line-length = 100
target-version = ["py38"]
include = '\.pyi?$'

[tool.isort]
profile = "black"
line_length = 100

[tool.ruff]
line-length = 100
target-version = "py38"
select = ["E", "F", "W", "I", "N", "UP", "B", "A"]
ignore = ["E203", "E501"]

[tool.ruff.per-file-ignores]
"__init__.py" = ["F401"]

[tool.mypy]
python_version = "3.8"
warn_return_any = true
warn_unused_configs = true
disallow_untyped_defs = false
disallow_incomplete_defs = false

[[tool.mypy.overrides]]
module = ["pysam.*", "Bio.*", "numpy.*", "pandas.*", "pybedtools.*", "pyfaidx.*"]
ignore_missing_imports = true

[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = "test_*.py"
python_functions = "test_*"
"""
    
    with open(pyproject_path, "w") as f:
        f.write(pyproject_content)
    
    print(f"✅ Created {pyproject_path}")
    return True


def create_pre_commit_config(project_root: Path, overwrite: bool = False) -> bool:
    """
    Create a pre-commit config file with hooks for code quality.
    
    Args:
        project_root: Root directory of the project
        overwrite: Whether to overwrite an existing file
    
    Returns:
        True if the file was created or updated, False otherwise
    """
    pre_commit_path = project_root / ".pre-commit-config.yaml"
    
    if pre_commit_path.exists() and not overwrite:
        print(f"⚠️ {pre_commit_path} already exists. Use --overwrite to replace it.")
        return False
    
    pre_commit_content = """repos:
-   repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v4.4.0
    hooks:
    -   id: trailing-whitespace
    -   id: end-of-file-fixer
    -   id: check-yaml
    -   id: check-added-large-files
    -   id: check-ast
    -   id: check-json
    -   id: check-merge-conflict
    -   id: debug-statements
    -   id: detect-private-key

-   repo: https://github.com/psf/black
    rev: 23.3.0
    hooks:
    -   id: black

-   repo: https://github.com/pycqa/isort
    rev: 5.12.0
    hooks:
    -   id: isort

-   repo: https://github.com/charliermarsh/ruff-pre-commit
    rev: v0.0.265
    hooks:
    -   id: ruff
        args: [--fix, --exit-non-zero-on-fix]

-   repo: https://github.com/asottile/pyupgrade
    rev: v3.4.0
    hooks:
    -   id: pyupgrade
        args: [--py38-plus]

-   repo: https://github.com/pre-commit/mirrors-mypy
    rev: v1.3.0
    hooks:
    -   id: mypy
        additional_dependencies: [types-requests]
        args: [--ignore-missing-imports]
"""
    
    with open(pre_commit_path, "w") as f:
        f.write(pre_commit_content)
    
    print(f"✅ Created {pre_commit_path}")
    return True


def main():
    """Main function to set up project configuration."""
    parser = argparse.ArgumentParser(
        description="Set up modern Python project configuration for hap.py"
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Overwrite existing files"
    )
    
    args = parser.parse_args()
    
    # Get project root directory (parent of the script)
    project_root = Path(__file__).parent.parent
    
    # Create pyproject.toml
    created_pyproject = create_pyproject_toml(project_root, args.overwrite)
    
    # Create pre-commit config
    created_pre_commit = create_pre_commit_config(project_root, args.overwrite)
    
    # Final message
    if created_pyproject or created_pre_commit:
        print("\n✅ Project configuration set up successfully!")
        print("\nNext steps:")
        print("1. Install dependencies: pip install -e .[dev]")
        print("2. Install pre-commit: pre-commit install")
        print("3. Run pre-commit on all files: pre-commit run --all-files")
    else:
        print("\n⚠️ No changes made. Use --overwrite to update existing files.")


if __name__ == "__main__":
    main()

# Pre-commit Hooks Guide for hap.py Migration

## Overview

This repository uses pre-commit hooks to automate code quality checks and assist with Python 3 migration. These tools help ensure consistent style, detect common issues, and automate many Python 2 to 3 conversions.

## Installation

```bash
# Install pre-commit
pip install pre-commit

# Install hooks
pre-commit install
```

## Configured Hooks

The following hooks are configured in .pre-commit-config.yaml:

1. **Basic checks**: Trailing whitespace, file endings, merge conflicts
2. **Black**: Code formatting
3. **Ruff**: Fast Python linter with auto-fixes
4. **isort**: Import sorting
5. **pyupgrade**: Automated Python syntax upgrading
6. **mypy**: Static type checking (optional)

## Hooks Usage During Migration

### Basic Workflow

```bash
# Run all hooks on changed files (automatic during commit)
git add file1.py file2.py
git commit -m "Message"  # Hooks run automatically

# Run all hooks manually on specific files
pre-commit run --files src/python/path/to/file.py

# Run all hooks on all files
pre-commit run --all-files
```

### Specific Hooks for Migration Tasks

```bash
# Format code with Black
pre-commit run black --files src/python/path/to/file.py

# Fix Python 2 to 3 syntax with pyupgrade
pre-commit run pyupgrade --files src/python/path/to/file.py

# Run linting with auto-fixes
pre-commit run ruff --files src/python/path/to/file.py

# Sort imports
pre-commit run isort --files src/python/path/to/file.py
```

### Recommended Migration Process

For each Python file being migrated:

1. Run initial 2to3 tool conversion
2. Run `pre-commit run pyupgrade --files file.py` to modernize Python 3 syntax
3. Run `pre-commit run ruff --files file.py` to detect and fix common issues
4. Run `pre-commit run black --files file.py` for consistent formatting
5. Manually address any remaining issues
6. Run `pre-commit run --all` as a final check

### Troubleshooting

If a hook fails:
1. Read the error message to understand the issue
2. Fix the issue manually if needed
3. Re-run the hook to verify the fix

If a hook causes undesired changes:
1. Check the hook configuration in .pre-commit-config.yaml
2. Adjust hook parameters or exclude patterns as needed

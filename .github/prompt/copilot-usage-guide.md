# GitHub Copilot Configuration

This file provides guidance on how to effectively use GitHub Copilot with the hap.py codebase.

## Effective Prompting

When using GitHub Copilot for this project:

1. **Be specific about bioinformatics context**: Mention that you're working with genomic data, VCF files, or variant calling when relevant
2. **Reference language version**: Specify Python 3 when requesting code suggestions
3. **Include type hints**: Provide type information in your requests
4. **Describe memory constraints**: Mention if you need code optimized for large genomic datasets
5. **Reference existing patterns**: Ask for code that follows patterns found elsewhere in the codebase

## Useful Prefixes for Inline Comments

- `# TODO:` - Mark areas that need future work
- `# FIXME:` - Identify bugs or issues to be fixed
- `# OPTIMIZE:` - Areas that need performance improvement
- `# HACK:` - Temporary solutions that should be revisited
- `# PYTHON3:` - Marks code specifically adjusted for Python 3 compatibility

## File-Level Hints

You can add file-level hints with special comments at the top of files:

```python
# language=python3
# domain=bioinformatics
# context=variant_calling
```

## Function-Level Prompting

When asking Copilot to generate or modify functions, provide:

1. The function signature with type hints
2. A clear docstring explaining the purpose
3. Any constraints or requirements
4. Example input/output if applicable

## Testing-Related Hints

When working on tests, use comments like:

```python
# Generate test cases for the variant normalization function
# Should test: SNPs, indels, complex variants, edge cases
```

## Debugging Assistance

When debugging, describe the error specifically:

```python
# Debug this IndexError: list index out of range when parsing VCF records
```

## Code Review Assistance

Use Copilot to help review code with comments like:

```python
# Review this function for Python 3 compatibility issues
# Check this algorithm for memory efficiency with large datasets
```

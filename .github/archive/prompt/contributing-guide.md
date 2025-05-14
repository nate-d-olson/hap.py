# Contributing to hap.py with GitHub Copilot

This guide provides instructions for contributing to the hap.py project using GitHub Copilot to enhance your development workflow.

## Getting Started as a Contributor

### Setting Up Your Environment

1. Fork the repository on GitHub
2. Clone your fork locally
3. Set up your development environment following the instructions in `.github/prompt/environment-setup.md`
4. Install the GitHub Copilot and GitHub Copilot Chat extensions in VS Code

### Understanding the Codebase

Use GitHub Copilot to help navigate and understand the codebase:

```
@copilot Give me an overview of how variant comparison works in the hap.py codebase
```

```
@copilot Explain the relationship between the Python and C++ components in this project
```

## Development Workflow

### 1. Choose an Issue to Work On

Browse existing issues or identify new areas for improvement. Use Copilot to help understand the issue:

```
@copilot Help me understand what's causing this Python 3 compatibility issue in the variant parser
```

### 2. Create a Feature Branch

```bash
git checkout -b feature/your-feature-name
```

### 3. Make Incremental Changes

Use GitHub Copilot to help with implementation:

**Step 1: Understand the problem**
```
@copilot Help me understand how haplotype matching is implemented in this file
```

**Step 2: Design a solution**
```
@copilot How should I approach updating this code to handle phased variants correctly in Python 3?
```

**Step 3: Implement the changes**
```
@copilot Help me update this function to use Python 3 string handling
```

### 4. Test Your Changes

Run the tests after each significant change:

```bash
cd /tmp/happy-build
src/sh/run_tests.sh
```

If tests fail, use Copilot to help debug:

```
@copilot This test is failing with this error message. How do I fix it?
```

### 5. Document Your Changes

Use Copilot to help write clear documentation:

```
@copilot Help me write Google-style docstrings for this function I've updated
```

```
@copilot How should I document this change in the RELEASES.md file?
```

### 6. Submit Your Pull Request

Create a comprehensive pull request description:

```
@copilot Help me write a pull request description for these Python 3 compatibility changes
```

## Contribution Guidelines

### Code Style

- Follow the style guidelines in `.github/instructions/`
- Use Copilot to help with formatting:

```
@copilot Reformat this code according to PEP 8 guidelines
```

### Adding Type Hints

Use Copilot to help add type hints during modernization:

```
@copilot Add appropriate type hints to this function following our project standards
```

### Testing Contributions

Always add or update tests for your changes:

```
@copilot Help me write a test case for this new variant normalization function
```

### Documentation Contributions

Update relevant documentation for your changes:

```
@copilot How should I update the documentation to reflect these changes to the VCF parser?
```

## Bioinformatics Best Practices

When contributing code that deals with genomic data:

- Ensure proper handling of chromosome names (str vs int, X/Y handling)
- Validate coordinates and reference bases
- Handle edge cases like complex variants and structural variants
- Consider memory efficiency for large genomic datasets

Use Copilot to implement these practices:

```
@copilot How should I validate chromosome coordinates in this function?
```

```
@copilot How can I make this variant comparison more memory efficient for whole genome data?
```

## Working with the Community

- Be respectful and constructive in pull request discussions
- Clearly explain your reasoning for implementation choices
- Be open to feedback and willing to make changes
- Credit others for their ideas and contributions

GitHub Copilot can help you communicate clearly:

```
@copilot Help me write a clear response to this code review comment about my implementation approach
```

## Further Resources

- [Python 3 Migration Guide](/.github/prompt/python-migration-help.md)
- [Bioinformatics Coding Standards](/.github/instructions/bioinformatics.instructions.md)
- [Troubleshooting Guide](/.github/prompt/troubleshooting.md)

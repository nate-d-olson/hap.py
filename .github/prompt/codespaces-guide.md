# GitHub Codespaces for hap.py Development

This guide provides instructions specifically for using GitHub Codespaces with GitHub Copilot for development on the hap.py project.

## Setting Up Codespaces

### Creating a New Codespace

1. Navigate to the GitHub repository (github.com/nate-d-olson/hap.py)
2. Click on the "Code" button, then the "Codespaces" tab
3. Click "Create codespace on main"
4. Wait for the environment to build

### Configuring the Codespace

When your Codespace loads, GitHub Copilot and other essential extensions will be pre-installed. The workspace settings will be automatically applied from the repository configuration.

## Building hap.py in Codespaces

Run the following commands in the terminal:

```bash
# Create a build directory
TEMP_BUILD="/tmp/happy-build"
python install.py $TEMP_BUILD

# Run tests to verify installation
cd $TEMP_BUILD
src/sh/run_tests.sh
```

## Using GitHub Copilot in Codespaces

### Chat Panel

1. Open the GitHub Copilot chat panel using the chat icon in the sidebar or using the keyboard shortcut `Ctrl+I` (or `Cmd+I` on macOS)
2. Start your prompt with "@copilot" to interact with Copilot

### Workspace-Specific Prompts

```
@copilot Help me understand the variant comparison logic in this bioinformatics tool
```

```
@copilot How can I optimize this VCF parsing function for large genomic files?
```

```
@copilot Fix Python 3 compatibility issues in this code
```

### Using #file Context

You can refer to specific files in your Copilot queries:

```
@copilot #file:src/python/haplotypes.py Explain how this module handles haplotype comparison
```

## Testing Changes

After making changes, rebuild and test your changes:

```bash
# Rebuild the project
cd /workspaces/hap.py
python install.py $TEMP_BUILD

# Run tests
cd $TEMP_BUILD
src/sh/run_tests.sh
```

## Persisting Your Work

Codespaces automatically commit your changes to a branch. To formalize your changes:

1. Create a new branch for your work
2. Commit your changes with a descriptive message
3. Push your branch to the repository
4. Create a pull request

## Terminal Commands for Common Tasks

### Python Linting

```bash
cd /workspaces/hap.py
black src/python
```

### Running Specific Tests

```bash
cd $TEMP_BUILD
src/sh/run_test.sh [specific_test_name]
```

### Checking Python Version Compatibility

```bash
cd /workspaces/hap.py
2to3 -p src/python/[specific_file.py]
```

## Sharing Your Codespace

If you need to share your work-in-progress with a collaborator:

1. Click on the Codespaces menu in the lower left corner
2. Select "Share Current Codespace"
3. Copy and share the provided link

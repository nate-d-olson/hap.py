# Environment Setup for hap.py Development

This guide will help you set up your development environment for working on the hap.py project with GitHub Copilot integration.

## Prerequisites

- Python 3.6+
- CMake 3.12+
- C++ compiler (gcc 7+ recommended)
- Git
- Visual Studio Code or GitHub Codespaces

## Local Development Setup

### 1. Clone the Repository

```bash
git clone https://github.com/nate-d-olson/hap.py.git
cd hap.py
```

### 2. Set Up Python Environment

Create a virtual environment:

```bash
python -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate
```

### 3. Install Dependencies

```bash
# Install Python dependencies
pip install -r happy.requirements.txt
```

### 4. Build the Project

```bash
# Create a temporary build directory
TEMP_BUILD="/tmp/happy-build"
python install.py $TEMP_BUILD
```

### 5. Run Tests

```bash
cd $TEMP_BUILD
src/sh/run_tests.sh
```

## GitHub Codespaces Setup

1. Open the GitHub repository in your browser
2. Click on "Code" button and select "Open with Codespaces"
3. Create a new codespace
4. Once the codespace is ready, open a terminal and follow steps 2-5 from the Local Development Setup

## Visual Studio Code Configuration

1. Install recommended extensions:
   - GitHub Copilot
   - GitHub Copilot Chat
   - Python
   - C/C++
   - CMake

2. Configure VS Code settings:
   - Open Command Palette (Ctrl+Shift+P / Cmd+Shift+P)
   - Search for "Preferences: Open Settings (JSON)"
   - Add the settings from `.github/instructions/vscode.instructions.md`

## Using GitHub Copilot for Development

### Enabling Copilot Features

1. Sign in to GitHub in VS Code
2. Ensure GitHub Copilot is enabled
3. Use Copilot Chat (Ctrl+I / Cmd+I) for conversational assistance
4. Accept or modify inline suggestions (Tab key)

### Useful Copilot Commands

- `/explain` - Explain the currently selected code
- `/tests` - Generate tests for the current function
- `/fix` - Suggest a fix for the current code
- `/optimize` - Suggest performance optimizations
- `/doc` - Generate documentation

## Common Development Tasks

### Converting Python 2 to Python 3

Use Copilot to help with conversion tasks:

```
@copilot Convert this Python 2 code to Python 3 following our project standards
```

### Adding Type Hints

```
@copilot Add type hints to this function following our project standards
```

### Optimizing Code for Performance

```
@copilot Help optimize this code for large genomic datasets
```

### Debugging Test Failures

```
@copilot Debug this test failure in the hap.py codebase
```

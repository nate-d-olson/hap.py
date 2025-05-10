---
applyTo: "**"
---
# Visual Studio Code Settings for hap.py Development

## Recommended Extensions

- **Python**: Microsoft's Python extension for VS Code
- **GitHub Copilot**: AI-assisted code suggestions
- **GitHub Copilot Chat**: Interactive coding assistance
- **C/C++**: Microsoft's C++ extension for VS Code
- **CMake**: CMake language support
- **Python Docstring Generator**: Generate Google-style docstrings
- **Pylance**: Python language server for improved type checking

## Workspace Settings

### Python Settings
```json
{
  "python.linting.enabled": true,
  "python.linting.flake8Enabled": true,
  "python.linting.pylintEnabled": false,
  "python.formatting.provider": "black",
  "python.formatting.blackArgs": ["--line-length", "88"],
  "editor.formatOnSave": true,
  "python.analysis.typeCheckingMode": "basic"
}
```

### C++ Settings
```json
{
  "C_Cpp.default.cppStandard": "c++11",
  "C_Cpp.default.includePath": [
    "${workspaceFolder}/src/c++/**",
    "${workspaceFolder}/external/**"
  ]
}
```

## Using GitHub Copilot

Enable inline suggestions and GitHub Copilot chat to get context-aware help while developing. 

To get the most out of Copilot for bioinformatics work:

1. Use comments to explain genomic concepts when writing code
2. Provide file format specs when working with specialized formats
3. Describe performance considerations for large genomic datasets
4. Reference specific biological concepts that might not be common knowledge

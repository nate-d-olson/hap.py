# Vibe Coding with GitHub Copilot for hap.py

This guide provides instructions for effectively using GitHub Copilot with Vibe coding for the hap.py project, either in GitHub Codespaces or Visual Studio Code locally.

## Setting Up Your Environment

### Visual Studio Code Setup
1. Install the [GitHub Copilot](https://marketplace.visualstudio.com/items?itemName=GitHub.copilot) extension
2. Install the [GitHub Copilot Chat](https://marketplace.visualstudio.com/items?itemName=GitHub.copilot-chat) extension
3. Optional: Install [Python](https://marketplace.visualstudio.com/items?itemName=ms-python.python) and [C/C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools) extensions

### GitHub Codespaces Setup
1. No additional setup required - Copilot extensions are pre-installed
2. Access Codespaces through GitHub repository page

## Using Vibe Coding

Vibe coding enables a more natural, conversational approach to programming with GitHub Copilot. Use these techniques for an optimal experience:

### Conversation Starters

Begin your conversations with GitHub Copilot Chat using these prompts:

- **Code Explanation**: "Explain how the variant comparison logic works in `src/python/haplotype_comparison.py`"
- **Code Modernization**: "Help me update this Python 2 code to Python 3 standards"
- **Bug Fixing**: "Debug why this variant normalization function fails with complex indels"
- **Feature Development**: "Design a function to calculate stratified precision/recall metrics"
- **Code Optimization**: "Optimize this function for handling large VCF files"

### File Context Commands

Use these commands to leverage repository context:

- **#file**: Reference specific files in your queries
  - Example: "Show me how to improve error handling in #file:src/python/vcf_preprocessing.py"

- **#codebase**: Ask about the broader codebase structure
  - Example: "#codebase What modules handle variant normalization?"

- **#workspace**: Get help with workspace-wide questions
  - Example: "#workspace How can I improve test coverage across the project?"

### Bioinformatics-Specific Commands

For genomics and bioinformatics-specific assistance:

- **#bioinformatics**: "Help me implement proper variant normalization for complex indels"
- **#vcf**: "Explain the best way to handle phase information in this variant comparison"
- **#performance**: "Optimize this code for processing large chromosome variant data"

## Best Practices for hap.py Development

1. **Incremental Migration**: Work in small, testable increments when migrating code
2. **Test-First Approach**: Run tests frequently to catch migration issues early
3. **Documentation**: Include detailed docstrings when modernizing code
4. **Context Provision**: Give Copilot context about bioinformatics concepts when needed
5. **Pattern Replication**: Ask Copilot to follow established patterns in the codebase

## Example Workflows

### Python 2 to 3 Migration Example
```
@copilot Help me migrate this Python 2 code to Python 3:

def process_variants(vcf_file):
    for line in open(vcf_file, "r"):
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 8:
            continue
        yield dict(zip(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"], parts[:8]))
```

### Debugging Test Failures
```
@copilot I'm getting this error in the test suite: 
"TypeError: 'dict_keys' object is not subscriptable" 
in src/python/haplotype_matching.py line 234. 
Help me fix this Python 3 compatibility issue.
```

### Performance Optimization
```
@copilot This VCF parsing function is taking too long with large files. 
Can you help optimize it for better performance with large genomic datasets?
```

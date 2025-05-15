# Add regex-based check for Python 2 style exceptions
def _check_exception_syntax_regex(self, content, filepath):
    """Check for old Python 2 style exception handling"""
    # Pattern to match Python 2 exception handling (except X, e:)
    pattern = r'except\s+([A-Za-z0-9_]+)\s*,\s*([A-Za-z0-9_]+)\s*:'
    
    for match in re.finditer(pattern, content):
        lineno = content[:match.start()].count('\n') + 1
        self.issues["except_syntax"].append({
            "file": filepath,
            "line": lineno,
            "message": "Old exception syntax (except Exception, e)"
        })

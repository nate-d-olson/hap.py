# Python 2 to 3 Migration Guidance

When working on Python 2 to 3 migration issues in the hap.py codebase, consider the following common patterns and their solutions:

## Current Priority Areas

The following modules require immediate attention (as of May 14, 2025):

1. **Haplo Module**: Focus on string/bytes handling in variant processing
2. **Tools Module**: Address file I/O and encoding issues
3. **Somatic Module**: Fix exception syntax and sequence processing

## Common Migration Issues

### Print Statements
- **Python 2**: `print "Hello"`
- **Python 3**: `print("Hello")`

### Division
- **Python 2**: `5/2` returns `2` (integer division)
- **Python 3**: `5/2` returns `2.5` (float division)
- **Solution**: Use `//` for integer division consistently

### Exception Handling
- **Python 2**: `except ValueError, e:`
- **Python 3**: `except ValueError as e:`

### Dictionaries
- **Python 2**: `.iteritems()`, `.iterkeys()`, `.itervalues()`
- **Python 3**: `.items()`, `.keys()`, `.values()` (which return views, not lists)

### Ranges
- **Python 2**: `xrange()` for iteration, `range()` creates lists
- **Python 3**: `range()` for iteration (returns generator-like object)

### Unicode vs Bytes
- **Python 2**: Strings are bytes, unicode is separate
- **Python 3**: Strings are unicode, bytes are separate
- **Solution**: Add explicit encoding/decoding where needed

### File I/O
- **Python 2**: File operations default to text mode
- **Python 3**: Binary mode needs explicit 'b' flag, text mode handles encodings
- **Solution**: Always specify mode explicitly: `open(file, 'rt')` or `open(file, 'rb')`

### Absolute Imports
- **Python 2**: Implicit relative imports
- **Python 3**: Explicit relative imports required
- **Solution**: Use `from . import module` instead of `import module`

### Removed Functions
- `cmp()`, `reduce()` (moved to functools), `execfile()`
- **Solution**: Replace with appropriate alternatives

## Helpful Libraries

- `six`: Compatibility library for Python 2 and 3
- `future`: Forward compatibility from Python 2
- `2to3`: Conversion tool (already applied, but useful for specific issues)

## Testing Approach

When fixing Python 3 compatibility issues:

1. Identify the root cause through error messages
2. Make minimal changes needed for compatibility
3. Test both functionality and performance
4. Watch for subtle differences in behavior (especially string/bytes handling)

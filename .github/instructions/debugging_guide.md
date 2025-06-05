# Debugging Guide for hap.py Integration Tests

## Common Issues and Solutions

### RTG Tool Detection Issues

**Symptoms:**
- "WARNING:root:Executable for rtg not found" messages
- "rtg: command not found" errors in tests

**Solutions:**
1. Ensure the RTG path is correctly set:
   ```bash
   # Check if the RTG executable exists
   ls external/rtg-tools-3.12.1/rtg

   # Create symbolic link if needed
   ln -s external/rtg-tools-3.12.1/rtg .

   # Add to PATH for tests
   export PATH=$PATH:$(pwd)/
   ```

2. Update tests to explicitly specify the RTG path:
   ```python
   # Use the rtg_executable fixture in tests
   def test_something(rtg_executable):
       # Use the fixture in test
       cmd = ["hap.py", "--engine=vcfeval", f"--engine-vcfeval-path={rtg_executable}"]
   ```

3. Check for RTG detection in `__init__.py`:
   ```python
   def init():
       # Also check our custom RTG location
       if os.path.exists("external/rtg-tools-3.12.1/rtg"):
           return
       if shutil.which("rtg"):
           return
       logging.warning("Executable for rtg not found")
   ```

### SDF Template Directory Issues

**Symptoms:**
- "directory already exists" errors from RTG format command
- Hanging tests due to filesystem issues

**Solutions:**
1. Use `tempfile.mkdtemp()` instead of `tempfile.NamedTemporaryFile()` for directories
2. Properly clean up temporary directories before creating new ones
3. Check if the test is crossing filesystem boundaries with temporary directories

### VCF Header Validation Errors

**Symptoms:**
- Missing FILTER field errors
- Duplicate FORMAT entries

**Solutions:**
1. Check VCF header validation in `vcf.py`:
   ```python
   def _check_header(header):
       # Ensure FILTER column is properly detected
       if "FILTER" not in header:
           # Look for it case-insensitively
           for h in header:
               if h.upper() == "FILTER":
                   header[header.index(h)] = "FILTER"
                   break
   ```

2. Deduplicate FORMAT entries before validation

### Missing Output Files

**Symptoms:**
- `roc.tsv` output file is missing
- Test fails because expected output isn't found

**Solutions:**
1. Verify that the output directory exists and is writable
2. Check that the command correctly specifies the output location
3. Ensure any preprocessing steps that create directories are working

## Debugging Process for Integration Tests

1. **Review Test Output:**
   ```bash
   pytest tests/integration/ -v | tee integration_test_output.txt
   ```

2. **Locate the First Failing Test:**
   Identify the first test that fails and focus on that.

3. **Examine Temporary Files:**
   Check temporary directories for partial outputs that might provide clues.

4. **Add Debug Logging:**
   ```python
   import logging
   logging.basicConfig(level=logging.DEBUG)
   ```

5. **Run Specific Test with Full Debug Info:**
   ```bash
   pytest tests/integration/specific_test.py::test_specific_function -v -s
   ```

6. **Manually Run Commands:**
   Extract the failing command and run it manually to observe direct output.

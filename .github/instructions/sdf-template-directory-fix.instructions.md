````markdown
# SDF Template Directory Fix Guide

This document provides a detailed explanation of the SDF template directory issue that causes many integration test failures, and how to fix it.

## Problem Description

When running RTG tools with the `format` command, it fails with the following error:

```
Error: The directory "/var/folders/.../vcfeval.sdf.xyz" already exists. Please remove it first or choose a different directory.
```

This happens because the current implementation in `vcfeval.py` uses `tempfile.NamedTemporaryFile` to generate a name for the SDF template directory, but this actually creates a file, not a directory. RTG's `format` command then fails because it expects to create the directory itself but finds a file already there.

## Root Cause Analysis

The problematic code in `vcfeval.py` looks something like this:

```python
try:
    with tempfile.NamedTemporaryFile(
        dir=args.scratch_prefix, prefix="vcfeval.sdf", suffix=".dir"
    ) as stf:
        template_dir = stf.name

    # Remove template dir if it already exists (RTG format will fail otherwise)
    if os.path.exists(template_dir):
        logging.warning(f"SDF template directory {template_dir} already exists. Removing it before running rtg format.")
        shutil.rmtree(template_dir)
    os.makedirs(template_dir, exist_ok=True)
```

The issue is:
1. `tempfile.NamedTemporaryFile` creates an actual file
2. When the `with` block ends, the file is deleted but not its name
3. The code then checks if directory exists (it doesn't, because a file was there)
4. Then it tries to create that directory
5. But RTG's `format` command tries to create the directory itself and fails

## Correct Implementation

The solution is to use `tempfile.mkdtemp()` which properly creates a unique temporary directory:

```python
try:
    # Use mkdtemp to create a unique directory for the SDF template
    template_dir = tempfile.mkdtemp(dir=args.scratch_prefix, prefix="vcfeval.sdf.")

    # No need to check existence or explicitly create the directory,
    # since mkdtemp guarantees a newly created, empty directory
    args.engine_vcfeval_template = template_dir
```

This approach is better because:
1. `mkdtemp` creates a directory, not a file
2. The directory name is guaranteed to be unique
3. There's no need for additional existence checks or creation steps
4. RTG's `format` command will work correctly with this directory

## Additional Considerations

When fixing this issue, also ensure:

1. **Cleanup**: The temporary directory should be properly cleaned up after use
2. **Error Handling**: Use try/finally or context managers to ensure cleanup
3. **Logging**: Add appropriate log messages to aid debugging
4. **Multiple Tests**: Ensure each test gets a unique directory

## Example Complete Fix

```python
try:
    # Create a unique temporary directory for the SDF template
    template_dir = tempfile.mkdtemp(dir=args.scratch_prefix, prefix="vcfeval.sdf.")
    args.engine_vcfeval_template = template_dir

    # Quote paths for shell safety
    quoted_engine = shlex.quote(args.engine_vcfeval)
    quoted_template = shlex.quote(args.engine_vcfeval_template)
    quoted_ref = shlex.quote(args.ref)

    runme = f"{quoted_engine} format -o {quoted_template} {quoted_ref}"

    logging.info(f"Creating SDF template with command: {runme}")
    process = subprocess.Popen(
        runme,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout, stderr = process.communicate()
    rc = process.returncode

    if rc != 0:
        error_msg = (
            f"Error running rtg tools. Return code was {rc}, "
            f"output: {stdout} / {stderr}"
        )
        logging.error(error_msg)
        raise subprocess.SubprocessError(error_msg)
except Exception as e:
    logging.error("Failed to create SDF template: %s", str(e))
    if template_dir and os.path.exists(template_dir):
        with contextlib.suppress(OSError):
            shutil.rmtree(template_dir)
    raise
```

## Testing the Fix

After implementing this fix, you can verify it works correctly by running:

```bash
pytest tests/integration/test_integration.py::test_vcfeval_integration -v
```

This specific test should pass without the "directory already exists" error.
````

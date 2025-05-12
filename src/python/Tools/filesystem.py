"""Core filesystem utilities."""

import os
import errno
from typing import Optional


def ensure_directory(path: str) -> None:
    """Create a directory and its parents if they don't exist.

    Args:
        path: The directory path to create

    Raises:
        OSError: If directory creation fails
        RuntimeError: If path exists but is not a directory
    """
    try:
        os.makedirs(os.path.abspath(path))
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if not os.path.isdir(path):
                raise RuntimeError(f"Failed to create directory: {path}")
        else:
            raise


def which(program: str) -> Optional[str]:
    """Find an executable in PATH.

    Args:
        program: Name of the program to find

    Returns:
        Full path to the program if found, None otherwise
    """

    def is_executable(fpath: str) -> bool:
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    # If path is provided, check only that
    fpath, _ = os.path.split(program)
    if fpath:
        return program if is_executable(program) else None

    # Check PATH
    for path in os.environ["PATH"].split(os.pathsep):
        path = path.strip('"')
        exe_path = os.path.join(path, program)
        if is_executable(exe_path):
            return exe_path

    return None


# Alias for backward compatibility
mkdir_p = ensure_directory

"""
Utility functions for test path management.
"""

import os
import shutil
import sys
from pathlib import Path


def get_project_root() -> Path:
    """Get the project root directory."""
    # Try to find the project root using multiple methods
    try:
        # First try using __file__ if available
        return Path(__file__).resolve().parent.parent
    except NameError:
        # If __file__ is not available (e.g., running from interactive shell)
        pass

    # Try to find root by looking for known files
    current = Path.cwd()
    while str(current) != "/":
        if (current / "setup.py").exists():
            return current
        current = current.parent

    raise RuntimeError("Could not find project root directory")


def setup_python_path():
    """Set up Python paths for importing modules."""
    project_root = get_project_root()
    src_path = project_root / "src"
    if not src_path.exists():
        raise RuntimeError(f"Source directory not found at {src_path}")

    if str(src_path) not in sys.path:
        sys.path.insert(0, str(src_path))
    os.environ["PYTHONPATH"] = str(src_path)

    # Add external dependencies to PYTHONPATH if they exist
    external_path = project_root / "external"
    if external_path.exists():
        os.environ["PYTHONPATH"] = f"{str(external_path)}:{os.environ['PYTHONPATH']}"


def get_example_data_dir() -> Path:
    """Get the example data directory path."""
    path = get_project_root() / "example"
    if not path.exists():
        raise RuntimeError(f"Example data directory not found at {path}")
    return path


def get_build_dir() -> Path:
    """Get the build directory path."""
    path = get_project_root() / "build"
    if not path.exists():
        raise RuntimeError(f"Build directory not found at {path}")
    return path


def get_external_dir() -> Path:
    """Get the external dependencies directory path."""
    path = get_project_root() / "external"
    if not path.exists():
        # Create external directory if it doesn't exist
        path.mkdir(parents=True, exist_ok=True)
        print(f"Created external dependencies directory at {path}")
    return path


def get_libexec_dir() -> Path:
    """Get the libexec directory path."""
    path = get_project_root() / "libexec"
    if not path.exists():
        # Create libexec directory if it doesn't exist
        path.mkdir(parents=True, exist_ok=True)
        print(f"Created libexec directory at {path}")
    return path


# Check for required external dependencies
def check_external_dependencies():
    """Check if required external dependencies are available."""
    project_root = get_project_root()

    # Check for other required tools
    required_tools = ["bcftools", "samtools", "pysam"]
    missing_tools = []
    for tool in required_tools:
        if not shutil.which(tool):
            missing_tools.append(tool)

    if missing_tools:
        print(f"Warning: Missing required tools: {', '.join(missing_tools)}")
        print("Some tests may fail due to missing dependencies.")


# Run dependency check when module is imported
check_external_dependencies()

# Set up paths when module is imported
setup_python_path()

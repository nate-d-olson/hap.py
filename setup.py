# Modern setup using scikit-build
import toml
from setuptools import find_packages, setup

# Load version from pyproject.toml
with open("pyproject.toml", encoding="utf-8") as f:
    pyproject_data = toml.load(f)
    version = pyproject_data["project"]["version"]

setup(
    name="hap.py",
    version=version,
    packages=find_packages(where="src/python"),
    package_dir={"": "src/python"},
    # Cython extensions are now built by CMake via scikit-build
    # No ext_modules needed here directly
    # Ensure CMakeLists.txt is configured to build and install Cython modules
    # into the correct location within the site-packages directory.
    # scikit-build will drive the CMake build process.
    # Other metadata (author, description, etc.) is in pyproject.toml
)

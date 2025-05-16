#!/usr/bin/env python3
"""
Script to verify the Python 3 build system compatibility.

This script checks:
1. CMake configuration files
2. Installation scripts
3. Required dependencies

Usage:
    python3 verify_py3_build.py [--fix] [--verbose]
"""

import argparse
import importlib
import importlib.util
import logging
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Build system files to check
CMAKE_FILES = [
    "CMakeLists.txt",
    "CMakeLists.txt.py3",
    "src/cmake/CythonSupport.cmake",
    "src/cmake/CythonSupport.py3.cmake",
]

INSTALL_SCRIPTS = [
    "install.py",
    "install_py3.py",
]

REQUIRED_DEPENDENCIES = [
    "numpy",
    "cython",
    "pysam",
    "pandas",
    "scipy",
]


def check_cmake_files(fix=False):
    """Check CMake configuration files for Python 3 compatibility issues"""
    issues_found = 0
    
    for cmake_file in CMAKE_FILES:
        if not os.path.exists(cmake_file):
            logger.warning(f"File not found: {cmake_file}")
            continue
        
        with open(cmake_file, encoding="utf-8") as f:
            content = f.read()
        
        # Check for potential Python 2 specific references
        py2_patterns = [
            r"PYTHON_VERSION 2",
            r"Python2",
            r"Python 2",
            r"python2",
            r"/usr/bin/python(?![3])"
        ]
        
        issues = []
        for pattern in py2_patterns:
            for match in re.finditer(pattern, content):
                line_num = content[:match.start()].count('\n') + 1
                issues.append((line_num, match.group(0)))
        
        if issues:
            issues_found += 1
            logger.info(f"Found potential Python 2 references in {cmake_file}:")
            for line_num, text in sorted(issues):
                logger.info(f"  Line {line_num}: {text}")
        
        # Check for FindPythonInterp/Libs (deprecated in new CMake) vs FindPython3
        if "find_package(PythonInterp" in content and "find_package(Python3" not in content:
            issues_found += 1
            logger.info(f"{cmake_file}: Uses deprecated FindPythonInterp instead of FindPython3")
            
            if fix:
                # Replace deprecated CMake Python find commands
                new_content = content.replace(
                    "find_package(PythonInterp", "find_package(Python3 COMPONENTS Interpreter"
                )
                new_content = new_content.replace(
                    "find_package(PythonLibs", "find_package(Python3 COMPONENTS Development"
                )
                
                # Update variables
                new_content = new_content.replace("PYTHON_EXECUTABLE", "Python3_EXECUTABLE")
                new_content = new_content.replace("PYTHON_LIBRARIES", "Python3_LIBRARIES")
                new_content = new_content.replace("PYTHON_INCLUDE_DIRS", "Python3_INCLUDE_DIRS")
                
                with open(cmake_file, "w", encoding="utf-8") as f:
                    f.write(new_content)
                logger.info(f"  Fixed CMake Python3 references in {cmake_file}")
    
    return issues_found


def check_install_scripts(fix=False):
    """Check installation scripts for Python 3 compatibility issues"""
    issues_found = 0
    
    for script in INSTALL_SCRIPTS:
        if not os.path.exists(script):
            logger.warning(f"File not found: {script}")
            continue
        
        with open(script, encoding="utf-8") as f:
            content = f.read()
        
        # Check shebang line
        if not content.startswith("#!/usr/bin/env python3"):
            issues_found += 1
            logger.info(f"{script}: Missing Python 3 shebang line")
            
            if fix:
                # Fix shebang line
                if content.startswith("#!/usr/bin/env python"):
                    new_content = content.replace("#!/usr/bin/env python", "#!/usr/bin/env python3", 1)
                    with open(script, "w", encoding="utf-8") as f:
                        f.write(new_content)
                    logger.info(f"  Fixed shebang line in {script}")
        
        # Check for Python 2 specific syntax or imports
        py2_patterns = [
            r"from __future__ import",
            r"except \w+ as \w+, \w+:",
            r"print (?!\()",
            r"xrange\(",
            r"\.iteritems\(\)",
            r"\.iterkeys\(\)",
            r"\.itervalues\(\)",
        ]
        
        issues = []
        for pattern in py2_patterns:
            for match in re.finditer(pattern, content):
                line_num = content[:match.start()].count('\n') + 1
                issues.append((line_num, match.group(0)))
        
        if issues:
            issues_found += 1
            logger.info(f"Found potential Python 2 syntax in {script}:")
            for line_num, text in sorted(issues):
                logger.info(f"  Line {line_num}: {text}")
    
    return issues_found


def check_dependencies():
    """Check for required Python dependencies"""
    missing_deps = []
    outdated_deps = []
    
    for dep in REQUIRED_DEPENDENCIES:
        try:
            # Try to import the module
            spec = importlib.util.find_spec(dep)
            if spec is None:
                missing_deps.append(dep)
                continue
            
            # Check version (where possible)
            if dep == "numpy":
                import numpy
                if numpy.__version__.startswith("1.") and int(numpy.__version__.split(".")[1]) < 16:
                    outdated_deps.append(f"{dep} (version {numpy.__version__}, need 1.16+)")
            elif dep == "cython":
                import Cython
                if hasattr(Cython, "__version__") and Cython.__version__.startswith("0.") and int(Cython.__version__.split(".")[1]) < 29:
                    outdated_deps.append(f"{dep} (version {Cython.__version__}, need 0.29+)")
        
        except ImportError:
            missing_deps.append(dep)
    
    if missing_deps:
        logger.warning(f"Missing required dependencies: {', '.join(missing_deps)}")
        logger.info("Install with: pip install " + " ".join(missing_deps))
    
    if outdated_deps:
        logger.warning(f"Outdated dependencies: {', '.join(outdated_deps)}")
        logger.info("Update with: pip install --upgrade " + " ".join(dep.split()[0] for dep in outdated_deps))
    
    return len(missing_deps) + len(outdated_deps)


def run_cmake_test():
    """Run a minimal CMake configuration to test Python 3 detection"""
    logger.info("Testing CMake Python 3 detection...")
    
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create a minimal CMakeLists.txt
        with open(temp_path / "CMakeLists.txt", "w", encoding="utf-8") as f:
            f.write("""
cmake_minimum_required(VERSION 3.12)
project(TestPython3 LANGUAGES C CXX)

# Test Python 3 detection
find_package(Python3 COMPONENTS Interpreter Development)
if(Python3_FOUND)
    message(STATUS "Python3 found: ${Python3_EXECUTABLE}")
    message(STATUS "Python3 version: ${Python3_VERSION}")
    message(STATUS "Python3 libraries: ${Python3_LIBRARIES}")
    message(STATUS "Python3 include dirs: ${Python3_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "Python3 not found")
endif()

# Test NumPy detection
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "NumPy include dir: ${NUMPY_INCLUDE_DIR}")

# Test Cython detection
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import Cython; print(Cython.__version__)"
    OUTPUT_VARIABLE CYTHON_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(CYTHON_VERSION)
    message(STATUS "Cython version: ${CYTHON_VERSION}")
else()
    message(WARNING "Cython not found")
endif()
            """)
        
        # Run CMake
        try:
            result = subprocess.run(
                ["cmake", "."],
                cwd=temp_dir,
                capture_output=True,
                text=True,
                check=True
            )
            logger.info("CMake test successful:")
            for line in result.stdout.splitlines():
                if "Python" in line or "NumPy" in line or "Cython" in line:
                    logger.info(f"  {line.strip()}")
            return 0
        except subprocess.CalledProcessError as e:
            logger.error(f"CMake test failed: {e}")
            logger.error(e.stderr)
            return 1


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Verify Python 3 build system compatibility"
    )
    parser.add_argument(
        "--fix", action="store_true", help="Fix common issues automatically"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Show verbose output"
    )
    parser.add_argument(
        "--skip-test", action="store_true", help="Skip CMake test"
    )
    args = parser.parse_args()
    
    logger.info("Checking Python 3 build system compatibility...")
    
    # Check CMake files
    cmake_issues = check_cmake_files(args.fix)
    if cmake_issues == 0:
        logger.info("✅ CMake files are compatible with Python 3")
    else:
        logger.warning(f"⚠️ Found {cmake_issues} issues in CMake files")
    
    # Check installation scripts
    script_issues = check_install_scripts(args.fix)
    if script_issues == 0:
        logger.info("✅ Installation scripts are compatible with Python 3")
    else:
        logger.warning(f"⚠️ Found {script_issues} issues in installation scripts")
    
    # Check dependencies
    dep_issues = check_dependencies()
    if dep_issues == 0:
        logger.info("✅ All required dependencies are installed")
    else:
        logger.warning(f"⚠️ Found {dep_issues} dependency issues")
    
    # Run CMake test
    if not args.skip_test:
        if run_cmake_test() == 0:
            logger.info("✅ CMake Python 3 detection test passed")
        else:
            logger.error("❌ CMake Python 3 detection test failed")
    
    total_issues = cmake_issues + script_issues + dep_issues
    if total_issues == 0:
        logger.info("\n🎉 All checks passed! Build system is compatible with Python 3.")
        return 0
    else:
        logger.warning(f"\n⚠️ Found {total_issues} issues. Please address them before proceeding.")
        if args.fix:
            logger.info("Some issues were automatically fixed. Please run the verification again.")
        else:
            logger.info("Run with --fix to automatically fix common issues.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
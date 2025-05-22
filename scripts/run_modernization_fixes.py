#!/usr/bin/env python3
"""
Run pre-commit hooks on the entire codebase and fix Python 2 artifacts.

This script:
1. Installs pre-commit if not already installed
2. Runs all hooks on all files
3. Identifies remaining Python 2 artifacts that need manual fixing
"""

import subprocess
import sys


def check_pre_commit_installed():
    """Check if pre-commit is installed and install if needed."""
    try:
        subprocess.run(
            ["pre-commit", "--version"],
            check=True, 
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        print("✅ pre-commit is already installed")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("⚠️ pre-commit not found, installing...")
        try:
            subprocess.run(
                [sys.executable, "-m", "pip", "install", "pre-commit"],
                check=True
            )
            print("✅ pre-commit installed successfully")
            return True
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed to install pre-commit: {e}")
            return False


def setup_pre_commit_hooks():
    """Install pre-commit hooks in the repository."""
    try:
        subprocess.run(
            ["pre-commit", "install"],
            check=True
        )
        print("✅ pre-commit hooks installed")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install pre-commit hooks: {e}")
        return False


def run_pre_commit_on_all_files():
    """Run pre-commit on all files in the repository."""
    print("\n📋 Running pre-commit on all files (this may take a while)...")
    try:
        result = subprocess.run(
            ["pre-commit", "run", "--all-files"],
            capture_output=True,
            text=True
        )
        
        # Output the result
        print(result.stdout)
        if result.stderr:
            print("Errors:", result.stderr)
            
        if result.returncode == 0:
            print("✅ All pre-commit hooks passed")
            return True
        else:
            print("⚠️ Some pre-commit hooks failed, manual fixes needed")
            return False
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running pre-commit: {e}")
        return False


def run_specific_hooks():
    """Run specific pre-commit hooks that are most likely to fix Python 2 artifacts."""
    print("\n📋 Running specific hooks to fix Python 2 artifacts...")
    
    hooks = ["black", "ruff", "isort", "pyupgrade"]
    
    for hook in hooks:
        print(f"\nRunning {hook}...")
        try:
            subprocess.run(
                ["pre-commit", "run", hook, "--all-files"],
                check=False
            )
        except subprocess.CalledProcessError:
            print(f"⚠️ {hook} reported issues that need fixing")


def scan_for_python2_artifacts(directory="src/python"):
    """Scan for remaining Python 2 artifacts that need manual fixing."""
    print("\n📋 Scanning for remaining Python 2 artifacts...")
    
    artifacts = {
        "print statements": ["grep", "-r", "print ", directory],
        "xrange": ["grep", "-r", "xrange", directory],
        "unicode literals": ["grep", "-r", "unicode_literals", directory],
        "old style except": ["grep", "-r", "except .*,", directory],
        "old style super": ["grep", "-r", "super(", directory],
        "dict.iteritems": ["grep", "-r", ".iteritems", directory],
    }
    
    found_artifacts = False
    
    for name, command in artifacts.items():
        print(f"\nChecking for {name}...")
        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True
            )
            
            if result.stdout:
                found_artifacts = True
                print(f"Found {name}:")
                for line in result.stdout.splitlines():
                    print(f"  {line}")
            else:
                print(f"✅ No {name} found")
                
        except subprocess.CalledProcessError:
            # grep returns 1 if no matches found
            print(f"✅ No {name} found")
            
    return not found_artifacts


def main():
    """Main function."""
    print("🚀 Running Python 3 modernization fixes...")
    
    # Check if pre-commit is installed
    if not check_pre_commit_installed():
        print("❌ Cannot continue without pre-commit")
        return 1
        
    # Setup pre-commit hooks
    if not setup_pre_commit_hooks():
        print("❌ Failed to setup pre-commit hooks")
        return 1
    
    # Run specific hooks first to fix common issues
    run_specific_hooks()
    
    # Run all pre-commit hooks
    run_pre_commit_on_all_files()
    
    # Scan for remaining Python 2 artifacts
    all_clean = scan_for_python2_artifacts()
    
    if all_clean:
        print("\n✅ No Python 2 artifacts found - codebase is clean!")
    else:
        print("\n⚠️ Some Python 2 artifacts remain - manual fixing required")
    
    print("\n📋 Next steps:")
    print("1. Fix any remaining Python 2 artifacts")
    print("2. Run 'pre-commit run --all-files' to verify all hooks pass")
    print("3. Implement Python equivalents for C++ components")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

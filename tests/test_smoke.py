import pytest
import subprocess
import sys
from happy import __version__


def test_import_happy():
    # Package imports without error
    import happy
    assert hasattr(happy, '__version__')


def test_version_format():
    # Version should follow semantic versioning (at least major.minor.patch)
    assert isinstance(__version__, str)
    parts = __version__.split('.')
    assert len(parts) >= 3
    assert all(part.isdigit() for part in parts)


def test_cli_help_hap():
    # Ensure CLI entrypoint 'hap' runs and shows help
    result = subprocess.run([sys.executable, '-m', 'happy.cli.hap', '--help'],
                            capture_output=True, text=True)
    assert result.returncode == 0
    output = result.stdout + result.stderr
    assert 'usage' in output.lower()


def test_cli_help_som():
    # Ensure CLI entrypoint 'som' runs and shows help
    result = subprocess.run([sys.executable, '-m', 'happy.cli.som', '--help'],
                            capture_output=True, text=True)
    assert result.returncode == 0
    output = result.stdout + result.stderr
    assert 'usage' in output.lower()

"""
Happy: CLI entry point package for hap.py commands.
"""

# Removed deprecated subcommands cnx and ftx

# Extend package path to include parent directory, enabling imports of Haplo as happy.Haplo
import os  # noqa: E402

__path__.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Expose subcommands (wrapped in try/except to avoid import errors during testing)
try:
    from .hap import main as hap_main  # noqa: F401
except ImportError:
    hap_main = None
try:
    from .pre import main as pre_main  # noqa: F401
except ImportError:
    pre_main = None
try:
    from .qfy import main as qfy_main  # noqa: F401
except ImportError:
    qfy_main = None

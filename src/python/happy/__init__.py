"""
Happy: CLI entry point package for hap.py commands.
"""
from .cnx import main as cnx_main  # noqa: F401
from .ftx import updateArgs as ftx_updateArgs  # noqa: F401

# Expose subcommands
from .hap import main as hap_main  # noqa: F401
from .ovc import main as ovc_main  # noqa: F401
from .pre import updateArgs as pre_updateArgs  # noqa: F401
from .qfy import updateArgs as qfy_updateArgs  # noqa: F401

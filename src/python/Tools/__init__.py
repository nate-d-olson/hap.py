"""Tools helper namespace – environment helpers & assorted utilities.

This modernised ``__init__`` avoids heavy side-effects during *import* while
retaining backwards-compatibility with historical code expecting
``Tools.init()`` and a ``version`` attribute.
"""

from __future__ import annotations

import contextlib
import logging
import os
import shutil
import sys
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _pkg_version
from pathlib import Path
from typing import Optional

from . import vcfextract

LOGGER = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Version handling – prefer the installed package metadata.
# ---------------------------------------------------------------------------

try:
    __version__: str = _pkg_version("happy")
except PackageNotFoundError:  # fallback for editable installs
    __version__ = "unknown"

# Historical fallback (Haplo.version) – suppress ImportError if the generated
# file is not present in editable / PyPI installs.
with contextlib.suppress(ImportError):  # pragma: no cover – legacy path
    from Haplo import version as _haplo_version_module  # type: ignore

    __version__ = getattr(_haplo_version_module, "__version__", __version__)

# Re-export for callers that did ``from Tools import version``
version = __version__  # type: ignore  # allow attribute-style access


# ---------------------------------------------------------------------------
# Lightweight environment bootstrap
# ---------------------------------------------------------------------------

_DEFAULT_EXECUTABLES = (
    "blocksplit",
    "hapenum",
    "dipenum",
    "hapcmp",
    "bcftools",
    "samtools",
)


def _which(prog: str) -> str | None:
    """Return full path to *prog* if it is executable in ``PATH``."""

    return shutil.which(prog)


# return diagnostic dictionary
def init(verbose: bool = False, *, return_info: bool = False):  # type: ignore  # noqa: D401
    """Initialise helper environment (PATH tweaks, diagnostics).

    import-time side-effects have been removed; call this once at runtime when
    real work starts.

    Parameters
    ----------
    verbose
        If ``True`` emit INFO diagnostics.  Otherwise DEBUG.
    return_info
        When ``True`` a mapping with diagnostic data is returned; otherwise
        the function returns *None*.
    """

    log_level = logging.INFO if verbose else logging.DEBUG
    LOGGER.setLevel(log_level)

    diagnostics = {
        "added_legacy_bin": False,
        "missing_exec": [],
    }

    project_base = Path(__file__).resolve().parents[3]
    legacy_bin = project_base / "bin"
    if legacy_bin.is_dir():
        os.environ["PATH"] = f"{legacy_bin}:{os.environ.get('PATH', '')}"
        diagnostics["added_legacy_bin"] = True
        LOGGER.debug("Added legacy helper binaries dir to PATH: %s", legacy_bin)

    for exe in _DEFAULT_EXECUTABLES:
        if _which(exe) is None:
            diagnostics["missing_exec"].append(exe)
            LOGGER.info(
                "Executable %s not found in PATH; some features may be unavailable", exe
            )

    return diagnostics if return_info else None


# ---------------------------------------------------------------------------
# Convenience helpers preserved from the legacy implementation
# ---------------------------------------------------------------------------


def defaultReference() -> str | None:  # noqa: N802  (keep public API name)
    """Return path to default reference FASTA if present on system."""

    candidates = [os.environ.get(env) for env in ("HGREF", "HG19") if env in os.environ]
    candidates.append("/opt/hap.py-data/hg19.fa")
    for path in candidates:
        if path and Path(path).exists():
            return path
    return None


# Public re-exports expected by callers (mostly numpy/pandas helpers in other
# sub-modules).  We keep the import guarded so minimal installations without
# heavy deps remain functional.

with contextlib.suppress(ImportError):  # pragma: no cover – optional deps
    import pandas as pd  # type: ignore

    __all__ = ["pd"]

# Fallback when pandas is not installed.
__all__ = list(__all__) if "__all__" in globals() else []

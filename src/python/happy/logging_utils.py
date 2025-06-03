"""Centralised logging configuration for the happy tool-suite."""

from __future__ import annotations

import logging
import sys
from pathlib import Path


def setup_logging(
    *, verbose: bool = False, quiet: bool = False, log_file: str | None = None
) -> None:
    """Configure root logger according to CLI flags.

    Parameters
    ----------
    verbose
        When ``True`` set level to :pydata:`logging.DEBUG`.
    quiet
        Suppress INFO messages (level becomes ``WARNING``) – ignored if
        *verbose* is also ``True``.
    log_file
        Optional path to a file that receives the full DEBUG stream in
        addition to stderr.
    """

    level = logging.INFO
    if verbose:
        level = logging.DEBUG
    elif quiet:
        level = logging.WARNING

    fmt = "%(asctime)s %(levelname)-8s | %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"

    logging.basicConfig(level=level, format=fmt, datefmt=datefmt, stream=sys.stderr)

    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(fmt, datefmt))
        logging.getLogger().addHandler(file_handler)

"""Logging utilities for Python 3."""

import logging
from typing import Any


class LoggingWriter:
    """Helper class to write output to the Python logging system.

    This class provides a file-like interface for writing to Python's
    logging system, making it suitable for capturing output from other
    tools or redirecting sys.stdout/stderr.
    """

    def __init__(self, level: int = logging.INFO) -> None:
        """Initialize logging writer.

        Args:
            level: Logging level to use (e.g. logging.INFO, logging.ERROR)
        """
        self.level = level

    def write(self, message: str) -> None:
        """Write a message to the log.

        Args:
            message: The message to log
        """
        message = message.replace("\n", "")
        if message:  # Don't log empty messages
            logging.log(self.level, message)

    def flush(self) -> None:
        """Flush the output (no-op for logging)."""
        pass

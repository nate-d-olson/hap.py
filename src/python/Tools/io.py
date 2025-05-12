"""Python 3 compatible input/output utilities."""

from typing import Union, BinaryIO, TextIO, Optional
import os


def open_file(
    filename: str, mode: str = "r", encoding: Optional[str] = None
) -> Union[BinaryIO, TextIO]:
    """Open files in a Python 3 compatible way.

    This function ensures proper handling of text encodings and binary modes.

    Args:
        filename: Path to the file to open
        mode: File open mode ('r', 'w', 'rb', 'wb', etc.)
        encoding: Text encoding to use. Defaults to UTF-8 for text mode

    Returns:
        A file object in binary or text mode

    Raises:
        IOError: If the file cannot be opened
    """
    if encoding is None:
        encoding = "utf-8"

    if "b" in mode:
        return open(filename, mode)
    else:
        return open(filename, mode, encoding=encoding)

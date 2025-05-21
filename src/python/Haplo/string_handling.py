"""
String handling utilities for Python 3 compatibility in Cython modules.

This module provides utility functions for handling string encoding/decoding
at the C++/Python boundary, which is a common source of issues in Python 3
migrations of Cython code.
"""

from typing import Optional, Union, TypeVar, overload

# Define a type variable for the return type of ensure_text_io
T = TypeVar('T', str, bytes)


def ensure_str(text: Union[str, bytes, None]) -> Optional[str]:
    """
    Ensure input is a string (unicode in Python 3), converting from bytes if necessary.

    This is useful for C++ interfaces that might return bytes in Python 3.

    Args:
        text: Text to convert, can be str, bytes, or None

    Returns:
        Unicode string, or None if input was None
    """
    if text is None:
        return None
    if isinstance(text, bytes):
        return text.decode("utf-8", errors="replace")
    return str(text)


def ensure_bytes(text: Union[str, bytes, None]) -> Optional[bytes]:
    """
    Ensure input is bytes, converting from string if necessary.

    This is useful for C++ interfaces that expect bytes in Python 3.

    Args:
        text: Text to convert, can be str, bytes, or None

    Returns:
        Bytes object, or None if input was None
    """
    if text is None:
        return None
    if isinstance(text, str):
        return text.encode("utf-8")
    return bytes(text)


# Use overload to define the specific input/output type combinations
@overload
def ensure_text_io(text: str, file_mode: str) -> Union[str, bytes]:
    ...


@overload
def ensure_text_io(text: bytes, file_mode: str) -> Union[str, bytes]:
    ...


def ensure_text_io(text: Union[str, bytes], file_mode: str) -> Union[str, bytes]:
    """
    Ensure text has the right type for the given file mode.

    Args:
        text: Text content to check and potentially convert
        file_mode: File mode ('r', 'rb', 'w', 'wb', etc.)

    Returns:
        Text in appropriate format for mode (str for text mode, bytes for binary)
    """
    is_binary_mode = "b" in file_mode

    if is_binary_mode and isinstance(text, str):
        return text.encode("utf-8")
    elif not is_binary_mode and isinstance(text, bytes):
        return text.decode("utf-8", errors="replace")

    return text

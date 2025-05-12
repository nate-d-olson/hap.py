"""Python 3 compatible BGZip file handling."""

import os
import subprocess
from typing import Optional, Union, BinaryIO, Any


class BGZipFile:
    """BGZip file helper with context manager support and enhanced Python 3 compatibility."""

    def __init__(self, filename: str, force: bool = False) -> None:
        """Initialize a BGZip file handler.
        Args:
            filename: Path to the output file
            force: Whether to overwrite existing files
        Raises:
            OSError: If bgzip executable is not found
            IOError: If file exists and force=False
            RuntimeError: If bgzip process fails to start
        """
        if os.path.exists(filename) and not force:
            raise IOError(f"File {filename} exists, use force=True to overwrite")

        if not any(
            os.access(os.path.join(p, "bgzip"), os.X_OK)
            for p in os.environ["PATH"].split(os.pathsep)
        ):
            raise OSError("bgzip executable not found in PATH")

        self.write_file: Optional[BinaryIO] = None
        self.zip_pipe: Optional[subprocess.Popen[bytes]] = None
        self._closed: bool = False
        self.name: str = filename

        try:
            self.write_file = open(filename, "wb")
            self.zip_pipe = subprocess.Popen(
                ["bgzip", "-f"],
                stdin=subprocess.PIPE,
                stdout=self.write_file,
                stderr=subprocess.PIPE,
                universal_newlines=False,
                bufsize=-1,
            )
            if self.zip_pipe.stdin is None:
                raise RuntimeError("Failed to open pipe to bgzip")
        except Exception as e:
            if self.write_file:
                self.write_file.close()
            raise RuntimeError(f"Failed to start bgzip process: {str(e)}") from e

    def __enter__(self) -> "BGZipFile":
        """Context manager entry."""
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> bool:
        """Context manager exit."""
        self.close()
        return False

    def close(self) -> None:
        """Close the BGZip file."""
        if not self._closed:
            try:
                if self.zip_pipe and self.zip_pipe.stdin:
                    self.zip_pipe.stdin.flush()
                    self.zip_pipe.stdin.close()
                if self.zip_pipe:
                    returncode = self.zip_pipe.wait(timeout=10)
                    if returncode != 0:
                        stderr = (
                            self.zip_pipe.stderr.read() if self.zip_pipe.stderr else b""
                        )
                        stderr = stderr.decode("utf-8", errors="replace")
                        raise RuntimeError(
                            f"bgzip failed with code {returncode}: {stderr}"
                        )
            finally:
                if self.write_file:
                    self.write_file.flush()
                    self.write_file.close()
                self._closed = True

    def write(self, data: Union[str, bytes]) -> None:
        """Write data to the BGZip file.
        Args:
            data: The data to write (either string or bytes)
        Raises:
            ValueError: If file is closed or pipe is not initialized
            TypeError: If data is neither string nor bytes
        """
        if self._closed:
            raise ValueError("I/O operation on closed file")

        if isinstance(data, str):
            data = data.encode("utf-8")
        elif not isinstance(data, bytes):
            raise TypeError(
                f"write() argument must be str or bytes, not {type(data).__name__}"
            )

        if self.zip_pipe and self.zip_pipe.stdin:
            self.zip_pipe.stdin.write(data)
        else:
            raise ValueError("BGZip pipe is not properly initialized")

    def flush(self) -> None:
        """Flush the BGZip file buffers."""
        if not self._closed and self.zip_pipe and self.zip_pipe.stdin:
            self.zip_pipe.stdin.flush()
            if self.write_file:
                self.write_file.flush()

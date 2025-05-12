"""Python 3 compatible BGZip file handling."""

import os
import subprocess
from typing import Optional, Union, BinaryIO


class BGZipFile:
    """BGZip file helper with context manager support and enhanced Python 3 compatibility."""

    def __init__(self, filename: str, force: bool = False) -> None:
        """Make a subprocess for bgzip.

        Args:
            filename: name of the output file
            force: true to overwrite if file exists

        Raises:
            OSError: if bgzip is not available
            IOError: if file exists and force=False
            RuntimeError: if bgzip process fails to start
        """
        if os.path.exists(filename) and not force:
            raise IOError(f"File {filename} exists, use force=True to overwrite")

        # Check if bgzip is available
        if not any(
            os.access(os.path.join(p, "bgzip"), os.X_OK)
            for p in os.environ["PATH"].split(os.pathsep)
        ):
            raise OSError("bgzip executable not found in PATH")

        self.write_file: Optional[IOBase] = None
        self.zip_pipe: Optional[subprocess.Popen] = None
        self._closed: bool = False
        self.name: str = filename

        try:
            self.write_file = open(filename, "wb")
            self.zip_pipe = subprocess.Popen(
                ["bgzip", "-f"],
                stdin=subprocess.PIPE,
                stdout=self.write_file,
                stderr=subprocess.PIPE,
                universal_newlines=False,  # Ensure binary mode
            )
            if self.zip_pipe.stdin is None:
                raise RuntimeError("Failed to open pipe to bgzip")
        except Exception as e:
            if self.write_file:
                self.write_file.close()
            raise RuntimeError(f"Failed to start bgzip process: {str(e)}") from e

    def __enter__(self) -> "BGZipFile":
        """Context manager entry.

        Returns:
            self: The BGZipFile instance
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        """Context manager exit.

        Args:
            exc_type: Exception type if an exception occurred
            exc_val: Exception value if an exception occurred
            exc_tb: Exception traceback if an exception occurred

        Returns:
            bool: False to not suppress exceptions
        """
        self.close()
        return False

    def close(self) -> None:
        """Close the BGZip file.

        Raises:
            RuntimeError: if bgzip process returns non-zero exit code
        """
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
            data: Data to write (str or bytes)

        Raises:
            ValueError: if file is closed or pipe is not initialized
            TypeError: if data is neither str nor bytes
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
        """Flush the BGZip file.

        This method flushes both the bgzip pipe and the underlying file.
        """
        if not self._closed and self.zip_pipe and self.zip_pipe.stdin:
            self.zip_pipe.stdin.flush()
            if self.write_file:
                self.write_file.flush()

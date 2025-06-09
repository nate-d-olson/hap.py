"""RTG Tools manager for automated setup and configuration."""

import logging
import os
import shutil
import subprocess
import tempfile
import urllib.request
import zipfile
from pathlib import Path

RTG_VERSION = "3.13"
RTG_DOWNLOAD_URL = f"https://github.com/RealTimeGenomics/rtg-tools/releases/download/{RTG_VERSION}/rtg-tools-{RTG_VERSION}-nojre.zip"


class RTGManager:
    """Manager for RTG Tools installation and execution."""

    def __init__(self):
        self.package_dir = Path(__file__).parent.parent  # hap_py directory
        self.external_dir = self.package_dir / "external"
        self.rtg_dir = self.external_dir / f"rtg-tools-{RTG_VERSION}"
        self.rtg_executable = self.rtg_dir / "rtg"

    def is_installed(self) -> bool:
        """Check if RTG tools is installed and executable."""
        return self.rtg_executable.exists() and os.access(self.rtg_executable, os.X_OK)

    def get_rtg_path(self) -> str:
        """Get path to RTG executable, installing if necessary."""
        # Check for environment variable override
        if "RTG_PATH" in os.environ and Path(os.environ["RTG_PATH"]).exists():
            logging.info(
                f"Using RTG tools from environment variable: {os.environ['RTG_PATH']}"
            )
            return os.environ["RTG_PATH"]

        # Check if RTG is already installed in our package
        if self.is_installed():
            return str(self.rtg_executable)

        # If not installed, install it now
        logging.info("RTG Tools not found, installing automatically...")
        self.install_rtg()

        if not self.is_installed():
            raise RuntimeError(
                "Failed to install RTG Tools. Please install manually and set RTG_PATH environment variable."
            )

        return str(self.rtg_executable)

    def install_rtg(self) -> None:
        """Download and install RTG tools."""
        # Create external directory if it doesn't exist
        self.external_dir.mkdir(parents=True, exist_ok=True)

        # Download RTG tools
        with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as temp_file:
            logging.info(f"Downloading RTG Tools v{RTG_VERSION}...")
            urllib.request.urlretrieve(RTG_DOWNLOAD_URL, temp_file.name)

            # Extract RTG tools
            logging.info("Extracting RTG Tools...")
            with zipfile.ZipFile(temp_file.name, "r") as zip_ref:
                zip_ref.extractall(self.external_dir)

        # Make RTG executable
        os.chmod(self.rtg_executable, 0o755)

        # Create auto-accept config to avoid interactive prompt
        self.create_rtg_config()

        logging.info(f"RTG Tools installed successfully to {self.rtg_dir}")

    def create_rtg_config(self) -> None:
        """Create RTG config file to avoid interactive setup."""
        config_path = self.rtg_dir / "rtg.cfg"
        with open(config_path, "w") as f:
            f.write(
                """# Auto-generated RTG config
RTG_TALKBACK=false
RTG_USAGE=false
RTG_JAVA=java
RTG_JAR={}
""".format(
                    self.rtg_dir / "RTG.jar"
                )
            )

    def run_rtg_command(self, *args) -> subprocess.CompletedProcess:
        """Run RTG command with proper environment setup."""
        rtg_path = self.get_rtg_path()
        cmd = [rtg_path] + list(args)

        # Set Java environment if needed
        env = os.environ.copy()
        if "JAVA_HOME" not in env and shutil.which("java") is None:
            raise RuntimeError(
                "Java is required for RTG Tools but not found in PATH. Please install Java."
            )

        return subprocess.run(cmd, env=env, check=True)


# Global instance for easy access
rtg_manager = RTGManager()


def get_rtg_path() -> str:
    """Get path to RTG executable for use in other modules."""
    return rtg_manager.get_rtg_path()


def run_rtg_vcfeval(
    truth_vcf, query_vcf, reference_sdf, output_dir, *args
) -> subprocess.CompletedProcess:
    """Run RTG vcfeval with the given parameters."""
    cmd_args = [
        "vcfeval",
        "-b",
        truth_vcf,
        "-c",
        query_vcf,
        "-t",
        reference_sdf,
        "-o",
        output_dir,
    ]
    cmd_args.extend(args)
    return rtg_manager.run_rtg_command(*cmd_args)

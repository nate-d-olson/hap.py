"""Installation utilities for hap.py package."""


def install():
    """Run post-install setup tasks."""
    try:
        # Import the RTG manager to trigger initial setup
        from ..external.rtg_manager import rtg_manager

        # Pre-download RTG tools during install if not present
        if not rtg_manager.is_installed():
            rtg_manager.install_rtg()
            print(f"RTG Tools installed successfully to {rtg_manager.rtg_dir}")
        else:
            print(f"Using existing RTG Tools at {rtg_manager.rtg_executable}")
    except Exception as e:
        print(f"Warning: Could not set up RTG Tools automatically: {str(e)}")
        print(
            "You may need to install RTG Tools manually and set RTG_PATH environment variable."
        )


if __name__ == "__main__":
    install()

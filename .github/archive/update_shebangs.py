#!/usr/bin/env python3
# update_shebangs.py - Update all shebang lines to use Python 3
#
# This script recursively updates the shebang lines in all Python files
# to ensure they use Python 3. It also adds encoding headers where missing.

import os
import sys
from pathlib import Path


def update_shebang(file_path):
    """Update shebang line in a file to use Python 3"""
    with open(file_path, encoding="utf-8", errors="replace") as f:
        content = f.read()

    # Check if shebang needs updating
    if content.startswith("#!/usr/bin/env python") and not content.startswith(
        "#!/usr/bin/env python3"
    ):
        content = content.replace("#!/usr/bin/env python", "#!/usr/bin/env python3", 1)
        # Add encoding header if missing
        if (
            "coding: utf-8" not in content[:100]
            and "# -*- coding:" not in content[:100]
        ):
            lines = content.split("\n")
            if len(lines) > 1 and lines[1].startswith("#"):
                lines.insert(2, "# -*- coding: utf-8 -*-")
            else:
                lines.insert(1, "# -*- coding: utf-8 -*-")
            content = "\n".join(lines)

        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)
        return True
    return False


def find_and_update_python_files(root_dir):
    """Recursively find and update Python files"""
    updated_files = []
    root_path = Path(root_dir)

    for path in root_path.glob("**/*.py"):
        if update_shebang(path):
            updated_files.append(str(path))

    return updated_files


if __name__ == "__main__":
    if len(sys.argv) > 1:
        root_dir = sys.argv[1]
    else:
        root_dir = os.path.dirname(os.path.abspath(__file__))

    updated_files = find_and_update_python_files(root_dir)

    print(f"Updated {len(updated_files)} files to use Python 3 shebang line.")
    if updated_files:
        print("Files updated:")
        for file in updated_files[:10]:
            print(f"  - {file}")
        if len(updated_files) > 10:
            print(f"  - ...and {len(updated_files) - 10} more files")

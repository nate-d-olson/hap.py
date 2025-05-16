#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
#
# 31/05/2017
#
# Collect session and run information
#

"""Module for collecting and reporting session and runtime information."""

import os
import platform
import sys
import time
from typing import Any, Dict

import Tools


def sessionInfo() -> Dict[str, Any]:
    """Return a dictionary with session and run information.

    Collects system information, Python environment details, and installed packages.

    Returns:
        Dict containing various session information and runtime details.
    """

    version = f"{Tools.version}"

    result: Dict[str, Any] = {
        "name": os.path.basename(sys.argv[0]),
        "timestamp": time.strftime("%a %b %d %X %Y"),
        "version": version,
        "runInfo": [{"key": "commandline", "value": " ".join(sys.argv)}],
        "uname": " / ".join(platform.uname()),
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "metadata": {
            "required": {
                "id": "haplotypes",
                "version": version,
                "module": f"{os.path.basename(sys.argv[0])}",
                "description": f"{sys.argv[0]} generated this JSON file via command line {' '.join(sys.argv)}",
            }
        },
        "environment": {str(k): str(v) for k, v in list(os.environ.items())},
    }

    # platform.dist() is removed in Python 3.8+, use distro module if available
    try:
        import distro

        result["dist"] = " / ".join(distro.linux_distribution())
    except ImportError:
        try:
            # Fallback for Python < 3.8
            result["dist"] = " / ".join(platform.dist())
        except AttributeError:
            result["dist"] = "Not available"

    # Mac version information
    mac_ver_info = platform.mac_ver()
    if mac_ver_info[0]:
        result["mac_ver"] = f"{mac_ver_info[0]} / {mac_ver_info[2]}"

    result["python_prefix"] = sys.prefix

    # Check for virtual environment
    # Python 3 uses sys.base_prefix instead of sys.real_prefix
    if hasattr(sys, "real_prefix") or (
        hasattr(sys, "base_prefix") and sys.base_prefix != sys.prefix
    ):
        result["python_virtualenv"] = True
        result["python_real_prefix"] = getattr(sys, "real_prefix", sys.base_prefix)

    # Collect system resources information
    try:
        import psutil

        result["cpus"] = psutil.cpu_count()
        result["logical_cpus"] = psutil.cpu_count(logical=True)
        cpu_freq = psutil.cpu_freq()
        if cpu_freq:
            result["cpu_freq"] = cpu_freq._asdict()
        result["memory"] = dict(psutil.virtual_memory()._asdict())
    except ImportError:
        pass

    # Collect installed packages
    try:
        # In Python 3, we need to use importlib.metadata instead of pip
        try:
            from importlib.metadata import distributions

            pip_packages = [str(dist) for dist in distributions()]
        except ImportError:
            # Fallback for Python < 3.8
            import pkg_resources

            pip_packages = [str(d) for d in pkg_resources.working_set]

        result["pip_packages"] = pip_packages
    except ImportError:
        pass

    return result

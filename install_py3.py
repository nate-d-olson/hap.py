#!/usr/bin/env python3
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
# hap.py installer (Python 3 compatible):
#
# * checks dependencies
# * builds code
# * makes virtualenv

import argparse
import fnmatch
import glob
import multiprocessing
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request as urllib_request


def check_python_version():
    """Check if the python version is sufficient"""
    if sys.version_info < (3, 6):
        raise Exception("You will need to run this with Python >= 3.6")


def create_python_environment(source_dir, args):
    """Create a Python runtime environment
    :return: shebang with path to the python executable
    """
    interp = args.python_interp
    pyver_output = subprocess.check_output(
        f"{interp} -c \"import sys; print(','.join(map(str, list(sys.version_info[0:3]))))\"",
        shell=True,
    ).strip()

    # Convert bytes to str if needed in Python 3
    if isinstance(pyver_output, bytes):
        pyver_output = pyver_output.decode("utf-8")

    pyver = tuple(map(int, pyver_output.split(",")))

    if pyver < (3, 6):
        raise Exception("Python >= 3.6 is required for installation.")

    # system python -- just return interp
    if args.python == "system":
        return "#!" + interp

    if not args.python_venv_dir:
        raise Exception("Please specify a virtualenv target installation directory.")

    if args.python_venv_dir_force:
        try:
            shutil.rmtree(args.python_venv_dir)
        except:
            pass

    if os.path.exists(args.python_venv_dir) and not args.python_venv_dir_update:
        raise Exception("The virtual environment directory already exists.")

    # Use venv module from Python 3 standard library instead of virtualenv
    virtualenv_cmd = [interp, "-m", "venv", args.python_venv_dir]
    print(f"Creating virtual environment: {' '.join(virtualenv_cmd)}", file=sys.stderr)
    subprocess.check_call(virtualenv_cmd)

    # install requirements
    ve_python = os.path.join(args.python_venv_dir, "bin", "python")
    ve_pip = os.path.join(args.python_venv_dir, "bin", "pip")

    # Install pip-tools for better dependency management
    cmds = [ve_pip, "install", "-U", "pip", "setuptools", "wheel"]
    print(f"Updating pip: {' '.join(cmds)}", file=sys.stderr)
    subprocess.check_call(cmds)

    # Install requirements from file
    requirements_file = os.path.join(source_dir, "happy.requirements.txt")
    cmds = [ve_pip, "install", "--no-cache-dir", "-r", requirements_file]
    print(f"Installing requirements: {' '.join(cmds)}", file=sys.stderr)
    subprocess.check_call(cmds)

    # Return the shebang for Python scripts
    return "#!" + os.path.join(args.python_venv_dir, "bin", "python")

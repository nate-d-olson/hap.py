# coding=utf-8
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

import copy
import os
import platform
import sys
import time

import distro  # Use distro package instead of platform.dist()
import Tools


def sessionInfo():
    """ Return a dictionary with session and run information
    """

    version = "%s" % Tools.version

    # Get Linux distribution info safely using distro module (platform.dist() is deprecated)
    try:
        linux_dist = " / ".join([distro.name(), distro.version(), distro.codename()])
    except:
        linux_dist = "Unknown"

    result = {'name': os.path.basename(sys.argv[0]),
              'timestamp': time.strftime("%a %b %d %X %Y"),
              'version': version,
              'runInfo': [{"key": "commandline", "value": " ".join(sys.argv)}],
              'uname': " / ".join(platform.uname()),
              'dist': linux_dist,  # Use new distro info
              'mac_ver': " / ".join([platform.mac_ver()[0], platform.mac_ver()[2]]),
              'python_implementation': platform.python_implementation(),
              'python_version': platform.python_version(),
              'metadata': {
                  "required": {
                      "id": "haplotypes",
                      'version': version,
                      "module": "%s" % os.path.basename(sys.argv[0]),
                      "description": "%s generated this JSON file via command line %s" % (
                          sys.argv[0], " ".join(sys.argv))}},
              'environment': {str(k): str(os.environ[k]) for k in list(os.environ.keys())}}

    result["python_prefix"] = sys.prefix
    if hasattr(sys, 'real_prefix'):
        result["python_virtualenv"] = True
        result["python_real_prefix"] = sys.real_prefix
    elif hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix:
        # Python 3 virtual environments use base_prefix instead of real_prefix
        result["python_virtualenv"] = True
        result["python_base_prefix"] = sys.base_prefix

    try:
        import psutil
        result["cpus"] = psutil.cpu_count()
        result["logical_cpus"] = psutil.cpu_count(True)
        result["cpu_freq"] = psutil.cpu_freq()
        result["memory"] = dict(psutil.virtual_memory().__dict__)
    except:
        pass

    try:
        # Use importlib-metadata for getting package info in Python 3
        try:
            from importlib.metadata import distributions
            pip_packages = [str(dist) for dist in distributions()]
        except ImportError:
            # Fall back to pip for older Python versions
            import pip
            try:
                # pip 10 and later
                from pip._internal.utils.misc import get_installed_distributions
                pip_packages = [str(i) for i in get_installed_distributions(local_only=True)]
            except ImportError:
                # pip 9 and earlier
                pip_packages = [str(i) for i in pip.get_installed_distributions(local_only=True)]
        
        result["pip_packages"] = pip_packages
    except Exception as e:
        result["pip_error"] = str(e)
        pass

    return result


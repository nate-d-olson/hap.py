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
    requirements_file = os.path.join(source_dir, "happy.requirements.py3.txt")
    if not os.path.exists(requirements_file):
        requirements_file = os.path.join(source_dir, "happy.requirements.txt")
        print(
            "Warning: Using Python 2 requirements file. Consider creating happy.requirements.py3.txt",
            file=sys.stderr,
        )

    cmds = [ve_pip, "install", "--no-cache-dir", "-r", requirements_file]
    print(f"Installing requirements: {' '.join(cmds)}", file=sys.stderr)
    subprocess.check_call(cmds)

    # Return the shebang for Python scripts
    return "#!" + os.path.join(args.python_venv_dir, "bin", "python")


def replace_shebang(filename, shebang):
    """Replace shebang line / reheader script files"""
    print(f"Fixing shebang line in {filename}", file=sys.stderr)

    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()

    with open(filename, "w", encoding="utf-8") as f:
        removed = False
        print(shebang, file=f)
        for i, l in enumerate(lines):
            if not removed and l.startswith("#!") and i < 10:
                removed = True
            else:
                f.write(l)


def build_haplotypes(source_dir, build_dir, args):
    if args.boost:
        boost_prefix = f"BOOST_ROOT={args.boost} "
    else:
        boost_prefix = ""

    config_command = (
        f"{source_dir}/configure.sh {args.configuration} {args.setup} {args.targetdir}"
    )

    if args.sge:
        config_command += " -DUSE_SGE=ON"

    if args.build_rtgtools:
        config_command += " -DBUILD_VCFEVAL=ON"
        if args.rtgtools_wrapper:
            if not os.path.exists(args.rtgtools_wrapper):
                raise Exception(
                    f"RTG-tools wrapper {args.rtgtools_wrapper} doesn't exist."
                )
            wrapper_path = os.path.abspath(args.rtgtools_wrapper).replace(' ', '\\ ')
            config_command += f"-DVCFEVAL_WRAPPER={wrapper_path}"

    to_run = f"{boost_prefix}cd {build_dir} && {boost_prefix} {config_command}"
    print(to_run, file=sys.stderr)
    subprocess.check_call(to_run, shell=True)

    setupscript = ""
    if args.setup != "auto":
        setupscript = (
            f" . {os.path.join(source_dir, 'src', 'sh', args.setup + '-setup.sh')} && "
        )

    setupscript += boost_prefix

    to_run = f"{setupscript}cd {build_dir} && {setupscript} make -j{args.processes}"
    print(to_run, file=sys.stderr)
    subprocess.check_call(to_run, shell=True)

    to_run = (
        f"{setupscript}cd {build_dir} && {setupscript} make -j{args.processes} install"
    )
    print(to_run, file=sys.stderr)
    subprocess.check_call(to_run, shell=True)


def test_haplotypes(source_dir, python_shebang, args):
    """Run the unit + integration tests"""
    to_run = f"cd {args.targetdir} && {os.path.join(source_dir, 'src', 'sh', 'run_tests.sh')}"
    print(to_run, file=sys.stderr)
    os.environ["PYTHON"] = python_shebang[2:]
    subprocess.check_call(to_run, shell=True)


def main():
    check_python_version()

    source_dir = os.path.abspath(os.path.dirname(__file__))

    parser = argparse.ArgumentParser("hap.py installer (Python 3)")
    parser.add_argument("targetdir", help="Target installation directory")

    parser.add_argument(
        "--sge-mode",
        dest="sge",
        action="store_true",
        default=False,
        help="Enable SGE mode, which will require an additional command "
        'line option "--force-interactive" to run interactively.',
    )

    parser.add_argument(
        "--python",
        dest="python",
        choices=["system", "virtualenv"],
        default="system",
        help="Which Python to use in the installation. 'virtualenv' "
        "will create a virtual environment in the folder "
        "specified with --python-virtualenv-dir",
    )

    parser.add_argument(
        "--python-interpreter",
        dest="python_interp",
        default=sys.executable,
        help="Python interpreter to use for the installed hap.py.",
    )

    parser.add_argument(
        "--python-virtualenv-update",
        dest="python_venv_dir_update",
        default=False,
        action="store_true",
        help="Update virtualenv if it already exists.",
    )

    parser.add_argument(
        "--python-virtualenv-force",
        dest="python_venv_dir_force",
        default=False,
        action="store_true",
        help="Force creating a virtualenv even if the target directory"
        " already exists. USE WITH CARE, THIS WILL REMOVE THE "
        "VIRTUALENV DIRECTORY!",
    )

    parser.add_argument(
        "--python-virtualenv-dir",
        dest="python_venv_dir",
        default="",
        help="Directory to install the virtualenv in.",
    )

    # C++ compile options
    setups = [
        os.path.basename(x).replace("-setup.sh", "")
        for x in glob.glob(os.path.join(source_dir, "src", "sh", "*-setup.sh"))
    ]

    setups.insert(0, "auto")

    parser.add_argument(
        "--configuration",
        dest="configuration",
        choices=["Debug", "Release", "RelWithDebInfo", "install"],
        default="Release",
        help="Build configuration (use Release if unsure).",
    )

    parser.add_argument(
        "--setup",
        dest="setup",
        choices=setups,
        default="auto",
        help="Build setup (or auto to use system-wide packages).",
    )

    parser.add_argument(
        "--boost-root", dest="boost", help="Where to find Boost.", default=""
    )

    parser.add_argument(
        "--scratch-path", dest="scratch_path", help="Where to build.", default="/tmp"
    )

    parser.add_argument(
        "--keep-scratch",
        dest="keep_scratch",
        help="Keep the scratch folder.",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "--with-rtgtools",
        dest="build_rtgtools",
        default=False,
        action="store_true",
        help="Get and build rtgtools. You need to have Java and Ant for this.",
    )

    parser.add_argument(
        "--rtgtools-wrapper",
        dest="rtgtools_wrapper",
        default=None,
        help="Wrapper script for rtgtools. This is optional, it is useful "
        "when the default version of Java must be replaced / the environment "
        "needs updating. There is an example in src/sh/rtg-wrapper.sh.",
    )

    parser.add_argument(
        "--build-processes",
        dest="processes",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of parallel processes to use for building.",
    )

    parser.add_argument(
        "--no-tests",
        dest="run_tests",
        default=True,
        action="store_false",
        help="Disable unit tests",
    )

    parser.add_argument(
        "--build-externals-only",
        dest="build_externals_only",
        default=False,
        action="store_true",
        help="Only build external dependencies",
    )

    args = parser.parse_args()

    args.targetdir = os.path.abspath(args.targetdir)

    if args.python == "virtualenv" and not args.python_venv_dir:
        raise Exception("Please specify a virtualenv target installation directory.")

    if args.python_venv_dir and not args.python_venv_dir.startswith("/"):
        args.python_venv_dir = args.targetdir

    if "LD_LIBRARY_PATH" in os.environ or "DYLD_LIBRARY_PATH" in os.environ:
        print(
            "WARNING: You have (DY)LD_LIBRARY_PATH set. Make sure these libraries are accessible "
            "in the same environment you will run in.",
            file=sys.stderr,
        )

    # fix dynamic linking
    if "LD_LIBRARY_PATH" in os.environ:
        os.environ["LD_RUN_PATH"] = os.environ["LD_LIBRARY_PATH"]

    # check/make Python environment
    python_shebang = create_python_environment(source_dir, args)

    if args.boost and not os.path.exists(args.boost):
        raise Exception("Boost directory doesn't exist.")

    # Check if external dependencies need to be built separately
    if args.build_externals_only:
        ext_script = os.path.join(source_dir, "external", "make_dependencies_py3.sh")
        if not os.path.exists(ext_script):
            ext_script = os.path.join(source_dir, "external", "make_dependencies.sh")
            print(
                "Warning: Using Python 2 external dependencies script. Consider creating make_dependencies_py3.sh",
                file=sys.stderr,
            )

        ext_command = f"cd {source_dir} && {ext_script} {args.targetdir}"
        print(f"Building external dependencies: {ext_command}", file=sys.stderr)
        subprocess.check_call(ext_command, shell=True)
        print("External dependencies built successfully", file=sys.stderr)
        return

    # build hap.py
    build_dir = tempfile.mkdtemp(prefix="build", dir=args.scratch_path)
    try:
        build_haplotypes(source_dir, build_dir, args)
    finally:
        if not args.keep_scratch:
            try:
                shutil.rmtree(build_dir)
            except Exception as e:
                print(
                    f"Warning: Failed to remove build directory: {e}", file=sys.stderr
                )

    # reheader Python files
    for root, _, filenames in os.walk(args.targetdir):
        for filename in fnmatch.filter(filenames, "*.py"):
            replace_shebang(os.path.join(root, filename), python_shebang)

    # Create symlink to Python 3 version of scripts if needed
    for py3_script in glob.glob(os.path.join(source_dir, "src", "python", "*.py.py3")):
        base_script = os.path.basename(py3_script[:-4])  # Remove .py3
        target_script = os.path.join(args.targetdir, "bin", base_script)
        if os.path.exists(target_script):
            print(
                f"Creating symlink for Python 3 version of {base_script}",
                file=sys.stderr,
            )
            os.rename(target_script, f"{target_script}.py2.bak")
            shutil.copy2(py3_script, target_script)
            replace_shebang(target_script, python_shebang)

    if args.run_tests:
        test_haplotypes(source_dir, python_shebang, args)


if __name__ == "__main__":
    main()

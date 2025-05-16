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

import os
import errno
import logging
import subprocess
from datetime import date
from typing import Optional, List, Union, IO

# noinspection PyUnresolvedReferences
try:
    import Haplo.version as vs  # pylint: disable=E0401,E0611
    version = vs.__version__
    has_sge = vs.has_sge
except ImportError:
    logging.warning("No version found. Please follow the installation instructions.")
    version = "unknown"
    has_sge = False


def defaultReference() -> Optional[str]:
    """Return the default reference file path if one is found.

    Returns:
        Path to the default reference file or None if not found.
    """
    to_try = ['/opt/hap.py-data/hg19.fa']
    try:
        to_try.insert(0, os.environ["HGREF"])
    except KeyError:
        pass

    try:
        to_try.insert(0, os.environ["HG19"])
    except KeyError:
        pass

    for x in to_try:
        if os.path.exists(x):
            return x
    logging.warning("No reference file found at default locations. You can set the environment"
                 " variable 'HGREF' or 'HG19' to point to a suitable Fasta file.")
    return None


def which(program: str) -> Optional[str]:
    """Find an executable in PATH.

    Args:
        program: Name of the program to find

    Returns:
        Full path to the executable or None if not found
    """
    def is_exe(xfpath):
        return os.path.isfile(xfpath) and os.access(xfpath, os.X_OK)

    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def init() -> None:
    """Initialize the environment with required paths and dependencies."""
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)

    base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
    paths = ["bin"]

    for p in paths:
        pp = os.path.join(base, p)
        if not os.path.exists(pp):
            raise Exception(f"Dependency path {pp} not found")
        os.environ["PATH"] = pp + os.pathsep + os.environ["PATH"]

    executables = [
        "blocksplit",
        "hapenum",
        "dipenum",
        "hapcmp",
        "xcmp",
        "bcftools",
        "samtools",
    ]

    for x in executables:
        if not which(x):
            raise Exception(f"Dependency {x} not found")

    os.environ['DYLD_LIBRARY_PATH'] = os.path.join(base, "lib")
    os.environ['LD_LIBRARY_PATH'] = os.path.join(base, "lib")


init()

# safely import here
# noinspection PyUnresolvedReferences
import pysam  # noqa: F401
# noinspection PyUnresolvedReferences
import pandas  # noqa: F401


class LoggingWriter:
    """ Helper class to write tracebacks to log file
    """
    def __init__(self, level: int):
        self.level = level

    def write(self, message: str) -> None:
        message = message.replace("\n", "")
        if message:
            logging.log(self.level, message)


def writeVCFHeader(filelike: IO, extrainfo: Union[str, List[str]] = "", chrprefix: str = "chr") -> None:
    """ Write a VCF header

    Args:
        filelike: File-like object to write to
        extrainfo: Additional INFO fields as string or list
        chrprefix: Prefix to add to chromosome names
    """
    header = ["CHROM",
              "POS",
              "ID",
              "REF",
              "ALT",
              "QUAL",
              "FILTER",
              "INFO",
              "FORMAT",
              "SIMPLE"]

    infos = ['##fileformat=VCFv4.1',
             '##reference=hg19',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']

    if extrainfo:
        if isinstance(extrainfo, list):
            infos += extrainfo
        else:
            infos += extrainfo.split("\n")

    contigs = [["1", "249250621"], ["2", "243199373"], ["3", "198022430"], ["4", "191154276"],
               ["5", "180915260"], ["6", "171115067"], ["7", "159138663"], ["8", "146364022"],
               ["9", "141213431"], ["10", "135534747"], ["11", "135006516"],
               ["12", "133851895"], ["13", "115169878"], ["14", "107349540"], ["15", "102531392"],
               ["16", "90354753"], ["17", "81195210"], ["18", "78077248"], ["19", "59128983"],
               ["20", "63025520"], ["21", "48129895"], ["22", "51304566"], ["X", "155270560"]]

    meta = [f"##fileDate={date.today().isoformat()}",
            "##source=HaploCompare",
            f"##source_version={version}"]

    for i in infos:
        filelike.write(i + "\n")
    for c in contigs:
        filelike.write(f'##contig=<ID={chrprefix}{c[0]},length={c[1]}>\n')
    for i in meta:
        filelike.write(i + "\n")

    filelike.write("#" + "\t".join(header) + "\n")


def mkdir_p(path: str) -> None:
    """ mkdir -p path

    Args:
        path: The directory to create
    """
    try:
        os.makedirs(os.path.abspath(path), exist_ok=True)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
    if not os.path.isdir(path):
        raise Exception(f"Failed to create directory {path}")


class BGZipFile:
    """ BGZip file helper
    """

    def __init__(self, filename: str, force: bool = False):
        """ Make a subprocess for bgzip

        Args:
            filename: name of the output file
            force: true to overwrite if file exists
        """
        if os.path.exists(filename) and not force:
            raise Exception(f"File {filename} exists, use force=True to overwrite")

        self.write_file = open(filename, "wb")
        zip_pipe = subprocess.Popen(["bgzip", "-f"],
                                    stdin=subprocess.PIPE,
                                    stdout=self.write_file,
                                    stderr=subprocess.PIPE,
                                    shell=True)
        self.zip_pipe = zip_pipe
        self.name = filename

    def close(self) -> None:
        """Close the file handle and subprocess."""
        self.zip_pipe.stdin.flush()
        self.zip_pipe.stdin.close()
        self.zip_pipe.wait()
        self.write_file.flush()
        self.write_file.close()

    def write(self, *args, **kwargs) -> None:
        """Write to the BGZipped file."""
        if isinstance(args[0], str):
            # Convert string to bytes for Python 3
            self.zip_pipe.stdin.write(args[0].encode('utf-8'))
        else:
            self.zip_pipe.stdin.write(*args, **kwargs)

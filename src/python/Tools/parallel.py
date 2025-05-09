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

import argparse
import gc
import logging
import multiprocessing
import os
import pickle
import subprocess
import sys
import tempfile
import traceback
from itertools import islice, repeat

from . import LoggingWriter


# Python 3 compatibility for file handling
def open_file(filename, mode="r"):
    """Helper function to open files in the correct mode for both text and binary."""
    if "b" in mode:
        return open(filename, mode)
    else:
        return open(filename, mode, encoding="utf-8")


POOL = None


def getPool(threads):
    """get / create pool"""
    global POOL
    if POOL:
        return POOL
    elif threads > 1:
        POOL = multiprocessing.Pool(threads)
        return POOL
    else:
        return None


def splitEvery(n, iterable):
    """split iterable into list blocks of size n"""
    if n is None:
        yield list(iterable)
    else:
        i = iter(iterable)
        piece = list(islice(i, n))
        while piece:
            yield piece
            piece = list(islice(i, n))


def unpickleSequentially(plist):
    """Unpickle and concatenate sequentially"""
    data = []
    while plist or data:
        if not data:
            fname = plist.pop(0)
            with open(fname, "rb") as f:  # Use binary mode for pickle
                data = pickle.load(f)
            os.unlink(fname)
        if data:
            yield data.pop(0)


def parMapper(arg):
    try:
        # garbage collect so we can reuse memory
        # when running on very large inputs
        gc.collect()
        return arg[1]["fun"](arg[0], *arg[1]["args"], **arg[1]["kwargs"])
    except Exception as e:
        logging.error("Exception when running %s:" % str(arg[1]['fun']))
        logging.error('-'*60)
        with LoggingWriter(logging.ERROR) as lw:
            traceback.print_exc(file=lw)
        logging.error('-'*60)
    except BaseException as e:
        logging.error("Exception when running %s:" % str(arg[1]['fun']))
        logging.error('-'*60)
        with LoggingWriter(logging.ERROR) as lw:
            traceback.print_exc(file=lw)
        logging.error('-'*60)
    return None


def runParallel(pool, fun, par, *args, **kwargs):
    """run a function in parallel on all elements in par

    :param pool: multiprocessing.Pool or None
    :param fun: a function
    :param par: a list of things to map to (each item is passed as the first argument to fun)
    :param args: more function arguments for fun
    :param kwargs: more function arguments for fun

    """
    if pool:
        result = pool.map(
            parMapper,
            list(zip(par, repeat({"fun": fun, "args": args, "kwargs": kwargs}))),
        )
    else:
        result = []
        for c in par:
            result.append(parMapper((c, {"fun": fun, "args": args, "kwargs": kwargs})))
    return result

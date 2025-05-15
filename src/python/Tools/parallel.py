#!/usr/bin/env python3
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
import sys
import gc
import argparse
import logging
import traceback
import subprocess
import multiprocessing
import pickle
import tempfile
from itertools import islice, repeat
from typing import Any, Callable, Dict, Iterator, List, Optional, Tuple, TypeVar, Union

from . import LoggingWriter

# Define type variables for generics
T = TypeVar('T')
R = TypeVar('R')

# Global pool for reuse
POOL = None


def getPool(threads: int) -> Optional[multiprocessing.Pool]:
    """ Get or create a process pool
    
    Args:
        threads: Number of threads to use
        
    Returns:
        Multiprocessing pool or None for single-threaded execution
    """
    global POOL
    if POOL:
        return POOL
    elif threads > 1:
        POOL = multiprocessing.Pool(threads)
        return POOL
    else:
        return None


def splitEvery(n: Optional[int], iterable: Iterator[T]) -> Iterator[List[T]]:
    """ Split iterable into list blocks of size n 
    
    Args:
        n: Size of each chunk, if None yields entire iterable as one chunk
        iterable: Input iterator to split
        
    Yields:
        Lists of items from the iterable with maximum size n
    """
    if n is None:
        yield list(iterable)
    else:
        i = iter(iterable)
        piece = list(islice(i, n))
        while piece:
            yield piece
            piece = list(islice(i, n))


def unpickleSequentially(plist: List[str]) -> Iterator[Any]:
    """ Unpickle and concatenate sequentially 
    
    Args:
        plist: List of pickle file paths
        
    Yields:
        Unpickled objects from the files
    """
    data = []
    plist_copy = list(plist)  # Create a copy to avoid modifying the input list
    
    while plist_copy or data:
        if not data:
            fname = plist_copy.pop(0)
            with open(fname, 'rb') as f:  # Open in binary mode for pickle
                data = pickle.load(f)
            os.unlink(fname)
        if data:
            yield data.pop(0)


def parMapper(arg: Tuple[T, Dict[str, Any]]) -> Optional[R]:
    """Wrapper function for parallel mapping
    
    Args:
        arg: Tuple with (item, function_info)
        
    Returns:
        Result of function call or None if exception occurred
    """
    try:
        # Garbage collect so we can reuse memory
        # when running on very large inputs
        gc.collect()
        return arg[1]['fun'](arg[0], *arg[1]['args'], **arg[1]['kwargs'])
    except Exception as e:
        logging.error(f"Exception when running {str(arg[1]['fun'])}:")
        logging.error('-'*60)
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        logging.error('-'*60)
    except BaseException as e:
        logging.error(f"Exception when running {str(arg[1]['fun'])}:")
        logging.error('-'*60)
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        logging.error('-'*60)
    return None


def runParallel(
    pool: Optional[multiprocessing.Pool], 
    fun: Callable[[T], R], 
    par: List[T], 
    *args: Any, 
    **kwargs: Any
) -> List[Optional[R]]:
    """ Run a function in parallel on all elements in par
    
    Args:
        pool: Multiprocessing.Pool or None for sequential execution
        fun: Function to apply to each item
        par: List of items to process (each item is passed as the first argument to fun)
        args: Additional positional arguments for fun
        kwargs: Additional keyword arguments for fun
    
    Returns:
        List of results in the same order as input items
    """
    # Create a list of (item, function_info) tuples for mapping
    func_info = {"fun": fun, "args": args, "kwargs": kwargs}
    
    if pool:
        # Use izip in Python 2, but in Python 3 zip is already lazy
        result = pool.map(parMapper, zip(par, repeat(func_info)))
    else:
        # Sequential execution
        result = []
        for item in par:
            result.append(parMapper((item, func_info)))
    
    return result

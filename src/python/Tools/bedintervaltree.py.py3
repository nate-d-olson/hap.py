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

import gzip
from bx.intervals.intersection import Interval
from bx.intervals.intersection import IntervalTree
from collections import defaultdict
from typing import List, Union, Callable, Optional, Dict, Any


class BedIntervalTree:
    """Reads in a BED file and converts it to an interval tree for searching"""
    
    def __init__(self):
        self.tree = defaultdict(IntervalTree)
        self.intCount = 0

        def mkzero():
            return int(0)
        self.count_by_label = defaultdict(mkzero)
        self.nt_count_by_label = defaultdict(mkzero)

    def __str__(self) -> str:
        return f"{self.intCount} intervals"

    def __repr__(self) -> str:
        return str(self)

    def _addEntryToTree(self, bedentry: List, label: str) -> None:
        """ Add a BED entry to the tree
        
        Args:
            bedentry: BED entry [chr, start, stop(, optional extra fields)]
            label: the label for the entry
        """
        chrom = bedentry[0]
        start = int(bedentry[1]) + 1
        end = int(bedentry[2]) + 1
        lbl = [label] + bedentry[3:]
        currInt = Interval(start, end, value=lbl, chrom=chrom)
        self.tree[chrom].add_interval(currInt)
        self.count_by_label[label] += 1
        self.nt_count_by_label[label] += end - start
        self.intCount += 1

    def intersect(self, chrom: str, start: int, end: int) -> List[Interval]:
        """ Return all overlapping intervals in chr:[start,end)
        
        Args:
            chrom: Chromosome
            start: start (1-based)
            end: end
            
        Returns:
            List of Interval objects
            
        Note:
            Intervals have a value associated, this value is an array -- the first column will be
            the label, followed by the bed columns
        """
        return self.tree[chrom].find(start, end)

    def countbases(self, chrom: Optional[str] = None, start: int = 0, end: int = 0, 
                  label: Optional[str] = None) -> int:
        """ Return the number of bases covered by intervals in chr:[start,end)
        
        Args:
            chrom: Chromosome
            start: start (1-based)
            end: end
            label: label
            
        Returns:
            Number of bases covered
        """
        if not chrom and label:
            return self.nt_count_by_label[label]
        elif not chrom and not label:
            return sum([self.nt_count_by_label[x] for x in self.nt_count_by_label.keys()])

        total_length = 0
        for x in self.tree[chrom].find(start, end):
            if not label or x.value == label:
                total_length += x.end - x.start
        return total_length

    def count(self, label: Optional[str] = None) -> int:
        """ Return number of records per label
        
        Args:
            label: string label
            
        Returns:
            Number of intervals which have the given label
        """
        if not label:
            return self.intCount
        else:
            return self.count_by_label[label]

    def addFromBed(self, bed_file: str, label: Union[str, Callable] = "fp", fixchr: bool = False) -> None:
        """ Add all intervals from a bed file, attaching a given label
        
        Args:
            bed_file: Bed File path
            label: Either a string label or a function to work on the bed columns
            fixchr: Fix chr prefix for contig names

        Notes:
            label can be something like this:

                def labeller(entry):
                    # return column 4
                    return entry[3]

            When None is passed, we'll use the first value in the bed column that comes up
            -> chr start end <this one>
        """
        if bed_file.endswith(".gz"):
            with gzip.open(bed_file, 'rt') as bed:
                self._process_bed_file(bed, label, fixchr)
        else:
            with open(bed_file, 'r') as bed:
                self._process_bed_file(bed, label, fixchr)
    
    def _process_bed_file(self, bed_file, label: Union[str, Callable], fixchr: bool) -> None:
        """Process the bed file content
        
        Args:
            bed_file: Open file object for reading
            label: Label or labeler function
            fixchr: Fix chr prefix for contig names
        """
        if callable(label):
            labeller = label
        elif label:
            labeller = lambda _: label
        else:
            labeller = lambda e: ",".join(map(str, e[3:]))

        for entry in bed_file:
            # split comma / semicolon entries
            fields = entry.replace(";", "\t").replace(",", "\t").strip().split("\t")
            if fixchr:
                fields[0] = "chr" + fields[0].replace("chr", "")

            # Apply the labeller function
            entry_label = labeller(fields)
            self._addEntryToTree(fields, entry_label)

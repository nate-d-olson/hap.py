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
from collections import defaultdict

from bx.intervals.intersection import Interval, IntervalTree


# Python 3 compatibility for file handling
def open_file(filename, mode="r"):
    """Helper function to open files in the correct mode for both text and binary."""
    if "b" in mode:
        return open(filename, mode)
    else:
        return open(filename, mode, encoding="utf-8")


class BedIntervalTree(object):
    def __init__(self):
        """reads in a BED file and converts it to an interval tree for searching"""
        self.tree = defaultdict(IntervalTree)
        self.intCount = 0

        def mkzero():
            return int(0)

        self.count_by_label = defaultdict(mkzero)
        self.nt_count_by_label = defaultdict(mkzero)

    def __str__(self):
        return str(self.intCount) + " intervals"

    def __repr__(self):
        return str(self)

    def _addEntryToTree(self, bedentry, label):
        """Add a BED entry to the tree
        :param bedentry: BED entry [chr, start, stop(, optional extra fields)]
        :type bedentry: list of int and string
        :param label: the label for the entry
        :type label: str
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

    def intersect(self, chrom, start, end):
        """Return all overlapping intervals in chr:[start,end)
        :param chrom: Chromosome
        :param start: start (1-based)
        :param end: end
        :rtype: list of Interval

        Intervals have a value associated, this value is an array -- the first column will be
        the label, followed by the bed columns

        """
        return self.tree[chrom].find(start, end)

    def countbases(self, chrom=None, start=0, end=0, label=None):
        """Return the number of bases covered by intervals in chr:[start,end)
        :param chrom: Chromosome
        :param start: start (1-based)
        :param end: end
        :param label: label
        :rtype: int

        Intervals have a value associated, this value is an array -- the first column will be
        the label, followed by the bed columns

        """
        if not chrom and label:
            return self.nt_count_by_label[label]
        elif not chrom and not label:
            return sum(
                [self.nt_count_by_label[x] for x in list(self.nt_count_by_label.keys())]
            )

        total_length = 0
        for x in self.tree[chrom].find(start, end):
            if not label or x.value == label:
                total_length += x.end - x.start
        return total_length

    def count(self, label=None):
        """Return number of records per label
        :param label: string label
        :return: number of intervals which have the given label
        """
        if not label:
            return self.intCount
        else:
            return self.count_by_label[label]

    def addFromBed(self, bed_file, label="fp", fixchr=False):
        """Add all intervals from a bed file, attaching a given label
        :param bed_file: Bed File
        :param label: either a string label or a function to work on the bed columns
        :param fixchr: fix chr prefix for contig names

        label can be something like this:

            def labeller(entry):
                # return column 4
                return entry[3]

        When None is passed, we'll use the first value in the bed column that comes up
        -> chr start end <this one>

        """
        if bed_file.endswith(".gz"):
            with gzip.open(
                bed_file, "rt", encoding="utf-8"
            ) as bed:  # Use text mode with encoding
                self._process_bed_file(bed, label, fixchr)
        else:
            with open_file(bed_file) as bed:
                self._process_bed_file(bed, label, fixchr)

    def _process_bed_file(self, bed, label, fixchr):
        """Process the bed file and add entries to the tree
        :param bed: Bed file object
        :param label: label or labeller function
        :param fixchr: fix chr prefix for contig names
        """
        if hasattr(label, "__call__"):
            labeller = label
        elif label:
            labeller = lambda _: label
        else:
            labeller = lambda e: ",".join(map(str, e[3:]))

        for entry in bed:
            # split comma / semicolon entries
            fields = entry.replace(";", "\t").replace(",", "\t").strip().split("\t")
            if fixchr:
                fields[0] = "chr" + fields[0].replace("chr", "")

            # we've made sure.
            # noinspection PyCallingNonCallable
            label = labeller(fields)

            self._addEntryToTree(fields, label)

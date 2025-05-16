#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt

"""
Module for computing ROC curves from genomic variant data tables.
Provides functionality for precision/recall calculations and
handling of variant classification data.
"""

import abc
import contextlib
import logging
import os
import subprocess
import tempfile
from typing import ClassVar, Dict, Optional, Type

import pandas


def tableROC(
    tbl: pandas.DataFrame,
    label_column: str,
    feature_column: str,
    filter_column: Optional[str] = None,
    filter_name: Optional[str] = None,
    roc_reversed: bool = False,
) -> pandas.DataFrame:
    """Compute ROC table from TP/FP/FN classification table.

    Args:
        tbl: Table with label and feature columns
        label_column: Column name which gives the label (TP/FP/FN)
        feature_column: Column name which gives the feature
        filter_column: Column that contains the filter fields
        filter_name: Column that contains the filter name
        roc_reversed: Reverse ROC behavior

    Returns:
        A pandas.DataFrame with TP/FP/FN/precision/recall columns.
    """
    with (
        tempfile.NamedTemporaryFile(delete=False) as tf1,
        tempfile.NamedTemporaryFile(delete=False) as tf2,
    ):
        try:
            fields = [feature_column, label_column]
            if filter_column:
                fields.append(filter_column)

            tbl[fields].to_csv(tf2.name, sep="\t", index=False)

            cmdline = f"roc -t {label_column} -v {feature_column} --verbose "
            if filter_column:
                cmdline += f" -f {filter_column}"
            if filter_name:
                cmdline += f" -n {filter_name}"
            if roc_reversed:
                cmdline += " -R 1"
            cmdline += f" -o {tf1.name} {tf2.name}"

            logging.info(f"Running {cmdline}")

            subprocess.check_call(cmdline, shell=True)
            try:
                result = pandas.read_table(tf1.name)
            except Exception:
                raise Exception("Cannot parse ROC output.")
            return result
        finally:
            with contextlib.suppress(Exception):
                os.unlink(tf1.name)

            with contextlib.suppress(Exception):
                os.unlink(tf2.name)


class ROC(abc.ABC):
    """ROC calculator base class."""

    classes: ClassVar[Dict[str, Type["ROC"]]] = {}
    features: ClassVar[Dict[str, str]] = {}

    def __init__(self) -> None:
        self.ftable = ""
        self.ftname: str = ""

    @abc.abstractmethod
    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        """Create ROC from feature table.

        Args:
            tbl: The input data table

        Returns:
            DataFrame containing ROC curve data
        """
        pass

    @classmethod
    def make(cls, cname: str) -> "ROC":
        """Factory method to create ROC instances.

        Args:
            cname: Name of the ROC class to instantiate

        Returns:
            An initialized ROC calculator instance
        """
        # Create an instance of the requested ROC class
        c = cls.classes[cname]()
        c.ftname = cls.features[cname]
        return c

    @classmethod
    def register(cls, name: str, feature_name: str) -> callable:
        """Register a new ROC class.

        Args:
            name: Name to register the class under
            feature_name: Name of the feature to use

        Returns:
            Decorator function that registers the class
        """

        def _wrap(subclass):
            cls.classes[name] = subclass
            cls.features[name] = feature_name
            return subclass

        return _wrap


@ROC.register("qual", "QUAL")
class QualROC(ROC):
    """QUAL-based ROC implementation."""

    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        """Create ROC from feature table using QUAL.

        Args:
            tbl: The input data table

        Returns:
            DataFrame containing ROC curve data
        """
        return tableROC(tbl, "type", self.ftname, filter_column="filter")


@ROC.register("vqslod", "VQSLOD")
class VQSLODROC(ROC):
    """VQSLOD-based ROC implementation."""

    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        """Create ROC from feature table using VQSLOD.

        Args:
            tbl: The input data table

        Returns:
            DataFrame containing ROC curve data
        """
        return tableROC(tbl, "type", self.ftname, roc_reversed=True)
        return tableROC(tbl, "type", self.ftname, roc_reversed=True)
        return tableROC(tbl, "type", self.ftname, roc_reversed=True)

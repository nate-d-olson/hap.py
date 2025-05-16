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

"""
Module for computing ROC curves from genomic variant data tables.
Provides functionality for precision/recall calculations and
handling of variant classification data.
"""

import abc
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

    with tempfile.NamedTemporaryFile(delete=False) as tf1, tempfile.NamedTemporaryFile(
        delete=False
    ) as tf2:
        tf1_name = tf1.name
        tf2_name = tf2.name
    try:
        fields = [feature_column, label_column]
        if filter_column:
            fields.append(filter_column)

        tbl[fields].to_csv(tf2_name, sep="\t", index=False)

        cmdline = "roc -t %s -v %s --verbose " % (label_column, feature_column)
        if filter_column:
            cmdline += " -f %s" % filter_column
        if filter_name:
            cmdline += " -n %s" % filter_name
        if roc_reversed:
            cmdline += " -R 1"
        cmdline += " -o %s %s" % (tf1_name, tf2_name)

        logging.info("Running %s" % cmdline)

        subprocess.check_call(cmdline, shell=True)
        try:
            result = pandas.read_csv(tf1_name, sep="\t")
        except Exception:
            raise Exception("Cannot parse ROC output.")
        return result
    finally:
        try:
            os.unlink(tf1_name)
        except Exception:
            pass
        try:
            os.unlink(tf2_name)
        except Exception:
            pass


class ROC(metaclass=abc.ABCMeta):
    """ROC calculator base class"""

    classes: ClassVar[Dict[str, Type["ROC"]]] = {}
    features: ClassVar[Dict[str, str]] = {}

    def __init__(self) -> None:
        self.ftable: str = ""
        self.ftname: str = ""

    @abc.abstractmethod
    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        """Create ROC from feature table

        Args:
            tbl: The input table with variant data

        Returns:
            DataFrame with ROC curve data
        """
        pass

    @classmethod
    def make(cls, cname: str) -> "ROC":
        """Create an instance of a ROC calculator

        Args:
            cname: Name of the ROC calculator to create

        Returns:
            An instance of the requested ROC calculator
        """
        # noinspection PyCallingNonCallable
        c = cls.classes[cname]()
        c.ftname = cls.features[cname]
        return c

    @classmethod
    def register(cls, name: str, ftname: str, cons: Type["ROC"]) -> None:
        """Register a ROC calculator

        Args:
            name: The name of the calculator
            ftname: The features/feature table name (will be accessible in the ftname attribute)
            cons: Class constructor
        """
        cls.classes[name] = cons
        cls.features[name] = ftname

    @classmethod
    def list(cls) -> list:
        """List all registered ROC calculators

        Returns:
            List of all registered ROC calculator names
        """
        return list(cls.classes.keys())


class StrelkaSNVRoc(ROC):
    """ROC calculator for Strelka SNVs"""

    def __init__(self) -> None:
        super().__init__()
        self.ftname: str = ""  # Added to resolve missing attribute error

    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        tbl.loc[tbl["NT"] != "ref", "QSS_NT"] = 0
        return tableROC(tbl, "tag", "QSS_NT", "FILTER", "QSS_ref")


ROC.register("strelka.snv.qss", "hcc.strelka.snv", StrelkaSNVRoc)


class StrelkaSNVVQSRRoc(ROC):
    """ROC calculator for Strelka SNVs (newer versions which use VQSR)"""

    def __init__(self) -> None:
        super().__init__()
        self.ftname: str = ""

    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        tbl.loc[tbl["NT"] != "ref", "VQSR"] = 0
        return tableROC(tbl, "tag", "VQSR", "FILTER", "LowQscore")


ROC.register("strelka.snv.vqsr", "hcc.strelka.snv", StrelkaSNVVQSRRoc)


class StrelkaSNVEVSRoc(ROC):
    """ROC calculator for Strelka SNVs (newer versions where VQSR is called EVS)"""

    def __init__(self) -> None:
        super().__init__()
        self.ftname: str = ""

    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        tbl.loc[tbl["NT"] != "ref", "EVS"] = 0
        return tableROC(tbl, "tag", "EVS", "FILTER", "LowEVS")


ROC.register("strelka.snv", "hcc.strelka.snv", StrelkaSNVEVSRoc)


class StrelkaIndelRoc(ROC):
    """ROC calculator for Strelka Indels"""

    def __init__(self) -> None:
        super().__init__()
        self.ftname: str = ""

    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        # fix QSI for NT != ref
        tbl.loc[tbl["NT"] != "ref", "QSI_NT"] = 0
        return tableROC(tbl, "tag", "QSI_NT", "FILTER", "QSI_ref")


ROC.register("strelka.indel", "hcc.strelka.indel", StrelkaIndelRoc)


class StrelkaIndelEVSRoc(ROC):
    """ROC calculator for Strelka Indels"""

    def __init__(self) -> None:
        super().__init__()
        self.ftname: str = ""

    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        # fix QSI for NT != ref
        return tableROC(tbl, "tag", "EVS", "FILTER", "LowEVS")


ROC.register("strelka.indel.evs", "hcc.strelka.indel", StrelkaIndelEVSRoc)


class Varscan2SNVRoc(ROC):
    """ROC calculator for Varscan2 SNVs"""

    def __init__(self) -> None:
        super().__init__()
        self.ftname: str = ""

    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        return tableROC(tbl, "tag", "SSC")


ROC.register("varscan2.snv", "hcc.varscan2.snv", Varscan2SNVRoc)


class Varscan2IndelRoc(ROC):
    """ROC calculator for Varscan2 Indels"""

    def __init__(self) -> None:
        super().__init__()
        self.ftname: str = ""

    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        return tableROC(tbl, "tag", "SSC")


ROC.register("varscan2.indel", "hcc.varscan2.indel", Varscan2IndelRoc)


class MutectSNVRoc(ROC):
    """ROC calculator for MuTect SNVs"""

    def __init__(self) -> None:
        super().__init__()
        self.ftname: str = ""

    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        return tableROC(tbl, "tag", "TLOD", "FILTER", "t_lod_fstar")


ROC.register("mutect.snv", "hcc.mutect.snv", MutectSNVRoc)


class MutectIndelRoc(ROC):
    """ROC calculator for MuTect Indels"""

    def __init__(self) -> None:
        super().__init__()
        self.ftname: str = ""

    def from_table(self, tbl: pandas.DataFrame) -> pandas.DataFrame:
        return tableROC(tbl, "tag", "TLOD", "FILTER", "t_lod_fstar")


ROC.register("mutect.indel", "hcc.mutect.indel", MutectIndelRoc)

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
#
# 9/9/2014
#
# Diploid VCF File Comparison
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import sys
import os
import argparse
import logging
import traceback
import multiprocessing
import pandas
import tempfile
import json
from typing import Dict, List, Optional, Any, Tuple, Union

scriptDir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(scriptDir, "..", "lib", "python3")))

import Tools


def parseArgs():
    """Parse the command line arguments"""
    parser = argparse.ArgumentParser("Haplotype Comparison ROC")

    parser.add_argument("roc", help="ROC data file")

    parser.add_argument("-o", "--reports-prefix", dest="prefix",
                      help="Prefix for report output files", default=None)

    parser.add_argument("-V", "--verbose", dest="verbose",
                      help="Increase verbosity of output", action="count", default=0)

    parser.add_argument("--logfile", dest="logfile",
                      help="Write logging information into file rather than to stderr",
                      default=None)

    parser.add_argument("-l", "--location", dest="location",
                      help="Comma-separated list of locations to restrict to.", default=None)

    parser.add_argument("-t", "--type", dest="type",
                      help="Variant types and regions to include.", default=None)

    parser.add_argument("-f", "--feature-table", nargs="+", dest="feature_table",
                      help="Features to annotate variants with and stratify on.",
                      default=[])

    parser.add_argument("--feature-filter", dest="feature_filter",
                      help="Feature pattern to limit feature selection to.",
                      default=None)

    parser.add_argument("--embed-subplots", dest="embed_subplots", help="Embed subplots in a big figure.",
                      default=False, action="store_true")

    parser.add_argument("--roc-filter", dest="roc_filter",
                      help="Comma-separated list of INFO fields to use.",
                      default=None)

    parser.add_argument("--roc-regions", dest="roc_regions",
                      help="Comma-separated list of regions to plot.",
                      default="*")

    parser.add_argument("--roc-kind", dest="roc_kind",
                      help="ROC kind: count-all (count all records), count-nonref, roc",
                      default="count-all")

    parser.add_argument("--stratification", dest="stratification",
                      help="List of stratifications, comma-separated.")

    parser.add_argument("--format", dest="format", help="Output format -- py plots into separate files",
                      default="html",
                      choices=["html", "png", "pdf", "py", "jpg", "all"])

    parser.add_argument("--width", dest="width", help="Width of plot in inches.", 
                      default=None, type=int)

    parser.add_argument("--height", dest="height", help="Height of plot in inches.", 
                      default=None, type=int)

    parser.add_argument("--nc", "--no-count-filtered", dest="count_filtered",
                      help="Don't count filtered variants in ROC.",
                      action="store_false", default=True)

    parser.add_argument("--chm", "--count-homref-matches", dest="homref_matches",
                      help="Also count homref matches.",
                      action="store_true", default=False)

    parser.add_argument("-L", "--legacy-output", dest="legacy",
                      help="Make legacy output.",
                      action="store_true", default=False)

    parser.add_argument("--title", dest="title", help="Plot title", default="ROC")
    parser.add_argument("--qq-max", dest="qq_max", help="Max value for QQ plots.",
                      type=float, default=None)
    parser.add_argument("--roc-lines", dest="roc_lines", help="Plot lines from spec file.",
                      default=None)

    parser.add_argument("--no-names", dest="no_names",
                      help="Don't write names into plots (when plotting multiple curves only).",
                      action="store_true", default=False)

    parser.add_argument("--logx", dest="logx",
                      help="Log scale for x axis",
                      action="store_true", default=False)

    parser.add_argument("--logy", dest="logy",
                      help="Log scale for y axis",
                      action="store_true", default=False)

    parser.add_argument("--qqnorm", dest="qqnorm",
                      help="Normalise QQ plot",
                      action="store_true", default=False)

    parser.add_argument("--xmax", dest="xmax", help="Maximum value for x axis.", 
                      type=float, default=None)
    
    parser.add_argument("--ymin", dest="ymin", help="Minimum value for y axis.", 
                      type=float, default=None)

    parser.add_argument("--ymax", dest="ymax", help="Maximum value for y axis.", 
                      type=float, default=None)

    parser.add_argument("--roc-params", dest="roc_params", nargs="+",
                      help="All remaining args are parsed as key=value strings and passed to hapyroc.",
                      default=[])

    result = parser.parse_args()

    # if no prefix was given, derive from the ROC file name
    if result.prefix is None:
        if result.roc.lower().endswith(".csv"):
            result.prefix = result.roc[:-4]
        elif result.roc.lower().endswith(".csv.gz"):
            result.prefix = result.roc[:-7]
        elif result.roc.lower().endswith(".h5") or result.roc.lower().endswith(".hdf"):
            result.prefix = result.roc[:-3]
        else:
            result.prefix = result.roc
    
    return result


def main():
    """Main function"""
    try:
        args = parseArgs()

        # how verbose?
        if args.verbose == 0:
            loglevel = logging.WARN
        elif args.verbose == 1:
            loglevel = logging.INFO
        else:
            loglevel = logging.DEBUG
        
        if args.logfile:
            logging.basicConfig(filename=args.logfile, level=loglevel)
        else:
            logging.basicConfig(level=loglevel)

        logging.getLogger("requests").setLevel(logging.WARNING)

        import Haplo.happyroc

        # Attempting to create a double instance of the plot can fail sometimes since 
        # matplotlib isn't really re-entrant
        # plotting_ok = True
        roc_args = {}

        # Set format
        if args.format in ["pdf", "png", "jpg"]:
            roc_args["output_format"] = args.format
        elif args.format == "all":
            # plot into all formats
            roc_args["output_format"] = ["pdf", "png", "jpg"]

        if args.width:
            roc_args["width"] = args.width
        
        if args.height:
            roc_args["height"] = args.height

        if args.logx:
            roc_args["logx"] = True
            
        if args.logy:
            roc_args["logy"] = True

        if args.qqnorm:
            roc_args["qqnorm"] = True

        if args.xmax is not None:
            roc_args["xmax"] = args.xmax
            
        if args.ymin is not None:
            roc_args["ymin"] = args.ymin

        if args.ymax is not None:
            roc_args["ymax"] = args.ymax

        if args.qq_max:
            roc_args["qq_max"] = args.qq_max

        if args.no_names:
            roc_args["no_names"] = True

        # RocPlot will do Agg plotting from the main thread 
        # this is safer than trying to plot using the main UI thread
        roc_args["agg"] = True
        if args.roc_lines:
            roc_args["line_specs"] = args.roc_lines
        
        for param in args.roc_params:
            k, v = param.split("=", 1)
            try:
                v = json.loads(v)
            except ValueError:
                pass
            roc_args[k] = v
            
        # determine files to output to
        roc_file = args.roc
        roc_filter = "QUAL"
        
        if args.roc_filter:
            roc_filter = args.roc_filter
            
        roc_regions = args.roc_regions.split(",")

        roc_type = "all"
        if args.type:
            roc_type = args.type.split(",")

        # if we have an HDF, try to use the locations from there
        locations = []
        if args.location:
            locations = args.location.split(",")
            if len(locations) == 1 and "*" == locations[0]:
                locations = ["*"]

        if args.stratification:
            strat = args.stratification.split(",")
            roc_args["stratification"] = strat
        
        if args.feature_table:
            ft = args.feature_table
            regions_bedfiles = []
            for f in ft:
                if os.path.exists(f):
                    regions_bedfiles.append(f)
                else:
                    raise Exception("Cannot find %s" % f)
                    
            if regions_bedfiles:
                import Haplo.happyroc
                fp = tempfile.NamedTemporaryFile(delete=False, suffix=".bed")
                try:
                    fp.close()
                    Haplo.happyroc.mergeFeatures(regions_bedfiles, fp.name)
                    roc_args["features_bedfile"] = fp.name
                except:
                    os.unlink(fp.name)
                    raise

        if args.feature_filter:
            roc_args["feature_filter"] = args.feature_filter

        if args.format == "py":
            # Use Python plots instead of HTML
            Haplo.happyroc.makeRocPlots(roc_file,
                                      roc_filter,
                                      args.prefix,
                                      args.title,
                                      roc_regions,
                                      locations,
                                      roc_type,
                                      args.count_filtered,
                                      args.homref_matches,
                                      roc_args,
                                      args.roc_kind)
        else:
            # default HTML plot
            Haplo.happyroc.makeRocHTML(roc_file,
                                     roc_filter,
                                     args.prefix,
                                     args.title,
                                     roc_regions,
                                     locations,
                                     roc_type,
                                     args.count_filtered,
                                     args.homref_matches,
                                     roc_args,
                                     args.roc_kind,
                                     args.embed_subplots)
        
        if args.legacy:
            # emit legacy style outputs
            Haplo.happyroc.makeLegacyOutput(roc_file, args.prefix)

        return 0
    except Exception as e:
        print(str(e))
        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

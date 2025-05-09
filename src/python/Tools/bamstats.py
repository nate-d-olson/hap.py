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

import logging

import numpy as np
import pandas
import pysam


# Python 3 compatibility for file handling
def open_file(filename, mode="r"):
    """Helper function to open files in the correct mode for both text and binary."""
    if "b" in mode:
        return open(filename, mode)
    else:
        return open(filename, mode, encoding="utf-8")


def bamStats(bamfile):
    """ Extract average depths + idxstats data from BAM file, return data frame """
    try:
        with pysam.AlignmentFile(bamfile, "rb") as samfile:
            istats = pysam.idxstats(bamfile)
            
            # Handle bytes vs. string in Python 3
            if isinstance(istats, bytes):
                istats = istats.decode('utf-8')
                
            result = []
            for x in istats.split("\n"):
                xs = x.replace("\n", "").split("\t")
                if len(xs) < 4:
                    logging.warn("Ignoring invalid stats line: %s" % x)
                    continue
                rec = {
                    "CHROM": xs[0],
                    "NT": int(xs[1]),
                    "MAPPED": int(xs[2]),
                    "UNMAPPED": int(xs[3]),
                    "READLEN": 0,
                    "COVERAGE": 0.0,
                }

                count = 0
                rls = 0.0
                try:
                    for read in samfile.fetch(xs[0]):
                        # In Python 3, pysam might return read.rlen or read.query_length
                        read_length = getattr(read, 'query_length', None)
                        if read_length is None:
                            read_length = getattr(read, 'rlen', 0)
                        
                        rls += float(read_length)
                        count += 1
                        if count > 10000:
                            break

                    rls /= count if count > 0 else 1.0
                    rec["READLEN"] = rls
                    rec["COVERAGE"] = float(rec["MAPPED"] * rec["READLEN"])/float(rec["NT"]) if rec["NT"] > 0 else 0.0
                except Exception as e:
                    logging.warn(f"Error processing {xs[0]}: {str(e)}")
                    pass
                result.append(rec)

            if result:
                result = pandas.DataFrame(result, columns=["CHROM", "NT", "MAPPED", "UNMAPPED", "READLEN", "COVERAGE"])
                min_readlen = np.min(result[result["READLEN"] > 0]["READLEN"])
                max_readlen = np.max(result[result["READLEN"] > 0]["READLEN"])
                if min_readlen != max_readlen:
                    logging.warn("Read lengths differ within the same BAM file: %s" % str(result["READLEN"].unique()))

                agg_result = [{
                    "CHROM": "TOTAL",
                    "NT": np.sum(result["NT"]),
                    "MAPPED": np.sum(result["MAPPED"]),
                    "UNMAPPED": np.sum(result["UNMAPPED"]),
                    "READLEN": (max_readlen + min_readlen) / 2.0,
                }]
                agg_result[-1]["COVERAGE"] = np.sum(result["MAPPED"].multiply(result["READLEN"])) / np.sum(result["NT"])

                auto_result = result[result["CHROM"].str.match(r"^(?:chr)?[0-9]+$")]
                if not auto_result.empty:
                    auto_readlen_pos = auto_result[auto_result["READLEN"] > 0]
                    if not auto_readlen_pos.empty:
                        min_readlen = np.min(auto_readlen_pos["READLEN"])
                        max_readlen = np.max(auto_readlen_pos["READLEN"])
                    else:
                        min_readlen = max_readlen = 0
                        
                    agg_result.append({
                        "CHROM": "AUTOSOME",
                        "NT": np.sum(auto_result["NT"]),
                        "MAPPED": np.sum(auto_result["MAPPED"]),
                        "UNMAPPED": np.sum(auto_result["UNMAPPED"]),
                        "READLEN": (max_readlen + min_readlen) / 2.0,
                    })
                    nt_sum = np.sum(auto_result["NT"])
                    if nt_sum > 0:
                        agg_result[-1]["COVERAGE"] = np.sum(auto_result["MAPPED"].multiply(auto_result["READLEN"])) / nt_sum
                    else:
                        agg_result[-1]["COVERAGE"] = 0.0

                return pandas.concat([result, pandas.DataFrame(agg_result, columns=["CHROM", "NT", "MAPPED", "UNMAPPED", "READLEN", "COVERAGE"])]).set_index(["CHROM"])
            else:
                return pandas.DataFrame(columns=["CHROM", "NT", "MAPPED", "UNMAPPED", "READLEN", "COVERAGE"])
    except Exception as e:
        logging.error(f"Failed to process BAM file {bamfile}: {str(e)}")
        return pandas.DataFrame(columns=["CHROM", "NT", "MAPPED", "UNMAPPED", "READLEN", "COVERAGE"])

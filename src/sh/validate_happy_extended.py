#!/usr/bin/env python

# Compare som.py stats.csv and hap.py extended.csv files

import sys
import argparse
import csv
import pprint as pp
import re
import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sompy-stats',
                        help='Path to som.py stats.csv file')
    parser.add_argument('--happy-extended',
                        help='Path to hap.py extended.csv file')
    parser.add_argument('--tolerance',
                        help='Tolerance (def 0.001)',
                        type=float,
                        default=0.001)
    args = parser.parse_args()
    return args

def eval_equal(metric_name, count_a, count_b, tol=0.001):
    try:
        count_a = int(count_a)
        count_b = int(count_b)
    except ValueError:
        try:
            count_a = float(count_a)
            count_b = float(count_b)
            if abs(count_a-count_b) < tol:
                return True
            logging.info("%s: abs(%s-%s)=%s (FAIL)" % (metric_name, count_a, count_b, abs(count_a-count_b)))
            return False
        except ValueError:
            pass
    # Test equality
    return count_a == count_b

def parse_sompy_stats(path):
    sompy_stats = csv.DictReader(open(path))
    result = dict()
    for s in sompy_stats:
        subset = s["type"]
        if re.match("indels.", subset):
            # store results per af_bin
            m = re.findall("indels.(\d+\.\d+)-(\d+\.\d+)", subset)[0]
            af_low = m[0][:4]
            af_high = m[1][:4]
            af_bin = "[%s,%s]" % (af_low, af_high) if af_low == "1.00" else "[%s,%s)" % (af_low, af_high)
            result["INDEL." + af_bin] = s

        if re.match("SNVs.", subset):
            # store results per af_bin
            m = re.findall("SNVs.(\d+\.\d+)-(\d+\.\d+)", subset)[0]
            af_low = m[0][:4]
            af_high = m[1][:4]
            af_bin = "[%s,%s]" % (af_low, af_high) if af_low == "1.00" else "[%s,%s)" % (af_low, af_high)
            result["SNP." + af_bin] = s

    return result

if __name__ == '__main__':
    args = parse_args()

    sompy_stats = parse_sompy_stats(path=args.sompy_stats)
    happy_extended = csv.DictReader(open(args.happy_extended))

    outcomes = dict(ALL=set(), PASS=set())
    for h in happy_extended:
        k = h["Type"] + "." + h["Subset"]
        try:
            s = sompy_stats[k]
        except KeyError:
            s = {"total.truth": 0, "tp": 0, "fn": 0, "total.query": 0, "fp": 0, "unk": 0}

        count_fields = [
            ("total.truth", "TRUTH.TOTAL"),
            ("tp", "TRUTH.TP"),
            ("fn", "TRUTH.FN"),
            ("total.query", "QUERY.TOTAL"),
            ("fp", "QUERY.FP"),
            ("unk", "QUERY.UNK"),
        ]

        if h["Filter"] == "ALL" and "ALL" not in k:
            for s_field, h_field in count_fields:
                outcomes["ALL"].add(eval_equal("%s: %s" % (k, s_field), s[s_field], h[h_field], tol=args.tolerance))

        if h["Filter"] == "PASS" and "ALL" not in k:
            for s_field, h_field in count_fields:
                outcomes["PASS"].add(eval_equal("%s: %s" % (k, s_field), s[s_field], h[h_field], tol=args.tolerance))

    logging.info("ALL: %s tests (%s failures)" % (len(outcomes["ALL"]), len([x for x in outcomes["ALL"] if not x])))
    logging.info("PASS: %s tests (%s failures)" % (len(outcomes["PASS"]), len([x for x in outcomes["PASS"] if not x])))

    retcode = 0
    if False in outcomes["ALL"] or False in outcomes["PASS"]:
        retcode = 1

    sys.exit(retcode)

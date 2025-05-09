#!/usr/bin/env python

# Compare som.py stats.csv and hap.py summary.csv files

import sys
import argparse
import csv
import pprint as pp
import logging

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s", level=logging.INFO
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sompy-stats", help="Path to som.py stats.csv file")
    parser.add_argument("--happy-summary", help="Path to hap.py summary.csv file")
    parser.add_argument(
        "--tolerance", help="Tolerance (def 0.001)", type=float, default=0.001
    )
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
            if abs(count_a - count_b) < tol:
                return True
            logging.info(
                "%s: abs(%s-%s)=%s (FAIL)"
                % (metric_name, count_a, count_b, abs(count_a - count_b))
            )
            return False
        except ValueError:
            pass
    # Test equality
    return count_a == count_b


if __name__ == "__main__":
    args = parse_args()
    sompy_stats = csv.DictReader(open(args.sompy_stats))
    happy_summary = csv.DictReader(open(args.happy_summary))

    # compare first row of som.py stats to PASS/ALL rows in hap.py summary
    s = next(sompy_stats)

    if s["type"] == "SNVs":
        vtype = "SNP"
    elif s["type"] == "indels":
        vtype = "INDEL"
    else:
        logging.error("Unrecognised variant type in som.py stats: %s" % s["type"])
        sys.exit(1)

    logging.info("Comparing %s counts..." % vtype)

    # get corresponding row(s) in hap.py summary
    h_rows = {"ALL": dict(), "PASS": dict()}
    for h in happy_summary:
        if h["Type"] == vtype and h["Filter"] in list(h_rows.keys()):
            if len(h_rows[h["Filter"]]) == 0:
                h_rows[h["Filter"]] = h

    outcomes = set()
    count_fields = [
        ("total.truth", "TRUTH.TOTAL"),
        ("tp", "TRUTH.TP"),
        ("fn", "TRUTH.FN"),
        ("total.query", "QUERY.TOTAL"),
        ("fp", "QUERY.FP"),
        ("unk", "QUERY.UNK"),
    ]
    for h_filter, h_row in list(h_rows.items()):
        logging.info("Filter=%s" % h_filter)
        for s_field, h_field in count_fields:
            outcomes.add(
                eval_equal(
                    "%s: %s" % (vtype, s_field),
                    s[s_field],
                    h_row[h_field],
                    tol=args.tolerance,
                )
            )

    logging.info(
        "%s tests (%s failures)" % (len(outcomes), len([x for x in outcomes if not x]))
    )

    if False in outcomes:
        sys.exit(1)
    sys.exit(0)

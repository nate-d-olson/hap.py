#!/usr/bin/env python3

# Compare two sompy feature tables

import sys
import pandas as pd
import numpy as np

eps = 1e-6


def csvread(filename):
    f = pd.read_csv(filename, na_values="")
    return f


def main():
    data1 = csvread(sys.argv[1])
    data2 = csvread(sys.argv[2])

    if len(data1) != len(data2):
        raise Exception("Table lengths differ: %i vs %i" % (len(data1), len(data2)))

    match = 0
    for i in range(0, len(data1)):
        r1 = data1.iloc[i]
        r2 = data2.iloc[i]

        # use exact string lookup & ignore NaN columns
        if (
            r1["CHROM"] == r2["CHROM"]
            and r1["POS"] == r2["POS"]
            and r1["tag"] == r2["tag"]
        ):
            match += 1
        else:
            print("ERROR at %i:" % i)
            print(r1.to_string())
            print(r2.to_string())
            print("")
            if match < 10:
                continue
            else:
                raise Exception("Records don't match.")

    print("All feature table rows match (%i)" % match)

    # match columns
    passed_columns = False
    try:
        float_cols = set([s for s in list(data1) if not s == s.upper()])

        match = 0
        failed = 0
        for i in range(0, len(data1)):
            r1 = data1.iloc[i]
            r2 = data2.iloc[i]
            for f in float_cols:
                if f in list(r1) and f in list(r2):
                    if np.isnan(r1[f]) and np.isnan(r2[f]):
                        match += 1
                    elif np.isnan(r1[f]) or np.isnan(r2[f]):
                        print(
                            "Value mismatch for %s: %s != %s"
                            % (f, str(r1[f]), str(r2[f]))
                        )
                    elif abs(r1[f] - r2[f]) < eps:
                        match += 1
                    else:
                        print("Value mismatch for %s: %f != %f" % (f, r1[f], r2[f]))
                        failed += 1

        if failed > 0:
            raise Exception(
                "Values in tables don't match (%i failures out of %i comparisons)"
                % (failed, match + failed)
            )
        passed_columns = True
        print("All values match (%i comparisons)" % match)
    finally:
        if not passed_columns:
            raise Exception("Columns don't match")


if __name__ == "__main__":
    main()

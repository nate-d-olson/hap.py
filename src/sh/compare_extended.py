#!/usr/bin/env python3

# Compare two extended csv files

import sys
import csv


def csvread(filename):
    f = open(filename)
    freader = csv.DictReader(f)
    data = {}
    label_columns = [
        "Type",
        "Subtype",
        "Subset",
        "Filter",
        "Genotype",
        "QQ.Field",
        "QQ",
    ]

    for l in freader:
        record = dict(l)
        label = "_".join([record[k] for k in label_columns])
        if "hap.py" in label:
            # ignore version line
            continue

        for k in record.keys():
            try:
                record[k] = int(record[k])
            except:
                try:
                    record[k] = float(record[k])
                except:
                    if record[k] == "." or record[k] == "":
                        record[k] = float("NaN")
        data[label] = record

    return data


def main():
    data1 = csvread(sys.argv[1])
    data2 = csvread(sys.argv[2])

    different_metrics = []

    keys1 = set(data1.keys())
    keys2 = set(data2.keys())

    if keys1 != keys2:
        print("Different keys:", file=sys.stderr)
        print("A - B = %s" % str(keys1 - keys2), file=sys.stderr)
        print("B - A = %s" % str(keys2 - keys1), file=sys.stderr)
        print("Only comparing common keys.", file=sys.stderr)

    for key in list(keys1.intersection(keys2)):
        r1_keys = set(data1[key].keys())
        r2_keys = set(data2[key].keys())

        if r1_keys != r2_keys:
            print("Different keys for %s:" % key, file=sys.stderr)
            print("A - B = %s" % str(r1_keys - r2_keys), file=sys.stderr)
            print("B - A = %s" % str(r2_keys - r1_keys), file=sys.stderr)
            print("Only comparing common fields.", file=sys.stderr)

        for field in list(r1_keys.intersection(r2_keys)):
            if field in [
                "Type",
                "Subtype",
                "Subset",
                "Filter",
                "Genotype",
                "QQ.Field",
                "QQ",
            ]:
                if data1[key][field] != data2[key][field]:
                    different_metrics.append(
                        (key, field, str(data1[key][field]), str(data2[key][field]))
                    )
            else:
                try:
                    a = float(data1[key][field])
                    b = float(data2[key][field])

                    if ("%.3g" % a) != ("%.3g" % b):
                        different_metrics.append(
                            (key, field, ("%.3g" % a), ("%.3g" % b), b - a)
                        )
                except:
                    if str(data1[key][field]) != str(data2[key][field]):
                        different_metrics.append(
                            (key, field, str(data1[key][field]), str(data2[key][field]))
                        )

    if different_metrics:
        print("ERROR -- Metric differences detected:", file=sys.stderr)
        print("-------------------------------------\n", file=sys.stderr)
        for m in different_metrics:
            if len(m) > 4:
                print("%s / %s: %s != %s difference: %f" % m, file=sys.stderr)
            else:
                print("%s / %s: %s != %s" % m, file=sys.stderr)
        print("-------------------------------------", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

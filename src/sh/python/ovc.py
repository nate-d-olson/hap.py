#!/usr/bin/env python3
"""Overlap check for blocksplit output intervals."""
import sys


def main():
    if len(sys.argv) != 2:
        print("Usage: ovc.py <intervals.bed>")
        sys.exit(1)
    path = sys.argv[1]
    intervals = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], parts[1], parts[2]
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            intervals.setdefault(chrom, []).append((start_i, end_i))
    # Check for overlaps within each chrom
    for chrom, ivs in intervals.items():
        ivs.sort()
        prev_end = -1
        for st, en in ivs:
            if st < prev_end:
                print(f"Overlap in {chrom}: start {st} < prev_end {prev_end}")
                sys.exit(1)
            prev_end = en
    sys.exit(0)


if __name__ == "__main__":
    main()

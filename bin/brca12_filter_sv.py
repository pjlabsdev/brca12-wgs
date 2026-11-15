#!/usr/bin/env python3
"""Restrict Delly calls to the SV search windows, as Manta does with --callRegions. Research use only.

Keeps a record when its position lies in a window (genes +/- sv_padding) and, for DEL/DUP/INV, also its
END lies in a window and the event spans at least --min-length bases (shorter events are small variants
of HaplotypeCaller/DeepVariant). BND records keep their partner wherever it is.

usage: brca12_filter_sv.py --windows sv.bed [--min-length 50] < in.vcf > out.vcf"""
import argparse
import re
import sys


def load(path):
    windows = {}
    with open(path) as handle:
        for line in handle:
            f = line.split("\t")
            if len(f) >= 3 and not line.startswith("#"):
                windows.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return windows


def inside(windows, chrom, pos):
    return any(s < pos <= e for s, e in windows.get(chrom, ()))


def keep(line, windows, min_length):
    f = line.rstrip("\n").split("\t")
    chrom, pos, info = f[0], int(f[1]), f[7]
    if not inside(windows, chrom, pos):
        return False
    svtype = re.search(r"(?:^|;)SVTYPE=([^;]+)", info)
    svtype = svtype.group(1) if svtype else ""
    if svtype == "BND":
        return True
    end = re.search(r"(?:^|;)END=(\d+)", info)
    end = int(end.group(1)) if end else pos
    return inside(windows, chrom, end) and end - pos >= min_length


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--windows", required=True)
    p.add_argument("--min-length", type=int, default=50)
    args = p.parse_args(argv)
    windows = load(args.windows)
    for line in sys.stdin:
        if line.startswith("#") or keep(line, windows, args.min_length):
            sys.stdout.write(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())

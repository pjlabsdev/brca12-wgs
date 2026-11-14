#!/usr/bin/env python3
"""Exon copy number against a reference cohort (see bin/brca12lib/exon_cnv.py). Research use only.

  brca12_exon_cnv.py --panel PANEL.json --write-bed regions.bed
      BED of the intervals whose depth is needed (exons +/- pad, control windows)
  samtools depth -a -Q 20 -b regions.bed BAM | brca12_exon_cnv.py --panel PANEL.json --out exon_cnv.json
      exon table and copy-number segments
A panel path whose file name starts with NO_ (or --disabled REASON) writes {"enabled": false}."""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from brca12lib import exon_cnv  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel", required=True)
    p.add_argument("--write-bed")
    p.add_argument("--depth", default="-", help="samtools depth output (default: stdin)")
    p.add_argument("--out")
    p.add_argument("--disabled", help="write a disabled result with this reason")
    args = p.parse_args(argv)
    if args.disabled or os.path.basename(args.panel).startswith("NO_"):
        if args.write_bed:
            open(args.write_bed, "w").close()
        if args.out:
            reason = args.disabled or "no reference panel for this assay"
            json.dump({"enabled": False, "reason": reason, "calls": [], "exons": []}, open(args.out, "w"), indent=1)
        return 0
    panel = exon_cnv.load_panel(args.panel)
    rows = exon_cnv.regions(panel)
    if args.write_bed:
        with open(args.write_bed, "w") as out:
            for chrom, start, end, name in rows:
                out.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    if args.out:
        handle = sys.stdin if args.depth == "-" else open(args.depth)
        depths = exon_cnv.mean_depths(handle, rows)
        result = exon_cnv.call(panel, depths)
        # depths and panel travel with the result so that the report can add allele balance from the
        # small-variant calls (brca12_report.py)
        result["depths"] = depths
        result["panelData"] = panel
        json.dump(result, open(args.out, "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())

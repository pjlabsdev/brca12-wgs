#!/usr/bin/env python3
"""Excel export (XLSX) of one brca12-wgs result.

    python3 bin/brca12_xlsx.py brca12_results.json brca12_report.xlsx

Sheets: overview with the research-use-only note, small variants (GRCh38 and GRCh37 positions,
HGVS, BIC names, ClinVar, ACMG/AMP hints, quality flags), CNV, SV, exon coverage, low-depth
intervals and tool versions. No third-party packages needed.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from brca12lib import tables, xlsx  # noqa: E402


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    with open(argv[0]) as handle:
        result = json.load(handle)
    xlsx.write(argv[1], tables.sample_sheets(result))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Join the long-read genotypes (ont_genotype.py) with the panel genotypes and the short-read
pipeline results of pathogenic/run_sv_carriers.sh for the same samples and sites.

Pipeline detection counts a Manta or small-variant call ("reported") and, for insertions, a
read-end insertion signature. Research use only.
"""
import csv
import sys


def main(ont_tsv, short_tsv, out_tsv):
    short = {(r["sample"], r["site"]): r for r in csv.DictReader(open(short_tsv), delimiter="\t")}
    rows = []
    for r in csv.DictReader(open(ont_tsv), delimiter="\t"):
        s = short.get((r["sample"], r["site"]))
        if s is None:
            continue
        panel = "present" if "1" in s["panel_genotype"] else "absent"
        pipeline = "call" if s["reported"] == "True" else ("read-end signal" if s["read_end_signal"] != "-" else "-")
        rows.append({"sample": r["sample"], "site": r["site"], "type": r["type"], "longReads": r["call"],
                     "supportingReads": f"{r['supportingReads']}/{r['spanningReads']}", "panel": panel,
                     "pipeline": pipeline})
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    agree = sum((r["longReads"] == "present") == (r["pipeline"] != "-") for r in rows)
    print(f"{len(rows)} sample-site pairs, pipeline agrees with long reads in {agree}")
    for r in rows:
        if (r["longReads"] == "present") != (r["pipeline"] != "-") or r["longReads"] != r["panel"]:
            print("\t".join(r.values()))


if __name__ == "__main__":
    main(*sys.argv[1:4])

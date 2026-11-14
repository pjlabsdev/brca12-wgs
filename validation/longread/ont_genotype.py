#!/usr/bin/env python3
"""Genotype structural-variant sites in long reads (1KG-ONT Vienna, Schloissnig et al. 2025).

Independent truth for the short-read validation: every site of a site table is checked in the
Oxford Nanopore reads of a sample. A read counts when it spans the site with a margin on both
sides; it supports the variant when its alignment carries an insertion (ALU, INS) or a deletion
(DEL) of at least half the expected length inside the window around the site.

A site is called present when at least --min-reads spanning reads support it and they make up at
least --min-fraction of the spanning reads. Research use only.
"""
import argparse
import csv
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

URL = "https://1kg-ont-vienna.s3.amazonaws.com/hg38/{sample}.hg38.cram"
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def read_table(path):
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("# ")]
    lines[0] = lines[0].lstrip("#")
    return list(csv.DictReader(lines, delimiter="\t"))


def events(pos, cigar):
    """Yield (kind, reference start, reference end, length) for insertions and deletions; return end."""
    out = []
    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op == "I":
            out.append(("INS", pos, pos, n))
        elif op in "DN":
            out.append(("DEL", pos, pos + n, n))
            pos += n
        elif op in "M=X":
            pos += n
    return out, pos


def genotype(sample, site, args):
    chrom, pos, end = site["chrom"], int(site["pos"]), int(site["end"])
    kind = "DEL" if site["type"] == "DEL" else "INS"
    need = max(50, int(site["length"]) // 2)
    lo, hi = pos - args.window, end + args.window
    cmd = ["samtools", "view", "-F", "0x904", "-q", str(args.mapq), "-T", args.reference,
           args.url.format(sample=sample), f"{chrom}:{lo}-{hi}"]
    for attempt in range(4):
        run = subprocess.run(cmd, capture_output=True, text=True)
        if run.returncode == 0:
            break
    else:
        return {"sample": sample, "site": site["id"], "status": "read error", "error": run.stderr.strip()[-200:]}
    spanning, sizes = 0, []
    for line in run.stdout.splitlines():
        f = line.split("\t", 6)
        start = int(f[3])
        evs, stop = events(start, f[5])
        if start > lo - args.margin or stop < hi + args.margin:
            continue
        spanning += 1
        size = sum(n for k, s, e, n in evs if k == kind and n >= args.min_piece and e >= lo and s <= hi)
        if size >= need:
            sizes.append(size)
    sizes.sort()
    fraction = round(len(sizes) / spanning, 2) if spanning else None
    if spanning < args.min_spanning:
        call = "not evaluable"
    elif len(sizes) >= args.min_reads and fraction >= args.min_fraction:
        call = "present"
    else:
        call = "absent"
    return {"sample": sample, "site": site["id"], "type": site["type"], "status": "ok", "call": call,
            "spanningReads": spanning, "supportingReads": len(sizes), "supportingFraction": fraction,
            "medianLength": sizes[len(sizes) // 2] if sizes else None}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--samples", required=True, help="file with one sample per line (first column)")
    p.add_argument("--sites", required=True, help="site table: id, chrom, pos, end, type, length")
    p.add_argument("--reference", required=True, help="GRCh38 FASTA with chr13 and chr17 (to decode CRAM)")
    p.add_argument("--out", required=True, help="output prefix (.json and .tsv)")
    p.add_argument("--url", default=URL)
    p.add_argument("--window", type=int, default=150, help="distance around the site in which the event may lie")
    p.add_argument("--margin", type=int, default=100, help="extra bases a read must cover beyond the window")
    p.add_argument("--mapq", type=int, default=10)
    p.add_argument("--min-piece", type=int, default=30, help="smallest single insertion/deletion that is summed")
    p.add_argument("--min-reads", type=int, default=2)
    p.add_argument("--min-fraction", type=float, default=0.10)
    p.add_argument("--min-spanning", type=int, default=4)
    p.add_argument("--threads", type=int, default=8)
    args = p.parse_args(argv)

    with open(args.samples) as fh:
        samples = [l.split()[0] for l in fh if l.strip() and not l.startswith("#")]
    sites = read_table(args.sites)
    jobs = [(s, site) for s in samples for site in sites]
    with ThreadPoolExecutor(args.threads) as pool:
        rows = list(pool.map(lambda j: genotype(j[0], j[1], args), jobs))
    with open(args.out + ".json", "w") as fh:
        json.dump({"source": "1KG-ONT Vienna, hg38 CRAM (s3://1kg-ont-vienna/hg38)",
                   "parameters": {k: getattr(args, k) for k in ("window", "margin", "mapq", "min_piece", "min_reads",
                                                                "min_fraction", "min_spanning")},
                   "rows": rows}, fh, indent=1)
    cols = ["sample", "site", "type", "status", "call", "spanningReads", "supportingReads", "supportingFraction", "medianLength"]
    with open(args.out + ".tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("-" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    bad = [r for r in rows if r["status"] != "ok"]
    print(f"{len(rows)} genotypes written, {len(bad)} read errors", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Frozen test of the known-allele genotyping (see frozen/PREREGISTRATION.md). Research use only.

For every sample of frozen/samples.tsv the SAM file SAMPLE.sam (samtools view of the two windows)
is run through bin/brca12_insertion_signals.py --gene-fasta, and the calls at the 126 bp tandem copy
(chr17:43128264) and the intronic BRCA2 Alu insertion (chr13:32348295) are compared with the 1KG-ONT
long-read genotypes of results/all_disc.*.

Usage: frozen_eval.py --gene-fasta genes.fa --sam-dir DIR --out results/frozen_test.tsv
"""
import argparse
import csv
import json
import math
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
REGIONS = "chr17\t43044294\t43170327\tBRCA1\nchr13\t32315076\t32400268\tBRCA2\n"
SITES = {"tandem126": ("ONT_BRCA1_INS_43128267", "chr17", 43128100, 43128450),
         "alu_brca2": ("ONT_BRCA2_INS_32348296", "chr13", 32348265, 32348330)}


def clopper_pearson(k, n, alpha=0.05):
    """Exact binomial confidence interval (bisection on the binomial tail)."""
    def tail(p, upper):
        terms = [math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(n + 1)]
        return sum(terms[k:]) if upper else sum(terms[:k + 1])
    def solve(f, target):
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if f(mid) < target else (lo, mid)
        return (lo + hi) / 2
    low = 0.0 if k == 0 else solve(lambda p: tail(p, True), alpha / 2)
    high = 1.0 if k == n else 1 - solve(lambda q: tail(1 - q, False), alpha / 2)
    return low, high


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--gene-fasta", required=True)
    p.add_argument("--sam-dir", required=True)
    p.add_argument("--samples", default=os.path.join(HERE, "frozen", "samples.tsv"))
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    disc = json.load(open(os.path.join(HERE, "results", "all_disc.json")))
    sites = {s["id"]: s for s in disc["sites"]}
    support = {}
    for r in csv.reader(open(os.path.join(HERE, "results", "all_disc.carriers.tsv")), delimiter="\t"):
        if r and not r[0].startswith("#"):
            support[(r[0], r[1])] = (int(r[2]), int(r[3]))
    samples = [l.split("\t")[0] for l in open(a.samples) if l.strip() and not l.startswith("#")]
    tmp = tempfile.mkdtemp()
    bed = os.path.join(tmp, "regions.bed")
    open(bed, "w").write(REGIONS)
    rows = []
    for s in samples:
        sam = os.path.join(a.sam_dir, s + ".sam")
        if not os.path.exists(sam):
            rows.append({"sample": s, "site": "", "status": "no short reads"})
            continue
        out = os.path.join(tmp, s + ".json")
        with open(sam) as fh:
            subprocess.run([sys.executable, os.path.join(REPO, "bin", "brca12_insertion_signals.py"), "--regions", bed,
                            "--gene-fasta", a.gene_fasta, "--out", out], stdin=fh, check=True)
        found = json.load(open(out))
        for allele, (site_id, chrom, lo, hi) in SITES.items():
            at = [r for r in found if r["chrom"] == chrom and lo < r["pos"] < hi]
            known = [r for r in at if "known allele" in r["signature"]]
            ev = known[0].get("evidence", []) if known else []
            truth = "carrier" if s in sites[site_id]["carriers"] else "non-carrier"
            spanning, supporting = support.get((s, site_id), ("", ""))
            rows.append({"sample": s, "site": allele, "status": "ok", "truth": truth,
                         "ontSpanning": spanning, "ontSupporting": supporting,
                         "called": "yes" if known else "no",
                         "realignment": "yes" if "allele realignment" in ev else "no",
                         "readEndSignature": "yes" if at and (at[0].get("clippedReadsLeftJunction") or
                                                              at[0].get("clippedReadsRightJunction")) else "no",
                         "evidence": ",".join(ev), "realignedReads": known[0].get("realignedReads", "") if known else "",
                         "zygosity": (known[0].get("zygosity") or {1: "heterozygous", 2: "homozygous"}.get(
                             known[0].get("haplotypeCopies"), "")) if known else ""})
    fields = ["sample", "site", "status", "truth", "ontSpanning", "ontSupporting", "called", "realignment",
              "readEndSignature", "evidence", "realignedReads", "zygosity"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    for allele in SITES:
        sub = [r for r in rows if r["site"] == allele]
        for col in ("called", "readEndSignature"):
            car = [r for r in sub if r["truth"] == "carrier"]
            non = [r for r in sub if r["truth"] == "non-carrier"]
            tp, fp = sum(r[col] == "yes" for r in car), sum(r[col] == "yes" for r in non)
            se, sp = clopper_pearson(tp, len(car)), clopper_pearson(len(non) - fp, len(non))
            print(f"{allele} {col}: sensitivity {tp}/{len(car)} (95 % CI {se[0]:.3f}-{se[1]:.3f}), "
                  f"specificity {len(non) - fp}/{len(non)} (95 % CI {sp[0]:.3f}-{sp[1]:.3f})")
        hom = [r for r in sub if r["truth"] == "carrier" and r["called"] == "yes" and r["ontSpanning"]]
        agree = sum((r["zygosity"] == "homozygous") == (r["ontSupporting"] / r["ontSpanning"] >= 0.9) for r in hom)
        print(f"{allele} zygosity agreement with long reads: {agree}/{len(hom)}")
    print("samples without short reads:", sum(r["status"] != "ok" for r in rows))


if __name__ == "__main__":
    main()

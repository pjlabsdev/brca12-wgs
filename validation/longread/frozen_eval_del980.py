#!/usr/bin/env python3
"""Addendum test of the 980 bp BRCA1 deletion (see frozen/ADDENDUM_del980.md). Research use only.

Usage: frozen_eval_del980.py --gene-fasta genes.fa --sam-dir DIR --out results/frozen_test_del980.tsv
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile

from frozen_eval import HERE, REGIONS, REPO, clopper_pearson

SITE = "ONT_BRCA1_DEL_43071728"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--gene-fasta", required=True)
    p.add_argument("--sam-dir", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    disc = json.load(open(os.path.join(HERE, "results", "all_disc.json")))
    carriers = next(s for s in disc["sites"] if s["id"] == SITE)["carriers"]
    support = {(r[0], r[1]): (int(r[2]), int(r[3])) for r in csv.reader(
        open(os.path.join(HERE, "results", "all_disc.carriers.tsv")), delimiter="\t") if r and not r[0].startswith("#")}
    samples = []
    for name in ("samples.tsv", "del980_extra_carriers.tsv"):
        samples += [l.split("\t")[0] for l in open(os.path.join(HERE, "frozen", name)) if l.strip() and not l.startswith("#")]
    tmp = tempfile.mkdtemp()
    bed = os.path.join(tmp, "regions.bed")
    open(bed, "w").write(REGIONS)
    rows = []
    for s in samples:
        out = os.path.join(tmp, s + ".json")
        with open(os.path.join(a.sam_dir, s + ".sam")) as fh:
            subprocess.run([sys.executable, os.path.join(REPO, "bin", "brca12_insertion_signals.py"), "--regions", bed,
                            "--gene-fasta", a.gene_fasta, "--out", out], stdin=fh, check=True)
        hit = [r for r in json.load(open(out)) if "deletion" in r["signature"] and r["chrom"] == "chr17"
               and abs(r["pos"] - 43071728) <= 50]
        span, sup = support.get((s, SITE), ("", ""))
        rows.append({"sample": s, "truth": "carrier" if s in carriers else "non-carrier", "ontSpanning": span,
                     "ontSupporting": sup, "called": "yes" if hit else "no",
                     "realignedReads": hit[0]["realignedReads"] if hit else "", "zygosity": hit[0]["zygosity"] if hit else ""})
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    car = [r for r in rows if r["truth"] == "carrier"]
    non = [r for r in rows if r["truth"] == "non-carrier"]
    tp, fp = sum(r["called"] == "yes" for r in car), sum(r["called"] == "yes" for r in non)
    se, sp = clopper_pearson(tp, len(car)), clopper_pearson(len(non) - fp, len(non))
    print(f"del980 carriers found {tp}/{len(car)} (95 % CI {se[0]:.3f}-{se[1]:.3f}), "
          f"non-carriers without call {len(non) - fp}/{len(non)} (95 % CI {sp[0]:.3f}-{sp[1]:.3f})")
    zy = [r for r in car if r["called"] == "yes" and r["ontSpanning"]]
    agree = sum((r["zygosity"] == "homozygous") == (r["ontSupporting"] / r["ontSpanning"] >= 0.9) for r in zy)
    print(f"del980 zygosity agreement with long reads: {agree}/{len(zy)}")


if __name__ == "__main__":
    main()

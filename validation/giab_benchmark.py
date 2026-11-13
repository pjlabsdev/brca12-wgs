#!/usr/bin/env python3
"""Compare brca12-wgs small-variant calls with a GIAB v4.2.1 truth set inside BRCA1/BRCA2.

Uses RTG vcfeval (haplotype-aware matching, genotype must agree) on the intersection of the
GIAB high-confidence regions with (a) the full gene loci and (b) the MANE exons +/- 20 bp.
Writes summary.tsv and summary.json with TP/FP/FN, precision, recall and F1 for SNVs and indels.
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys


def read_bed(path, names=None):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if names and (len(f) < 4 or not any(f[3].startswith(n) for n in names)):
                continue
            rows.append((f[0], int(f[1]), int(f[2])))
    return rows


def merge(rows):
    out = []
    for c, s, e in sorted(rows):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1] = (c, out[-1][1], max(out[-1][2], e))
        else:
            out.append((c, s, e))
    return out


def intersect(a, b):
    a, b, out = merge(a), merge(b), []
    for c, s, e in a:
        for c2, s2, e2 in b:
            if c == c2 and s2 < e and s < e2:
                out.append((c, max(s, s2), min(e, e2)))
    return merge(out)


def write_bed(path, rows):
    with open(path, "w") as handle:
        for c, s, e in rows:
            handle.write(f"{c}\t{s}\t{e}\n")
    return sum(e - s for _, s, e in rows)


def count(vcf, kind):
    if not os.path.exists(vcf):
        return 0
    out = subprocess.run(["bcftools", "view", "-H", "-v", kind, vcf], capture_output=True, text=True, check=True)
    return sum(1 for line in out.stdout.splitlines() if line.strip())


def wilson(k, n, z=1.96):
    if n == 0:
        return None, None
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return round(centre - half, 4), round(centre + half, 4)


def evaluate(args, label, regions, outdir):
    bed = os.path.join(outdir, f"{label}.eval.bed")
    bases = write_bed(bed, regions)
    run = os.path.join(outdir, f"vcfeval_{label}")
    shutil.rmtree(run, ignore_errors=True)        # vcfeval refuses an existing output directory
    subprocess.run(["rtg", "vcfeval", "-b", args.truth, "-c", args.calls, "-t", args.sdf,
                    "-e", bed, "-o", run, "--threads", str(args.threads)],
                   check=True, capture_output=True, text=True)
    rows = []
    for kind, name in (("snps", "SNV"), ("indels", "INDEL")):
        tp_base = count(f"{run}/tp-baseline.vcf.gz", kind)
        tp_call = count(f"{run}/tp.vcf.gz", kind)
        fp = count(f"{run}/fp.vcf.gz", kind)
        fn = count(f"{run}/fn.vcf.gz", kind)
        recall = tp_base / (tp_base + fn) if tp_base + fn else None
        precision = tp_call / (tp_call + fp) if tp_call + fp else None
        f1 = (2 * precision * recall / (precision + recall)) if precision and recall else None
        rows.append({"sample": args.sample, "region": label, "type": name, "evaluatedBases": bases,
                     "truthTotal": tp_base + fn, "TP": tp_base, "FP": fp, "FN": fn,
                     "recall": None if recall is None else round(recall, 4),
                     "recall95CI": wilson(tp_base, tp_base + fn),
                     "precision": None if precision is None else round(precision, 4),
                     "precision95CI": wilson(tp_call, tp_call + fp),
                     "F1": None if f1 is None else round(f1, 4)})
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--calls", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--truth-bed", required=True)
    p.add_argument("--regions", required=True, help="assets/brca12_regions.bed")
    p.add_argument("--exons", required=True, help="assets/brca12_exons.bed")
    p.add_argument("--sdf", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    conf = read_bed(args.truth_bed)
    genes = read_bed(args.regions, names=("BRCA1", "BRCA2"))
    exons = [(c, max(0, s - 20), e + 20) for c, s, e in read_bed(args.exons, names=("BRCA1", "BRCA2"))]
    rows = evaluate(args, "gene_loci", intersect(conf, genes), args.outdir)
    rows += evaluate(args, "exons_pm20", intersect(conf, exons), args.outdir)
    covered = sum(e - s for _, s, e in intersect(conf, genes))
    total = sum(e - s for _, s, e in merge(genes))
    with open(os.path.join(args.outdir, "summary.json"), "w") as handle:
        json.dump({"sample": args.sample, "geneBases": total, "geneBasesInGiabConfident": covered,
                   "results": rows}, handle, indent=2)
    keys = ["sample", "region", "type", "evaluatedBases", "truthTotal", "TP", "FP", "FN", "recall", "precision", "F1"]
    with open(os.path.join(args.outdir, "summary.tsv"), "w") as handle:
        handle.write("\t".join(keys) + "\n")
        for r in rows:
            handle.write("\t".join("" if r[k] is None else str(r[k]) for k in keys) + "\n")
    for r in rows:
        print("\t".join("" if r[k] is None else str(r[k]) for k in keys))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Random exon-affecting deletions and inversions for the frozen read-pair SV test (deterministic, seed below).

Deletions: 40, size log-uniform 50 bp to 20 kb, each removes at least part of a coding exon.
Inversions: 20, size log-uniform 1 to 20 kb, each contains a whole coding exon.
Exons: MANE Select, assets/brca12_exons.bed. Events lie inside the gene region, 2 kb from its ends.
Output: events.tsv in the format of ../pathogenic/spike.py (one heterozygous event per sample).
Usage: design.py EXON_BED OUT_TSV
"""
import math
import random
import sys

SEED = 202611151
REGIONS = {"BRCA1": ("chr17", 43044294, 43170327), "BRCA2": ("chr13", 32315076, 32400268)}


def main():
    exons = []
    for line in open(sys.argv[1]):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        exons.append((f[3].split("_")[0], f[0], int(f[1]) + 1, int(f[2])))
    rng = random.Random(SEED)
    rows = []
    for kind, n, lo, hi in (("LDEL", 40, 50, 20000), ("LINV", 20, 1000, 20000)):
        for i in range(1, n + 1):
            while True:
                gene, chrom, es, ee = rng.choice(exons)
                size = int(round(math.exp(rng.uniform(math.log(lo), math.log(hi)))))
                if kind == "LINV":
                    if size < ee - es + 200:
                        continue
                    start = rng.randint(ee - size + 100, es - 100)
                else:
                    # overlap the exon by at least 1 bp
                    start = rng.randint(es - size + 1, ee)
                end = start + size - 1
                _, rs, re_ = REGIONS[gene]
                if rs + 2000 < start and end < re_ - 2000:
                    break
            tag = "DEL" if kind == "LDEL" else "INV"
            what = "deletion" if kind == "LDEL" else "inversion"
            rows.append(f"{tag}_{i:02d}\t{gene}\t{kind}\t{chrom}\t{start}\t{end}\t{size}\t"
                        f"simulated {what}, random breakpoints, affects {gene} exon at {es}-{ee}\n")
    with open(sys.argv[2], "w") as out:
        out.write("#sample\tgene\tkind\tchrom\tstart\tend\tsize\tnote\n")
        out.writelines(rows)


if __name__ == "__main__":
    main()

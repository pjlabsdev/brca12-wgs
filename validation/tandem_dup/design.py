#!/usr/bin/env python3
"""Random exon-spanning tandem duplications for the frozen read-pair test (deterministic, seed below).

Size log-uniform between 1 and 20 kb; each event contains at least one coding exon of BRCA1 or BRCA2
(MANE Select, assets/brca12_exons.bed) and lies inside the gene region. Output: events.tsv in the
format of ../pathogenic/spike.py (one heterozygous event per sample).
Usage: design.py EXON_BED OUT_TSV
"""
import math
import random
import sys

SEED = 202611142
N = 40
REGIONS = {"BRCA1": ("chr17", 43044294, 43170327), "BRCA2": ("chr13", 32315076, 32400268)}


def main():
    exons = []
    for line in open(sys.argv[1]):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        exons.append((f[3].split("_")[0], f[0], int(f[1]) + 1, int(f[2])))
    rng = random.Random(SEED)
    with open(sys.argv[2], "w") as out:
        out.write("#sample\tgene\tkind\tchrom\tstart\tend\tsize\tnote\n")
        for i in range(1, N + 1):
            while True:
                gene, chrom, es, ee = rng.choice(exons)
                size = int(round(math.exp(rng.uniform(math.log(1000), math.log(20000)))))
                if size < ee - es + 200:
                    continue
                start = rng.randint(ee - size + 100, es - 100)
                end = start + size - 1
                _, rs, re_ = REGIONS[gene]
                if rs + 2000 < start and end < re_ - 2000:
                    break
            out.write(f"TDUP_{i:02d}\t{gene}\tLDUP\t{chrom}\t{start}\t{end}\t{size}\t"
                      f"simulated tandem duplication, random breakpoints, contains {gene} exon at {es}-{ee}\n")


if __name__ == "__main__":
    main()

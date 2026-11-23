#!/usr/bin/env python3
"""Builds the annotation tables that brca12-wgs 2.0 carries in the binary (assets/*.tsv.gz).

  clinvar_brca12.tsv.gz  ClinVar GRCh38 VCF records in the BRCA1/BRCA2 gene regions:
                         chrom pos ref alt variation_id clnsig review_status
  gnomad_brca12.tsv.gz   gnomAD records (v2.1.1 exomes lifted to GRCh38, v3.1.2 genomes):
                         chrom pos ref alt dataset present faf95 population non_cancer
                         The FAF is chosen exactly as brca12-wgs 1.1 (bin/brca12lib/evidence.py):
                         highest faf95 over populations, founder/mixed populations skipped,
                         non-cancer subset preferred where the dataset has one.

usage: build_annotation_tables.py CLINVAR.vcf.gz GNOMAD_EXOMES.vcf.gz GNOMAD_GENOMES.vcf.gz OUTDIR
"""
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
from brca12lib import evidence  # noqa: E402

REGIONS = {"13": (32315000, 32400300), "17": (43044000, 43170400)}


def in_region(chrom, pos):
    c = chrom[3:] if chrom.startswith("chr") else chrom
    r = REGIONS.get(c)
    return r is not None and r[0] <= pos <= r[1], "chr" + c


def clinvar(path, out):
    date = "unknown"
    rows = []
    with gzip.open(path, "rt") as h:
        for line in h:
            if line.startswith("##fileDate="):
                date = line.strip().split("=", 1)[1]
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 8)
            ok, chrom = in_region(f[0], int(f[1]))
            if not ok or f[4] in (".", "") or re.search(r"[^ACGTN,]", f[4]):
                continue
            info = evidence._info(f[7])
            sig = str(info.get("CLNSIG", "-"))
            rev = str(info.get("CLNREVSTAT", "-"))
            for alt in f[4].split(","):
                rows.append((chrom, int(f[1]), f[3], alt, f[2], sig, rev))
    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
    with gzip.open(out, "wt", compresslevel=9) as o:
        o.write(f"# ClinVar GRCh38 {date}, BRCA1/BRCA2 gene regions\n")
        for r in rows:
            o.write("\t".join(map(str, r)) + "\n")
    return len(rows), date


def gnomad(files, out):
    rows = []
    for label, path in files:
        with gzip.open(path, "rt") as h:
            for line in h:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t", 8)
                ok, chrom = in_region(f[0], int(f[1]))
                if not ok:
                    continue
                info = evidence._info(f[7])
                ac_keys = [k for k in info if re.fullmatch(r"(non_cancer_)?AC|AC_non_cancer", k)]
                present = any((evidence._float(info[k]) or 0) > 0 for k in ac_keys) or not ac_keys
                cands = []
                for key, value in info.items():
                    tokens = key.lower().replace("-", "_").split("_")
                    if "faf95" not in tokens or "joint" in tokens or evidence._SEX & set(tokens):
                        continue
                    pops = [t for t in tokens if t not in evidence._NOT_A_POPULATION]
                    pop = pops[-1] if pops else "max"
                    if pop in evidence.SKIP_POPULATIONS:
                        continue
                    faf = evidence._float(value)
                    if faf is not None:
                        cands.append(("cancer" in tokens, faf, pop))
                nc, faf, pop = (False, None, "-")
                if cands:
                    pool = [c for c in cands if c[0]] or cands
                    nc, faf, pop = max(pool, key=lambda c: c[1])
                if not present and not faf:
                    continue
                for alt in f[4].split(","):
                    rows.append((chrom, int(f[1]), f[3], alt, label, int(present), "-" if faf is None else repr(faf), pop, int(nc)))
    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3], r[4]))
    with gzip.open(out, "wt", compresslevel=9) as o:
        o.write("# gnomAD " + ", ".join(l for l, _ in files) + "; FAF95 as brca12-wgs 1.1\n")
        for r in rows:
            o.write("\t".join(map(str, r)) + "\n")
    return len(rows)


if __name__ == "__main__":
    cv, ex, ge, outdir = sys.argv[1:5]
    n, date = clinvar(cv, os.path.join(outdir, "clinvar_brca12.tsv.gz"))
    m = gnomad([("gnomAD v2.1.1 exomes", ex), ("gnomAD v3.1.2 genomes", ge)], os.path.join(outdir, "gnomad_brca12.tsv.gz"))
    print(f"ClinVar {date}: {n} alleles; gnomAD: {m} records")

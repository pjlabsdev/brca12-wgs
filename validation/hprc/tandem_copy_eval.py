#!/usr/bin/env python3
"""Tandem-copy signature at the common 126 bp insertion upstream of BRCA1 exon 1 (chr17:43128264):
clipped reads at both unit borders, the marker base of the extra copy (chr17:43128365 C) and the
linked haplotype (tag SNPs 3 to 21 kb downstream).

Runs bin/brca12_insertion_signals.py --gene-fasta on the BAM of every given pipeline output
directory (or on SAM files SAMPLE.sam, samtools view of the gene regions, with --sam) and compares
the result at the site with a truth table (sample, carrier/reference). Signatures elsewhere in the
genes are counted as well. Research use only.

Usage: tandem_copy_eval.py --reference ref.fa --truth truth.tsv --out result.tsv RUN_DIR...
       tandem_copy_eval.py --gene-fasta genes.fa --truth truth.tsv --out result.tsv --sam SAMPLE.sam...
"""
import argparse
import csv
import glob
import json
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SITE = ("chr17", 43128100, 43128450)
REGIONS = "chr17\t43044294\t43170327\tBRCA1\nchr13\t32315076\t32400268\tBRCA2\n"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--reference")
    p.add_argument("--gene-fasta", help="samtools faidx output of both gene regions (instead of --reference)")
    p.add_argument("--sam", action="store_true", help="inputs are SAM files SAMPLE.sam, not run directories")
    p.add_argument("--truth", required=True, help="TSV: sample, truth (carrier/reference)")
    p.add_argument("--out", required=True)
    p.add_argument("runs", nargs="+")
    a = p.parse_args(argv)
    truth = {r[0]: r[1] for r in csv.reader(open(a.truth), delimiter="\t") if r and not r[0].startswith("#")}
    tmp = tempfile.mkdtemp()
    bed, fa = os.path.join(tmp, "regions.bed"), os.path.join(tmp, "genes.fa")
    open(bed, "w").write(REGIONS)
    if a.gene_fasta:
        fa = a.gene_fasta
    else:
        with open(fa, "w") as fh:
            subprocess.run(["samtools", "faidx", a.reference, "chr17:43044295-43170327", "chr13:32315077-32400268"],
                           stdout=fh, check=True)
    rows = []
    for run in a.runs:
        sample = os.path.basename(run.rstrip("/"))
        if a.sam:
            sample = sample[:-4] if sample.endswith(".sam") else sample
        if sample not in truth or any(r["sample"] == sample for r in rows):
            continue
        out = os.path.join(tmp, sample + ".json")
        cmd = [sys.executable, os.path.join(REPO, "bin", "brca12_insertion_signals.py"), "--regions", bed,
               "--gene-fasta", fa, "--out", out]
        if a.sam:
            with open(run) as fh:
                subprocess.run(cmd, stdin=fh, check=True)
        else:
            bam = glob.glob(os.path.join(run, "alignment", "*.markdup.bam"))[0]
            view = subprocess.Popen(["samtools", "view", "-L", bed, bam], stdout=subprocess.PIPE)
            subprocess.run(cmd, stdin=view.stdout, check=True)
            view.wait()
        found = [r for r in json.load(open(out)) if "tandem copy" in r["signature"]]
        at = [r for r in found if r["chrom"] == SITE[0] and SITE[1] < r["pos"] < SITE[2]]
        known = [r for r in at if "known allele" in r["signature"]]
        rows.append({"sample": sample, "truth": truth[sample], "signature": "yes" if at else "no",
                     "knownAllele": "yes" if known else "no",
                     "evidence": ",".join(known[0]["evidence"]) if known else "",
                     "haplotypeTags": known[0].get("haplotypeTags", "") if known else "",
                     "haplotypeCopies": known[0].get("haplotypeCopies", "") if known else "",
                     "length": at[0]["duplicatedLength"] if at else "",
                     "reads": f"{at[0]['clippedReadsLeftJunction']}+{at[0]['clippedReadsRightJunction']}" if at else "",
                     "markerReads": at[0].get("markerReads", "") if at else "",
                     "elsewhere": len(found) - len(at)})
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for t in ("carrier", "reference"):
        sub = [r for r in rows if r["truth"] == t]
        print(f"{t}: signature in {sum(r['signature'] == 'yes' for r in sub)} of {len(sub)}, "
              f"known allele (2 of 3 signals) in {sum(r['knownAllele'] == 'yes' for r in sub)}")
    print("signatures elsewhere:", sum(r["elsewhere"] for r in rows))


if __name__ == "__main__":
    main()

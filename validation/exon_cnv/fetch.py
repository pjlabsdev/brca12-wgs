"""Collect per-exon read intervals (MAPQ >= 20, no duplicates/secondary/supplementary) for BRCA1/BRCA2
exons and depth in control windows from a 1000 Genomes 30x CRAM. Research use only.

usage: fetch.py REF.fa EXONS.bed OUTDIR SAMPLE RUN PREFIX
Output OUTDIR/SAMPLE.json: {"exons": {name: [[start, end, pairhash], ...]}, "controls": {name: mean depth}}"""
import json, os, sys, zlib
import pysam

S3 = "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage"
REF, BED, OUT, SAMPLE, RUN, PREFIX = sys.argv[1:7]
PAD = 150
CONTROLS = {  # single-copy windows on the same chromosomes, away from the two loci
    "c17_TP53": ("chr17", 7668000, 7672000), "c17_q12": ("chr17", 39700000, 39704000),
    "c17_q21b": ("chr17", 50000000, 50004000), "c17_q23": ("chr17", 61000000, 61004000),
    "c13_FLT3": ("chr13", 28000000, 28004000), "c13_RB1": ("chr13", 48300000, 48304000),
    "c13_q31": ("chr13", 90000000, 90004000), "c13_q33": ("chr13", 100000000, 100004000),
}
exons = [l.split("\t")[:4] for l in open(BED) if l.strip() and not l.startswith("#")]
out = os.path.join(OUT, SAMPLE + ".json")
if os.path.exists(out):
    sys.exit(0)
ok = lambda r: not (r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail
                    or r.mapping_quality < 20)
res = {"sample": SAMPLE, "exons": {}, "controls": {}}
with pysam.AlignmentFile(f"{S3}/{PREFIX}/{RUN}/{SAMPLE}.final.cram", "rc", reference_filename=REF) as f:
    for chrom, s, e, name in exons:
        s, e = int(s), int(e)
        res["exons"][name] = [[r.reference_start, r.reference_end, zlib.crc32(r.query_name.encode()) & 0xffff]
                              for r in f.fetch(chrom, max(0, s - PAD), e + PAD) if ok(r)]
    for name, (chrom, s, e) in CONTROLS.items():
        tot = 0
        for r in f.fetch(chrom, s, e):
            if ok(r):
                tot += max(0, min(r.reference_end, e) - max(r.reference_start, s))
        res["controls"][name] = tot / (e - s)
tmp = out + ".tmp"
json.dump(res, open(tmp, "w"))
os.rename(tmp, out)

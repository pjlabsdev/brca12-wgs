"""Depth per exon window (exon +/- 200 bp, reads overlapping it, MAPQ >= 20, no duplicates) and control
windows from a 1000 Genomes 30x CRAM, for the real-data screen. Research use only.
usage: screen_fetch.py REF.fa EXONS.bed OUTDIR SAMPLE RUN PREFIX"""
import json, os, sys
import pysam
sys.path.insert(0, sys.argv[7] if len(sys.argv) > 7 else ".")
from fetch_controls import CONTROLS
S3 = "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage"
REF, BED, OUT, SAMPLE, RUN, PREFIX = sys.argv[1:7]
PAD = 200
out = os.path.join(OUT, SAMPLE + ".json")
if os.path.exists(out):
    sys.exit(0)
ok = lambda r: not (r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail or r.mapping_quality < 20)
res = {"sample": SAMPLE, "depth": {}, "controls": {}}
with pysam.AlignmentFile(f"{S3}/{PREFIX}/{RUN}/{SAMPLE}.final.cram", "rc", reference_filename=REF) as f:
    for l in open(BED):
        if not l.strip() or l.startswith("#"):
            continue
        chrom, s, e, name = l.split("\t")[:4]
        s, e = int(s) - PAD, int(e) + PAD
        tot = sum(max(0, min(r.reference_end, e) - max(r.reference_start, s)) for r in f.fetch(chrom, s, e) if ok(r))
        res["depth"][name] = tot / (e - s)
    for name, c in CONTROLS.items():
        tot = sum(max(0, min(r.reference_end, c["end"]) - max(r.reference_start, c["start"]))
                  for r in f.fetch(c["chrom"], c["start"], c["end"]) if ok(r))
        res["controls"][name] = tot / (c["end"] - c["start"])
tmp = out + ".tmp"; json.dump(res, open(tmp, "w")); os.rename(tmp, out)

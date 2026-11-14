"""Per-exon reads with the alleles they carry at the sample's phased heterozygous SNVs (1000 Genomes
phased panel), for haplotype-aware simulation of exon deletions and duplications. Research use only.

usage: fetch_alleles.py REF.fa EXONS.bed PANEL_VCF_DIR OUTDIR SAMPLE RUN PREFIX
PANEL_VCF_DIR holds the .tbi files of the 1kGP phased panel (chr13, chr17); the VCFs are read remotely.
Output OUTDIR/SAMPLE.json:
  {"het": {exon: [[pos1, ref, alt, hap_of_alt], ...]},
   "exons": {exon: [[start, end, pairhash, {pos1: base}], ...]}, "controls": {...}}"""
import json, os, sys, zlib
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fetch_controls import CONTROLS  # noqa: E402

S3 = "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage"
VCF = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
       "20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
REF, BED, TBI, OUT, SAMPLE, RUN, PREFIX = sys.argv[1:8]
PAD = 200
out = os.path.join(OUT, SAMPLE + ".json")
if os.path.exists(out):
    sys.exit(0)
exons = [l.split("\t")[:4] for l in open(BED) if l.strip() and not l.startswith("#")]
ok = lambda r: not (r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail
                    or r.mapping_quality < 20)
vcfs = {c: pysam.VariantFile(VCF.format(c), index_filename=os.path.join(TBI, os.path.basename(VCF.format(c)) + ".tbi"))
        for c in ("chr13", "chr17")}
res = {"sample": SAMPLE, "het": {}, "exons": {}, "controls": {}}
with pysam.AlignmentFile(f"{S3}/{PREFIX}/{RUN}/{SAMPLE}.final.cram", "rc", reference_filename=REF) as f:
    for chrom, s, e, name in exons:
        s, e = int(s) - PAD, int(e) + PAD
        het = []
        for v in vcfs[chrom].fetch(chrom, s, e):
            if len(v.ref) != 1 or not v.alts or len(v.alts) != 1 or len(v.alts[0]) != 1:
                continue
            smp = v.samples[SAMPLE]
            gt = smp["GT"]
            if gt in ((0, 1), (1, 0)) and smp.phased:
                het.append([v.pos, v.ref, v.alts[0], 0 if gt == (1, 0) else 1])
        res["het"][name] = het
        sites = {h[0] - 1: h[0] for h in het}
        reads = []
        for r in f.fetch(chrom, max(0, s - 150), e + 150):
            if not ok(r):
                continue
            alleles = {}
            if sites:
                q = r.query_sequence
                for qp, rp in r.get_aligned_pairs(matches_only=True):
                    if rp in sites and r.query_qualities[qp] >= 13:
                        alleles[sites[rp]] = q[qp]
            reads.append([r.reference_start, r.reference_end, zlib.crc32(r.query_name.encode()) & 0xffff, alleles])
        res["exons"][name] = reads
    for name, c in CONTROLS.items():
        tot = 0
        for r in f.fetch(c["chrom"], c["start"], c["end"]):
            if ok(r):
                tot += max(0, min(r.reference_end, c["end"]) - max(r.reference_start, c["start"]))
        res["controls"][name] = tot / (c["end"] - c["start"])
tmp = out + ".tmp"
json.dump(res, open(tmp, "w"))
os.rename(tmp, out)

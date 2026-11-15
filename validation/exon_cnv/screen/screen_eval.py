"""Real-data screen: cohort calls/hints in unmodified 1000 Genomes samples and read-level check of
each hit by the rule in CONFIRMATION_RULE.md (fixed before the screen). Research use only.

usage: screen_eval.py REPO REF.fa SAMPLES.tsv DATADIR OUT.tsv"""
import json, os, statistics, sys
import pysam

REPO, REF, SAMPLES, DATA, OUT = sys.argv[1:6]
sys.path.insert(0, os.path.join(REPO, "bin"))
from brca12lib import exon_cnv  # noqa: E402

S3 = "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage"
panel = exon_cnv.load_panel(os.path.join(REPO, "assets", "exon_cnv_panel_1kg30x.json"))
ok = lambda r: not (r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail
                    or r.mapping_quality < 20)


def check(sample, run, prefix, seg):
    """Read-level evidence for one segment (CONFIRMATION_RULE.md)."""
    ex = [panel["exons"][n] for n in seg["exons"]]
    chrom = ex[0]["chrom"]
    s0, e0 = min(e["start"] for e in ex), max(e["end"] for e in ex)
    lo, hi = s0 - 20000, e0 + 20000
    f = pysam.AlignmentFile(f"{S3}/{prefix}/{run}/{sample}.final.cram", "rc", reference_filename=REF)
    B = 500
    bins = {}
    inserts, longpairs, outward, clips, outpairs = [], [], 0, {}, []
    for r in f.fetch(chrom, lo, hi):
        if not ok(r):
            continue
        for p in range(max(r.reference_start, lo), min(r.reference_end, hi)):
            bins[(p - lo) // B] = bins.get((p - lo) // B, 0) + 1
        if r.is_paired and r.is_read1 and not r.mate_is_unmapped and r.next_reference_name == chrom:
            tl = abs(r.template_length)
            if 0 < tl < 2000:
                inserts.append(tl)
            # orientation: FR = leftmost read forward, mate reverse
            left_fwd = (not r.is_reverse) if r.reference_start <= r.next_reference_start else (not r.mate_is_reverse)
            right_rev = r.mate_is_reverse if r.reference_start <= r.next_reference_start else r.is_reverse
            if left_fwd and right_rev:
                longpairs.append((min(r.reference_start, r.next_reference_start), max(r.reference_start, r.next_reference_start), tl))
            elif (not left_fwd) and (not right_rev) and r.is_reverse != r.mate_is_reverse:
                outward += 1
                outpairs.append((min(r.reference_start, r.next_reference_start), max(r.reference_start, r.next_reference_start)))
        c = r.cigartuples
        if c[0][0] == 4 and c[0][1] >= 15:
            clips[r.reference_start // 5 * 5] = clips.get(r.reference_start // 5 * 5, 0) + 1
        if c[-1][0] == 4 and c[-1][1] >= 15:
            clips[r.reference_end // 5 * 5] = clips.get(r.reference_end // 5 * 5, 0) + 1
    seg_bins = [k for k in bins if s0 - 200 <= lo + k * B + B and lo + k * B <= e0 + 200]
    out_bins = [k for k in bins if k not in seg_bins]
    # depth over the segment's exons +/- 200 bp vs median of bins outside
    inside = statistics.mean(bins[k] / B for k in seg_bins) if seg_bins else 0
    outside = statistics.median(bins[k] / B for k in out_bins) if out_bins else 0
    ratio = inside / outside if outside else None
    med = statistics.median(inserts) if inserts else 400
    deletion = seg["copyNumber"] < 2
    if deletion:
        bp_pairs = sum(1 for a, b, tl in longpairs if tl >= med + 500 and a < s0 and b > e0 - 1)
    else:
        bp_pairs = outward
    edge_clips = max([v for k, v in clips.items() if (abs(k - s0) <= 2000 or abs(k - e0) <= 2000)] or [0])
    depth_ok = ratio is not None and (ratio <= 0.70 if deletion else ratio >= 1.30)
    bp_ok = bp_pairs >= 3 or edge_clips >= 3
    status = "confirmed" if depth_ok and bp_ok else "depth-supported, breakpoint not seen" if depth_ok else "not confirmed"
    # post hoc, descriptive only (not part of the fixed rule): outward pairs whose mates sit at the two ends
    # of the segment (one within 3 kb before or inside its start, the other inside or within 3 kb after its end)
    local_out = sum(1 for a, b in outpairs if s0 - 3000 <= a <= e0 and s0 <= b <= e0 + 3000 and b - a >= 300)
    return {"chrom": chrom, "start": s0 + 1, "end": e0, "depthRatio": round(ratio, 2) if ratio else None,
            "breakpointPairs": bp_pairs, "maxEdgeClips": edge_clips, "medianInsert": med, "status": status,
            "localOutwardPairs": local_out}


rows = []
n = 0
for line in open(SAMPLES):
    sample, run, prefix = [x.strip() for x in line.rstrip("\n").split("\t")[:3]]
    path = os.path.join(DATA, sample + ".json")
    if not os.path.exists(path):
        continue
    n += 1
    d = json.load(open(path))
    depths = dict(d["depth"]); depths.update(d["controls"])
    res = exon_cnv.call(panel, depths)
    for tier in ("calls", "hints"):
        for seg in res.get(tier, []):
            ev = check(sample, run, prefix, seg)
            rows.append([sample, tier[:-1], seg["gene"], seg["type"], ",".join(seg["affectedExons"]), seg["quality"],
                         seg["meanLog2Deviation"], f"{ev['chrom']}:{ev['start']}-{ev['end']}", ev["depthRatio"],
                         ev["breakpointPairs"], ev["maxEdgeClips"], ev["status"], ev["localOutwardPairs"]])
            print("\t".join(map(str, rows[-1])), flush=True)
with open(OUT, "w") as o:
    o.write("sample\ttier\tgene\ttype\texons\tquality\tmean_log2\texon_span\tdepth_ratio\tbreakpoint_pairs\tmax_edge_clips\tstatus\tposthoc_local_outward_pairs\n")
    for r in rows:
        o.write("\t".join(map(str, r)) + "\n")
print(f"screened {n} samples; {len(rows)} segments", file=sys.stderr)

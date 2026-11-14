"""Exon copy number against a reference cohort. Research use only.

Per exon: mean depth over the exon +/- `pad` bases (MAPQ >= 20, no duplicates), divided by the
median depth of single-copy control windows of the same sample, log2. A panel of normal genomes of
the same assay (assets/exon_cnv_panel_1kg30x.json: 80 1000 Genomes 30x samples) gives each exon a
centre and a spread. A hidden Markov model over the exons of each gene (genomic order, copy-number
states 0-4) turns the deviations into segments; a segment is reported when its quality, -10 log10 of
the posterior probability that its exons have two copies, reaches `minQuality`.

Validation (validation/exon_cnv): rules frozen before the test; on 60 new 1000 Genomes samples with
simulated events 569/600 deletions and 209/360 duplications, no call in the 60 unmodified genomes."""
import json
import math
import statistics

STATES = (0, 1, 2, 3, 4)
SHIFT = {0: -3.0, 1: -1.0, 2: 0.0, 3: math.log2(1.5), 4: 1.0}
PRIOR = {0: 1e-5, 1: 5e-4, 2: 1 - 1e-3, 3: 5e-4 - 2e-5, 4: 1e-5}


def load_panel(path):
    with open(path) as handle:
        panel = json.load(handle)
    for key in ("exons", "controls", "pad", "pChange", "minQuality"):
        if key not in panel:
            raise ValueError(f"exon CNV panel {path}: missing '{key}'")
    return panel


def regions(panel):
    """BED rows (0-based) of all intervals whose depth is needed: exons +/- pad, control windows."""
    pad = panel["pad"]
    rows = [(e["chrom"], max(0, e["start"] - pad), e["end"] + pad, name) for name, e in panel["exons"].items()]
    rows += [(c["chrom"], c["start"], c["end"], name) for name, c in panel["controls"].items()]
    return sorted(rows, key=lambda r: (r[0], r[1]))


def mean_depths(depth_lines, rows):
    """Mean depth per interval from `samtools depth -a` lines (chrom, 1-based pos, depth)."""
    by_chrom = {}
    for chrom, start, end, name in rows:
        by_chrom.setdefault(chrom, []).append((start, end, name))
    sums = {name: 0 for *_, name in rows}
    seen = {name: 0 for *_, name in rows}
    for line in depth_lines:
        fields = line.split("\t")
        if len(fields) < 3 or fields[0] not in by_chrom:
            continue
        pos = int(fields[1]) - 1
        depth = int(fields[2])
        for start, end, name in by_chrom[fields[0]]:
            if start <= pos < end:
                sums[name] += depth
                seen[name] += 1
    return {name: (sums[name] / seen[name] if seen[name] else None) for name in sums}


def _logpdf(x, mu, sd):
    return -0.5 * ((x - mu) / sd) ** 2 - math.log(sd * math.sqrt(2 * math.pi))


def _lse(values):
    m = max(values)
    return m + math.log(sum(math.exp(v - m) for v in values))


def hmm(ys, sds, p_change):
    """Forward-backward and Viterbi over exons. Returns per-exon posteriors and the Viterbi path."""
    n = len(ys)
    emis = [{c: _logpdf(ys[i], SHIFT[c], sds[i] if c != 0 else 1.0) for c in STATES} for i in range(n)]
    trans = {(a, b): math.log(1 - p_change) if a == b else math.log(p_change * PRIOR[b] / (1 - PRIOR[a]))
             for a in STATES for b in STATES}
    start = {c: math.log(PRIOR[c]) + emis[0][c] for c in STATES}
    fw = [start]
    for i in range(1, n):
        fw.append({b: emis[i][b] + _lse([fw[-1][a] + trans[(a, b)] for a in STATES]) for b in STATES})
    bw = [None] * n
    bw[-1] = {c: 0.0 for c in STATES}
    for i in range(n - 2, -1, -1):
        bw[i] = {a: _lse([trans[(a, b)] + emis[i + 1][b] + bw[i + 1][b] for b in STATES]) for a in STATES}
    total = _lse(list(fw[-1].values()))
    post = [{c: math.exp(fw[i][c] + bw[i][c] - total) for c in STATES} for i in range(n)]
    v, back = [dict(start)], []
    for i in range(1, n):
        row, br = {}, {}
        for b in STATES:
            a = max(STATES, key=lambda s: v[-1][s] + trans[(s, b)])
            row[b] = v[-1][a] + trans[(a, b)] + emis[i][b]
            br[b] = a
        v.append(row)
        back.append(br)
    path = [max(STATES, key=lambda c: v[-1][c])]
    for br in reversed(back):
        path.append(br[path[-1]])
    return post, path[::-1]


def call(panel, depths):
    """Exon table and copy-number segments for one sample from interval mean depths."""
    controls = [depths.get(name) for name in panel["controls"]]
    controls = [d for d in controls if d]
    result = {"enabled": True, "panel": panel.get("description"), "panelSize": panel.get("size"),
              "pad": panel["pad"], "minQuality": panel["minQuality"], "exons": [], "calls": []}
    if len(controls) < max(3, len(panel["controls"]) // 2):
        result.update(enabled=False, reason="too few control windows with depth")
        return result
    norm = statistics.median(controls)
    result["normaliserDepth"] = round(norm, 2)
    genes = {}
    for name, e in panel["exons"].items():
        genes.setdefault(e["gene"], []).append((e["start"], name, e))
    for gene, rows in genes.items():
        rows.sort(key=lambda r: r[0])
        rows = [r for r in rows if depths.get(r[1]) is not None]
        if not rows:
            continue
        ys = [math.log2(max(depths[n], 0.05) / norm) - e["median"] for _, n, e in rows]
        sds = [e["sd"] for _, _, e in rows]
        post, path = hmm(ys, sds, panel["pChange"])
        for (_, n, e), y, sd, p, s in zip(rows, ys, sds, post, path):
            result["exons"].append({"gene": gene, "exon": n, "label": e.get("label", n), "log2Deviation": round(y, 3),
                                    "z": round(y / sd, 2), "copyNumber": s, "pTwoCopies": round(p[2], 4)})
        i = 0
        while i < len(rows):
            if path[i] == 2:
                i += 1
                continue
            j = i
            while j + 1 < len(rows) and path[j + 1] == path[i]:
                j += 1
            p2 = max(1e-300, max(post[k][2] for k in range(i, j + 1)))
            quality = min(99, int(-10 * math.log10(p2)))
            if quality >= panel["minQuality"]:
                first, last = rows[i][2], rows[j][2]
                result["calls"].append({
                    "gene": gene, "copyNumber": path[i],
                    "type": "Deletion" if path[i] < 2 else "Duplikation",
                    "chrom": first["chrom"], "start": min(first["start"], last["start"]) + 1,
                    "end": max(first["end"], last["end"]),
                    "exons": [rows[k][1] for k in range(i, j + 1)],
                    "affectedExons": [rows[k][2].get("label", rows[k][1]) for k in range(i, j + 1)],
                    "quality": quality,
                    "meanLog2Deviation": round(sum(ys[i:j + 1]) / (j - i + 1), 3),
                })
            i = j + 1
    return result


def build_panel(samples, exons, controls, pad, p_change, min_quality, min_sd=0.07, description=None):
    """samples: list of {interval name: mean depth}; exons: {name: {chrom, start, end, gene, label}}."""
    ratios = []
    for depths in samples:
        norm = statistics.median(depths[c] for c in controls)
        ratios.append({n: math.log2(max(depths[n], 0.05) / norm) for n in exons})
    out = {}
    for name, e in exons.items():
        values = [r[name] for r in ratios]
        med = statistics.median(values)
        mad = statistics.median(abs(x - med) for x in values) * 1.4826
        out[name] = dict(e, median=round(med, 4), sd=round(max(mad, min_sd), 4))
    return {"description": description, "size": len(samples), "pad": pad, "pChange": p_change,
            "minQuality": min_quality, "exons": out, "controls": controls}

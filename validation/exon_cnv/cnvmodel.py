"""Exon copy number against a reference cohort (prototype, research use only).

Per sample and exon: mean depth (MAPQ >= 20) divided by the sample's median depth in single-copy
control windows -> log2 ratio. A panel of normal genomes of the same assay gives per exon a robust
centre (median) and spread (MAD). A test sample's deviation per exon is modelled by a hidden Markov
model over the exons of each gene in genomic order with copy-number states 1, 2 and 3 (and 0/4 for
completeness); Gaussian emissions with the exon's panel spread. Segments with a state other than 2
are reported with a quality: -10 log10 of the posterior probability that every exon in the segment
has copy number 2."""
import math
import statistics

STATES = [0, 1, 2, 3, 4]
SHIFT = {0: -3.0, 1: -1.0, 2: 0.0, 3: math.log2(1.5), 4: 1.0}
PRIOR = {0: 1e-5, 1: 5e-4, 2: 1 - 1e-3, 3: 5e-4 - 2e-5, 4: 1e-5}
P_CHANGE = 0.02


def exon_depth(intervals, start, end):
    """Mean depth over [start, end) from read intervals [[s, e, ...], ...]."""
    tot = 0
    for s, e, *_ in intervals:
        tot += max(0, min(e, end) - max(s, start))
    return tot / (end - start)


PAD = 100   # depth over exon +/- PAD: breakpoints of exon-level events lie in the introns


def log2_ratios(sample, exons, pad=None):
    pad = PAD if pad is None else pad
    norm = statistics.median(sample["controls"].values())
    out = {}
    for chrom, s, e, name in exons:
        d = exon_depth(sample["exons"][name], s - pad, e + pad)
        out[name] = math.log2(max(d, 0.05) / norm)
    return out, norm


def build_panel(panel_ratios, min_sd=0.07):
    names = panel_ratios[0].keys()
    ref = {}
    for n in names:
        v = [r[n] for r in panel_ratios]
        med = statistics.median(v)
        mad = statistics.median(abs(x - med) for x in v) * 1.4826
        ref[n] = (med, max(mad, min_sd))
    return ref


def _norm_logpdf(x, mu, sd):
    return -0.5 * ((x - mu) / sd) ** 2 - math.log(sd * math.sqrt(2 * math.pi))


def hmm(ys, sds, p_change=None):
    p_change = P_CHANGE if p_change is None else p_change
    """Forward-backward over exons; returns per-exon posterior {state: p} and Viterbi path."""
    n = len(ys)
    ls = [{c: math.log(PRIOR[c])} for c in STATES]
    emis = [{c: _norm_logpdf(ys[i], SHIFT[c], sds[i] if c != 0 else 1.0) for c in STATES} for i in range(n)]
    trans = {(a, b): math.log(1 - p_change) if a == b else math.log(p_change * PRIOR[b] / (1 - PRIOR[a]))
             for a in STATES for b in STATES}

    def lse(xs):
        m = max(xs)
        return m + math.log(sum(math.exp(x - m) for x in xs))
    fw = [{c: math.log(PRIOR[c]) + emis[0][c] for c in STATES}]
    for i in range(1, n):
        fw.append({b: emis[i][b] + lse([fw[-1][a] + trans[(a, b)] for a in STATES]) for b in STATES})
    bw = [None] * n
    bw[-1] = {c: 0.0 for c in STATES}
    for i in range(n - 2, -1, -1):
        bw[i] = {a: lse([trans[(a, b)] + emis[i + 1][b] + bw[i + 1][b] for b in STATES]) for a in STATES}
    total = lse([fw[-1][c] for c in STATES])
    post = [{c: math.exp(fw[i][c] + bw[i][c] - total) for c in STATES} for i in range(n)]
    # Viterbi
    v = [{c: math.log(PRIOR[c]) + emis[0][c] for c in STATES}]
    back = []
    for i in range(1, n):
        row, br = {}, {}
        for b in STATES:
            a = max(STATES, key=lambda a: v[-1][a] + trans[(a, b)])
            row[b] = v[-1][a] + trans[(a, b)] + emis[i][b]; br[b] = a
        v.append(row); back.append(br)
    path = [max(STATES, key=lambda c: v[-1][c])]
    for br in reversed(back):
        path.append(br[path[-1]])
    return post, path[::-1]


def call(ratios, panel, exons, min_q=20):
    """Segments of copy number != 2 per gene. exons: BED rows in genomic order."""
    calls = []
    genes = {}
    for chrom, s, e, name in exons:
        genes.setdefault(name.split("_")[0], []).append((chrom, s, e, name))
    for gene, rows in genes.items():
        rows.sort(key=lambda r: r[1])
        ys = [ratios[r[3]] - panel[r[3]][0] for r in rows]
        sds = [panel[r[3]][1] for r in rows]
        post, path = hmm(ys, sds)
        i = 0
        while i < len(rows):
            if path[i] == 2:
                i += 1; continue
            j = i
            while j + 1 < len(rows) and path[j + 1] == path[i]:
                j += 1
            p2 = max(1e-300, max(post[k][2] for k in range(i, j + 1)))
            q = min(99, int(-10 * math.log10(p2)))
            if q >= min_q:
                calls.append({"gene": gene, "cn": path[i], "exons": [rows[k][3] for k in range(i, j + 1)],
                              "q": q, "meanShift": round(sum(ys[i:j + 1]) / (j - i + 1), 3)})
            i = j + 1
    return calls


def baseline_hints(sample, exons, low=0.65, high=1.35):
    """Current pipeline rule: exon depth / chromosome median outside [0.65, 1.35]."""
    norm = statistics.median(sample["controls"].values())
    out = []
    for chrom, s, e, name in exons:
        r = exon_depth(sample["exons"][name], s, e) / norm
        if r < low or r > high:
            out.append((name, "loss" if r < low else "gain"))
    return out

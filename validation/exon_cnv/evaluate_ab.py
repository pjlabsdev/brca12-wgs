"""Exon deletions/duplications simulated on one real haplotype, depth-only vs depth + allele balance.
Research use only. All events are simulated (not patient samples).

Reads are assigned to a haplotype by the allele they carry at the sample's phased heterozygous SNVs
(1000 Genomes phased panel); read pairs without such a site get a random haplotype from a hash of
the read name. A deletion removes the reads of one haplotype inside the event, a duplication adds a
copy of them; allele counts at heterozygous sites therefore shift as in a real event. As in the
pipeline (HaplotypeCaller calls), a site counts as heterozygous only with >= 2 reads of each allele
and a minor allele fraction >= 0.15.

usage: evaluate_ab.py EXONS.bed DATADIR PANEL.tsv EVAL.tsv SEED [--min-q Q] [--no-ab] [--out FILE]"""
import argparse, json, math, os, random, statistics, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bin"))
sys.path.insert(0, HERE)
from brca12lib import exon_cnv, coverage  # noqa: E402
from fetch_controls import CONTROLS  # noqa: E402

ap = argparse.ArgumentParser()
for a in ("bed", "data", "panel", "eval"):
    ap.add_argument(a)
ap.add_argument("seed", type=int)
ap.add_argument("--min-q", type=int, default=10)
ap.add_argument("--p-change", type=float, default=0.05)
ap.add_argument("--pad", type=int, default=200)
ap.add_argument("--no-ab", action="store_true")
ap.add_argument("--min-allele-reads", type=int, default=2)
ap.add_argument("--min-maf", type=float, default=0.15)
ap.add_argument("--out")
a = ap.parse_args()

exons = {}
for line in open(a.bed):
    if line.strip() and not line.startswith("#"):
        c, s, e, n = line.split("\t")[:4]
        g, _, label = coverage.exon_label(n)
        exons[n] = {"chrom": c, "start": int(s), "end": int(e), "gene": g, "label": label}
genes = {}
for n, e in sorted(exons.items(), key=lambda x: (x[1]["chrom"], x[1]["start"])):
    genes.setdefault(e["gene"], []).append(n)
ids = lambda f: [l.split("\t")[0] for l in open(f) if l.strip()]
load = lambda s: json.load(open(os.path.join(a.data, s + ".json")))


def depth(reads, s, e):
    return sum(max(0, min(r[1], e) - max(r[0], s)) for r in reads) / (e - s)


def depths_of(d):
    out = {n: depth(d["exons"][n], x["start"] - a.pad, x["end"] + a.pad) for n, x in exons.items()}
    out.update(d["controls"])
    return out


def hap_of(read, het):
    """haplotype (0/1) of a read from its alleles, else None."""
    votes = []
    for pos, ref, alt, hap_alt in het:
        b = read[3].get(str(pos), read[3].get(pos))
        if b == alt:
            votes.append(hap_alt)
        elif b == ref:
            votes.append(1 - hap_alt)
    return votes[0] if votes and len(set(votes)) == 1 else None


def simulate(d, sel, kind, hap):
    out = {"exons": dict(d["exons"]), "controls": d["controls"], "het": d["het"]}
    for n in sel:
        reads, het = d["exons"][n], d["het"][n]
        pair = {}
        for r in reads:
            h = hap_of(r, het)
            if h is not None:
                pair[r[2]] = h
        def h_of(r):
            h = hap_of(r, het)
            return h if h is not None else pair.get(r[2], r[2] % 2)
        if kind == "DEL":
            out["exons"][n] = [r for r in reads if h_of(r) != hap]
        else:
            out["exons"][n] = reads + [r for r in reads if h_of(r) == hap]
    return out


def allele_counts(d):
    res = {}
    for n in exons:
        sites = []
        for pos, ref, alt, _ in d["het"].get(n, []):
            r_, a_ = 0, 0
            for read in d["exons"][n]:
                b = read[3].get(str(pos), read[3].get(pos))
                r_ += b == ref
                a_ += b == alt
            if min(r_, a_) >= a.min_allele_reads and min(r_, a_) / (r_ + a_) >= a.min_maf:
                sites.append((r_, a_))
        if sites:
            res[n] = sites
    return res


# panel from the panel samples (same window and normaliser as the pipeline)
panel = exon_cnv.build_panel([depths_of(load(s)) for s in ids(a.panel) if os.path.exists(os.path.join(a.data, s + ".json"))],
                             exons, CONTROLS, a.pad, a.p_change, a.min_q)
KINDS = [("DEL", 1), ("DEL", 2), ("DEL", 3), ("DEL", 5), ("DEL", "all"), ("DUP", 1), ("DUP", 2), ("DUP", 3)]
FIXED = [("BRCA1", "DEL", ["BRCA1_ex01", "BRCA1_ex02"], "PSEUDO_DEL"), ("BRCA1", "DUP", ["BRCA1_ex01", "BRCA1_ex02"], "PSEUDO_DUP")]
KNOWN_REAL = {"NA18949": {"BRCA1_ex14", "BRCA1_ex15"}}
rng = random.Random(a.seed)
per, fp, rows = {}, {"null_calls": 0, "other_calls": 0, "null_samples": 0}, []


def run(d):
    return exon_cnv.call(panel, depths_of(d), None if a.no_ab else allele_counts(d))["calls"]


for sid in [s for s in ids(a.eval) if os.path.exists(os.path.join(a.data, s + ".json"))]:
    base = load(sid)
    real = KNOWN_REAL.get(sid, set())
    c0 = [c for c in run(base) if not set(c["exons"]) & real]
    fp["null_samples"] += 1
    fp["null_calls"] += len(c0)
    rows.append((sid, "none", "", ";".join(f"{c['gene']}:{c['copyNumber']}:{len(c['exons'])}ex:q{c['quality']}" for c in c0)))
    for gene, names in genes.items():
        scen = [(k, n, f"{k}{n}") for k, n in KINDS] + [(k, sel, lab) for g, k, sel, lab in FIXED if g == gene]
        for kind, n, label in scen:
            if isinstance(n, list):
                sel = n
            elif n == "all":
                sel = names
            else:
                i = rng.randrange(0, len(names) - n + 1)
                sel = names[i:i + n]
            hap = rng.randrange(2)
            calls = run(simulate(base, sel, kind, hap))
            want = (lambda c: c < 2) if kind == "DEL" else (lambda c: c > 2)
            hit = [c for c in calls if c["gene"] == gene and want(c["copyNumber"]) and set(c["exons"]) & set(sel)]
            other = [c for c in calls if not (c["gene"] == gene and set(c["exons"]) & set(sel)) and not set(c["exons"]) & real]
            fp["other_calls"] += len(other)
            pk = per.setdefault(label, {"n": 0, "found": 0, "exact": 0})
            pk["n"] += 1
            pk["found"] += bool(hit)
            pk["exact"] += any(set(c["exons"]) == set(sel) for c in hit)
            rows.append((sid, f"{gene}:{label}", ",".join(sel), ";".join(f"{c['gene']}:{c['copyNumber']}:{len(c['exons'])}ex:q{c['quality']}" for c in calls)))

mode = "depth only" if a.no_ab else "depth + allele balance"
print(f"{mode}; panel {panel['size']}; eval samples {fp['null_samples']}; min_q {a.min_q}; simulated events, not patient samples")
dels = sum(v["found"] for k, v in per.items() if k.startswith("DEL")); ndel = sum(v["n"] for k, v in per.items() if k.startswith("DEL"))
dups = sum(v["found"] for k, v in per.items() if k.startswith("DUP")); ndup = sum(v["n"] for k, v in per.items() if k.startswith("DUP"))
print(f"deletions {dels}/{ndel}  duplications {dups}/{ndup}")
for k, v in per.items():
    print(f"  {k:10s} n={v['n']:3d} found {v['found']:3d} exact {v['exact']:3d}")
print("false calls:", fp)
if a.out:
    with open(a.out, "w") as o:
        o.write("sample\tevent\tevent_exons\tcalls\n")
        for r in rows:
            o.write("\t".join(r) + "\n")

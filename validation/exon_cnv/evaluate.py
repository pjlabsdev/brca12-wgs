"""Simulated exon deletions/duplications in real 1000 Genomes reads, cohort HMM vs current rule.
Research use only. Simulated events: read pairs inside the event are removed (deletion) or
duplicated (duplication) for one of two pseudo-haplotypes chosen by a hash of the read name, so
coverage profile, GC and mappability effects stay those of the real sample.

usage: evaluate.py EXONS.bed DATADIR PANEL.tsv EVAL.tsv SEED [--min-q Q] [--out FILE]"""
import argparse, copy, json, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cnvmodel as M

ap = argparse.ArgumentParser()
ap.add_argument("bed"); ap.add_argument("data"); ap.add_argument("panel"); ap.add_argument("eval")
ap.add_argument("seed", type=int); ap.add_argument("--min-q", type=int, default=20); ap.add_argument("--out")
ap.add_argument("--pad", type=int); ap.add_argument("--p-change", type=float); ap.add_argument("--prior-event", type=float)
a = ap.parse_args()
if a.pad is not None: M.PAD = a.pad
if a.p_change is not None: M.P_CHANGE = a.p_change
if a.prior_event is not None:
    pe = a.prior_event; M.PRIOR.update({1: pe, 3: pe, 0: pe / 50, 4: pe / 50}); M.PRIOR[2] = 1 - sum(v for k, v in M.PRIOR.items() if k != 2)

exons = [(c, int(s), int(e), n) for c, s, e, n, *_ in
         (l.rstrip("\n").split("\t") for l in open(a.bed) if l.strip() and not l.startswith("#"))]
load = lambda s: json.load(open(os.path.join(a.data, s + ".json")))
ids = lambda f: [l.split("\t")[0] for l in open(f) if l.strip()]
panel_ids = [s for s in ids(a.panel) if os.path.exists(os.path.join(a.data, s + ".json"))]
panel = M.build_panel([M.log2_ratios(load(s), exons)[0] for s in panel_ids])

genes = {}
for row in sorted(exons, key=lambda r: (r[0], r[1])):
    genes.setdefault(row[3].split("_")[0], []).append(row[3])
# real events found in the samples themselves (read-level confirmed); calls on them are not false
KNOWN_REAL = {"NA18949": {"BRCA1_ex14", "BRCA1_ex15"}}
# fixed scenarios: BRCA1 exon 1-2 deletion/duplication through the pseudogene BRCA1P1 (chr17:43123829-43167978)
FIXED = [("BRCA1", "DEL", ["BRCA1_ex01", "BRCA1_ex02"], "PSEUDO_DEL"), ("BRCA1", "DUP", ["BRCA1_ex01", "BRCA1_ex02"], "PSEUDO_DUP")]
KINDS = [("DEL", 1), ("DEL", 2), ("DEL", 3), ("DEL", 5), ("DEL", "all"), ("DUP", 1), ("DUP", 2), ("DUP", 3)]


def apply(sample, gene, kind, n, rng):
    s = copy.deepcopy(sample)
    names = genes[gene]
    if isinstance(n, list):
        sel = n
    elif n == "all":
        sel = names
    else:
        i = rng.randrange(0, len(names) - n + 1); sel = names[i:i + n]
    hap = rng.randrange(2)
    for name in sel:
        reads = s["exons"][name]
        if kind == "DEL":
            s["exons"][name] = [r for r in reads if r[2] % 2 != hap]
        else:
            s["exons"][name] = reads + [r for r in reads if r[2] % 2 == hap]
    return s, sel


rng = random.Random(a.seed)
res = {"cohort": {"tp": 0, "exact": 0, "fn": 0}, "baseline": {"tp": 0, "fn": 0}}
per_kind = {}
fp = {"cohort_null_calls": 0, "baseline_null_hints": 0, "null_samples": 0, "cohort_other_gene_calls": 0}
rows = []
eval_ids = [s for s in ids(a.eval) if os.path.exists(os.path.join(a.data, s + ".json"))]
for sid in eval_ids:
    base = load(sid)
    r0, _ = M.log2_ratios(base, exons)
    real = KNOWN_REAL.get(sid, set())
    c0 = [c for c in M.call(r0, panel, exons, a.min_q) if not set(c["exons"]) & real]
    b0 = [h for h in M.baseline_hints(base, exons) if h[0] not in real]
    fp["null_samples"] += 1; fp["cohort_null_calls"] += len(c0); fp["baseline_null_hints"] += len(b0)
    rows.append((sid, "none", "", ";".join(f"{c['gene']}:{c['cn']}:{len(c['exons'])}ex:q{c['q']}" for c in c0), len(b0)))
    for gene in genes:
        for kind, n, label in [(k, n, f"{k}{n}") for k, n in KINDS] + [(k, sel, lab) for g, k, sel, lab in FIXED if g == gene]:
            s, sel = apply(base, gene, kind, n, rng)
            r, _ = M.log2_ratios(s, exons)
            calls = M.call(r, panel, exons, a.min_q)
            want = (lambda c: c < 2) if kind == "DEL" else (lambda c: c > 2)
            hit = [c for c in calls if c["gene"] == gene and want(c["cn"]) and set(c["exons"]) & set(sel)]
            other = [c for c in calls if not (c["gene"] == gene and set(c["exons"]) & set(sel)) and not set(c["exons"]) & real]
            fp["cohort_other_gene_calls"] += len(other)
            k = label
            pk = per_kind.setdefault(k, {"n": 0, "cohort": 0, "exact": 0, "baseline": 0})
            pk["n"] += 1
            if hit:
                res["cohort"]["tp"] += 1; pk["cohort"] += 1
                if any(set(c["exons"]) == set(sel) for c in hit):
                    res["cohort"]["exact"] += 1; pk["exact"] += 1
            else:
                res["cohort"]["fn"] += 1
            bh = [h for h in M.baseline_hints(s, exons) if h[0] in sel and h[1] == ("loss" if kind == "DEL" else "gain")]
            if bh:
                res["baseline"]["tp"] += 1; pk["baseline"] += 1
            else:
                res["baseline"]["fn"] += 1
            rows.append((sid, f"{gene}:{k}", ",".join(sel), ";".join(f"{c['gene']}:{c['cn']}:{len(c['exons'])}ex:q{c['q']}" for c in calls), len(bh)))

print("panel", len(panel_ids), "eval samples", len(eval_ids), "min_q", a.min_q)
print("simulated events (simulated, not patient samples):", res)
for k, v in per_kind.items():
    print(f"  {k:6s} n={v['n']:3d} cohort {v['cohort']:3d} (exact {v['exact']:3d})  current rule {v['baseline']:3d}")
print("false calls:", fp)
if a.out:
    with open(a.out, "w") as o:
        o.write("sample\tevent\tevent_exons\tcohort_calls\tcurrent_rule_hits\n")
        for r in rows:
            o.write("\t".join(map(str, r)) + "\n")

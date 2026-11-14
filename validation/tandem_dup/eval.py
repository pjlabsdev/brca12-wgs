#!/usr/bin/env python3
"""Evaluate the frozen read-pair tandem duplication test (PREREGISTRATION.md).
Usage: eval.py OUTDIR  -> results/tandem_dup_test.tsv and results/tandem_dup_summary.txt"""
import glob
import gzip
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pathogenic", "founder"))
from eval import clopper_pearson  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = {}
for line in open(os.path.join(HERE, "events.tsv")):
    if not line.startswith("#"):
        f = line.split("\t")
        TRUTH[f[0]] = (f[3], int(f[4]), int(f[5]))
OTHER = {"SPIKE_L3": ("chr17", 43079800, 43085799), "SPIKE_L5": ("chr13", 32369000, 32373999)}


def match(chrom, s, e, truth, tol=1000):
    return truth and chrom == truth[0] and abs(s - truth[1]) <= tol and abs(e - truth[2]) <= tol


def ours(path):
    return [r for r in json.load(open(path)) if "tandem duplication" in r.get("signature", [])]


def manta(path, truth):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[6] != "PASS":
                continue
            info = f[7]
            end = int(info.split("END=")[1].split(";")[0]) if "END=" in info else int(f[1])
            if match(f[0], int(f[1]), end, truth) or (truth and f[0] == truth[0] and
                                                       (abs(int(f[1]) - truth[1]) <= 1000 or abs(int(f[1]) - truth[2]) <= 1000)):
                return info.split("SVTYPE=")[1].split(";")[0]
    return "no"


def main():
    out = sys.argv[1]
    rows = []
    for path in sorted(glob.glob(os.path.join(out, "res", "*.json"))):
        t = os.path.basename(path)[:-5]
        name, f = t.rsplit(".", 1)
        truth = TRUTH.get(name) or OTHER.get(name)
        calls = ours(path)
        hit = [c for c in calls if match(c["chrom"], c["pos"], c["end"], truth)]
        rows.append({"test": t, "sample": name, "set": "simulated", "reads": {"f1": "all", "f2": "1/2", "f3": "1/3"}[f],
                     "truth": "carrier" if name in TRUTH else ("other duplication" if name in OTHER else "non-carrier"),
                     "size": truth[2] - truth[1] + 1 if truth else 0,
                     "called": "yes" if hit else "no", "falseCalls": len(calls) - len(hit),
                     "support": f"{hit[0]['outwardPairs']} pairs / {hit[0]['splitReads']} split" if hit else "",
                     "manta": manta(os.path.join(out, "res", t + ".manta.vcf.gz"), truth) if truth else "n/a"})
    for path in sorted(glob.glob(os.path.join(out, "real", "*.json"))):
        s = os.path.basename(path)[:-5]
        calls = ours(path)
        rows.append({"test": s, "sample": s, "set": "real 1000 Genomes", "reads": "all", "truth": "non-carrier (to be checked)",
                     "size": 0, "called": "no", "falseCalls": len(calls),
                     "support": ";".join(f"{c['chrom']}:{c['pos']}-{c['end']} {c['outwardPairs']}p/{c['splitReads']}s" for c in calls),
                     "manta": "not run"})
    cols = list(rows[0])
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "tandem_dup_test.tsv"), "w") as h:
        h.write("\t".join(cols) + "\n")
        for r in rows:
            h.write("\t".join(str(r[c]) for c in cols) + "\n")

    def rate(sub, key):
        k = sum(1 for r in sub if r[key] not in ("no", "n/a"))
        lo, hi = clopper_pearson(k, len(sub))
        return f"{k}/{len(sub)} ({100 * lo:.1f}-{100 * hi:.1f} %)"
    car = [r for r in rows if r["truth"] == "carrier"]
    lines = ["Read-pair tandem duplication scan, frozen test. Simulated sets are SIMULATED reads (real HG002",
             "positions, fragment sizes and base qualities); the 1000 Genomes samples are real reads.", ""]
    lines.append(f"simulated carriers: read pairs {rate(car, 'called')}, Manta {rate(car, 'manta')}")
    for reads in ("all", "1/2", "1/3"):
        s2 = [r for r in car if r["reads"] == reads]
        lines.append(f"  reads {reads}: read pairs {rate(s2, 'called')}, Manta {rate(s2, 'manta')}")
    for lo_, hi_ in ((1000, 3000), (3000, 8000), (8000, 20001)):
        s2 = [r for r in car if lo_ <= r["size"] < hi_]
        lines.append(f"  size {lo_}-{hi_ - 1} bp: read pairs {rate(s2, 'called')}, Manta {rate(s2, 'manta')}")
    sim = [r for r in rows if r["set"] == "simulated"]
    lines.append(f"false calls in simulated samples (all depths): {sum(r['falseCalls'] for r in sim)} in {len(sim)} analyses")
    real = [r for r in rows if r["set"] != "simulated"]
    lines.append(f"real genomes with a call: {sum(1 for r in real if r['falseCalls'])} of {len(real)} "
                 "(each checked against ONT long reads, see README)")
    text = "\n".join(lines) + "\n"
    open(os.path.join(HERE, "results", "tandem_dup_summary.txt"), "w").write(text)
    print(text)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Evaluate the frozen BRCA1 exon 13 deletion test (see ADDENDUM_del3835.md).

Usage: eval_dup13.py OUTDIR  (written by run_del3835.sh and realneg_del3835.sh)
Writes results/del3835_test.tsv and results/del3835_test_summary.txt.
"""
import glob
import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval import clopper_pearson  # noqa: E402

START, END = 43080384, 43084218


def ours(path):
    rows = json.load(open(path))
    call = next((r for r in rows if "deletion" in r.get("signature", []) and r.get("clinicalSignificance") == "pathogenic"
                 and abs(r["pos"] - START) <= 50 and abs(r["end"] - END) <= 50), None)
    return call


def manta(path):
    hits = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            info = f[7]
            end = int(info.split("END=")[1].split(";")[0]) if "END=" in info else None
            near = f[0] == "chr17" and (abs(int(f[1]) - START) <= 50 or abs(int(f[1]) - END) <= 50 or
                                        (end is not None and abs(end - END) <= 50))
            if near and f[6] == "PASS":
                hits.append(info.split("SVTYPE=")[1].split(";")[0])
    return ",".join(sorted(set(hits)))


def main():
    out = sys.argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    rows = []
    for path in sorted(glob.glob(os.path.join(out, "res", "*.json"))):
        t = os.path.basename(path)[:-5]
        name, f = t.rsplit(".", 1)
        call = ours(path)
        rows.append({"test": t, "sample": name, "set": "simulated", "reads": {"f1": "all", "f2": "1/2", "f3": "1/3"}[f],
                     "depthAtSite": open(os.path.join(out, "res", t + ".depth")).read().strip() or "0",
                     "truth": "carrier" if name.startswith("DEL3835_") else "non-carrier",
                     "catalogCall": "yes" if call else "no", "realignedReads": call["realignedReads"] if call else "",
                     "zygosity": call["zygosity"] if call else "",
                     "manta": manta(os.path.join(out, "res", t + ".manta.vcf.gz")) or "no"})
    for path in sorted(glob.glob(os.path.join(out, "real", "*.json"))):
        s = os.path.basename(path)[:-5]
        call = ours(path)
        dpf = os.path.join(out, "real", s + ".depth")
        rows.append({"test": s, "sample": s, "set": "real 1000 Genomes", "reads": "all",
                     "depthAtSite": (open(dpf).read().strip() if os.path.exists(dpf) else "") or "0",
                     "truth": "non-carrier (no known carrier)", "catalogCall": "yes" if call else "no",
                     "realignedReads": call["realignedReads"] if call else "",
                     "zygosity": call["zygosity"] if call else "", "manta": "not run"})
    cols = list(rows[0])
    os.makedirs(os.path.join(here, "results"), exist_ok=True)
    with open(os.path.join(here, "results", "del3835_test.tsv"), "w") as handle:
        handle.write("\t".join(cols) + "\n")
        for r in rows:
            handle.write("\t".join(str(r[c]) for c in cols) + "\n")

    def rate(sub, key, any_=False):
        k = sum(1 for r in sub if (r[key] != "no" if any_ else r[key] == "yes"))
        lo, hi = clopper_pearson(k, len(sub))
        return f"{k}/{len(sub)} ({100 * lo:.1f}-{100 * hi:.1f} %)"
    lines = ["BRCA1 exon 13 deletion (3,835 bp), frozen test. Simulated sets are SIMULATED reads (real HG002 positions,",
             "fragment sizes and base qualities); the 1000 Genomes samples are real reads.", ""]
    car = [r for r in rows if r["truth"] == "carrier"]
    neg = [r for r in rows if r["truth"] == "non-carrier"]
    real = [r for r in rows if r["set"] != "simulated"]
    if car:
        lines.append(f"simulated carriers: catalog {rate(car, 'catalogCall')}, Manta {rate(car, 'manta', True)}")
        for reads in ("all", "1/2", "1/3"):
            s2 = [r for r in car if r["reads"] == reads]
            dps = sorted(int(r["depthAtSite"]) for r in s2)
            lines.append(f"  reads {reads} (median depth {dps[len(dps) // 2]}): catalog {rate(s2, 'catalogCall')}, "
                         f"Manta {rate(s2, 'manta', True)}")
        het = sum(1 for r in car if r["zygosity"] == "heterozygous")
        lines.append(f"  zygosity heterozygous (truth): {het}/{sum(1 for r in car if r['catalogCall'] == 'yes')} calls")
    if neg:
        lines.append(f"simulated non-carriers: catalog {rate(neg, 'catalogCall')}, Manta {rate(neg, 'manta', True)}")
        l3 = [r for r in neg if r["sample"] == "SPIKE_L2"]
        lines.append(f"  of these SPIKE_L2 (3,835 bp exon deletion, synthetic breakpoints 116 bp away): catalog {rate(l3, 'catalogCall')}")
    if real:
        lines.append(f"real non-carriers: catalog {rate(real, 'catalogCall')}")
    text = "\n".join(lines) + "\n"
    open(os.path.join(here, "results", "del3835_test_summary.txt"), "w").write(text)
    print(text)


if __name__ == "__main__":
    main()

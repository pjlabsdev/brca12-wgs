#!/usr/bin/env python3
"""Evaluate the frozen BRCA2 c.156_157insAlu test (see PREREGISTRATION.md).

Usage: eval.py OUTDIR  (OUTDIR as written by run.sh and realneg.sh)
Writes results/founder_test.tsv and prints the summary (results/founder_test_summary.txt).
"""
import glob
import gzip
import json
import os
import sys
from math import comb

SITE = 32319165
CARRIER_PREFIX = ("FOUNDER_", "SPIKE_L1")


def clopper_pearson(k, n, alpha=0.05):
    if n == 0:
        return (float("nan"), float("nan"))

    def cdf(x, p):
        return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(x + 1))

    def solve(f):
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if f(mid) else (lo, mid)
        return (lo + hi) / 2
    low = 0.0 if k == 0 else solve(lambda p: 1 - cdf(k - 1, p) < alpha / 2)
    high = 1.0 if k == n else solve(lambda p: cdf(k, p) > alpha / 2)
    return low, high


def ours(path):
    rows = json.load(open(path))
    near = [r for r in rows if r["chrom"].replace("chr", "") == "13" and abs(r["pos"] - SITE) <= 30]
    call = next((r for r in near if r.get("clinicalSignificance") == "pathogenic"), None)
    clips = any("Alu" in r.get("signature", []) and
                (r.get("clippedReadsLeftJunction", 0) + r.get("clippedReadsRightJunction", 0)) > 0 for r in near)
    return call, clips


def manta(path):
    if not os.path.exists(path):
        return ""
    hits = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[0] == "chr13" and abs(int(f[1]) - SITE) <= 50 and f[6] == "PASS":
                hits.append(f[7].split("SVTYPE=")[1].split(";")[0])
    return ",".join(sorted(set(hits)))


def main():
    out = sys.argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    rows = []
    for path in sorted(glob.glob(os.path.join(out, "res", "*.json"))):
        t = os.path.basename(path)[:-5]
        name, depth_set = t.rsplit(".", 1)
        call, clips = ours(path)
        dp = open(os.path.join(out, "res", t + ".depth")).read().strip() or "0"
        rows.append({"test": t, "sample": name, "set": "simulated", "reads": {"f1": "all", "f2": "1/2", "f3": "1/3"}[depth_set],
                     "depthAtSite": dp, "truth": "carrier" if name.startswith(CARRIER_PREFIX) else "non-carrier",
                     "catalogCall": "yes" if call else "no",
                     "realignedReads": call["realignedReads"] if call else "",
                     "zygosity": call["zygosity"] if call else "",
                     "readEndSignature": "yes" if clips else "no",
                     "manta": manta(os.path.join(out, "res", t + ".manta.vcf.gz")) or "no"})
    for path in sorted(glob.glob(os.path.join(out, "real", "*.json"))):
        s = os.path.basename(path)[:-5]
        call, clips = ours(path)
        dpf = os.path.join(out, "real", s + ".depth")
        dp = open(dpf).read().strip() if os.path.exists(dpf) else ""
        rows.append({"test": s, "sample": s, "set": "real 1000 Genomes", "reads": "all", "depthAtSite": dp or "0",
                     "truth": "non-carrier (no known carrier)", "catalogCall": "yes" if call else "no",
                     "realignedReads": call["realignedReads"] if call else "", "zygosity": call["zygosity"] if call else "",
                     "readEndSignature": "yes" if clips else "no", "manta": "not run"})
    cols = list(rows[0])
    os.makedirs(os.path.join(here, "results"), exist_ok=True)
    with open(os.path.join(here, "results", "founder_test.tsv"), "w") as handle:
        handle.write("\t".join(cols) + "\n")
        for r in rows:
            handle.write("\t".join(str(r[c]) for c in cols) + "\n")
    lines = ["BRCA2 c.156_157insAlu, frozen test. Simulated carriers and spike-in negatives are SIMULATED reads",
             "(real HG002 positions, fragment sizes and base qualities); the 1000 Genomes samples are real reads.", ""]

    def rate(sub, key, val="yes"):
        k = sum(1 for r in sub if (r[key] != "no" if val == "any" else r[key] == val))
        lo, hi = clopper_pearson(k, len(sub))
        return f"{k}/{len(sub)} ({100 * lo:.1f}-{100 * hi:.1f} %)"
    for label, sub in (("simulated carriers", [r for r in rows if r["truth"] == "carrier" and r["sample"].startswith("FOUNDER_")]),
                       ("SPIKE_L1 (genomic AluY copy)", [r for r in rows if r["sample"] == "SPIKE_L1"]),
                       ("simulated non-carriers", [r for r in rows if r["truth"] == "non-carrier"]),
                       ("real non-carriers", [r for r in rows if r["set"] != "simulated"])):
        if not sub:
            continue
        lines.append(f"{label}: catalog call {rate(sub, 'catalogCall')}, read-end signature {rate(sub, 'readEndSignature')}"
                     + (f", Manta {rate(sub, 'manta', 'any')}" if sub[0]["set"] == "simulated" else ""))
        if sub[0]["truth"] == "carrier":
            for reads in ("all", "1/2", "1/3"):
                s2 = [r for r in sub if r["reads"] == reads]
                if s2:
                    dps = sorted(int(r["depthAtSite"]) for r in s2)
                    lines.append(f"  reads {reads} (median depth {dps[len(dps) // 2]}): catalog {rate(s2, 'catalogCall')}, "
                                 f"read-end {rate(s2, 'readEndSignature')}, Manta {rate(s2, 'manta', 'any')}")
            het = sum(1 for r in sub if r["zygosity"] == "heterozygous")
            lines.append(f"  zygosity heterozygous (truth): {het}/{sum(1 for r in sub if r['catalogCall'] == 'yes')} calls")
    if any(r["sample"].startswith("FOUNDER_") for r in rows):
        lines.append("")
        lines.append("simulated carriers by element divergence from the consensus (all depths):")
        for d in ("00", "02", "05", "08", "10"):
            s2 = [r for r in rows if r["sample"].startswith(f"FOUNDER_d{d}_")]
            if s2:
                lines.append(f"  {int(d)} %: catalog {rate(s2, 'catalogCall')}, read-end {rate(s2, 'readEndSignature')}, "
                             f"Manta {rate(s2, 'manta', 'any')}")
    text = "\n".join(lines) + "\n"
    open(os.path.join(here, "results", "founder_test_summary.txt"), "w").write(text)
    print(text)


if __name__ == "__main__":
    main()

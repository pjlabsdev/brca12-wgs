#!/usr/bin/env python3
"""Find all insertions and deletions of 50 bp or more in BRCA1 and BRCA2 in the long reads of the
1KG-ONT Vienna cohort and genotype every sample at every site found.

Step 1 reads the two gene regions of each sample once and keeps, per read, its span and its large
insertions/deletions (cached as JSON per sample, so the run can be resumed).
Step 2 clusters events of the same kind that lie within --cluster bp of each other and whose
lengths differ by less than 30 %. A cluster becomes a site when at least one sample carries it.
Step 3 genotypes every sample at every site with the rule of ont_genotype.py: present when at
least --min-reads spanning reads support it and they make up at least --min-fraction.

The result is a truth set from reads alone; it is not a curated callset. Research use only.
"""
import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

from ont_genotype import URL, events

GENES = (("chr17", 43044295, 43170327, "BRCA1"), ("chr13", 32315077, 32400268, "BRCA2"))


def fetch(sample, args):
    path = os.path.join(args.cache, sample + ".json")
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    reads = []
    for chrom, lo, hi, gene in GENES:
        cmd = ["samtools", "view", "-F", "0x904", "-q", str(args.mapq), "-T", args.reference,
               args.url.format(sample=sample), f"{chrom}:{lo}-{hi}"]
        for attempt in range(4):
            run = subprocess.run(cmd, capture_output=True, text=True)
            if run.returncode == 0:
                break
        else:
            print(f"{sample}: read error {run.stderr.strip()[-200:]}", file=sys.stderr)
            return None
        for line in run.stdout.splitlines():
            f = line.split("\t", 6)
            start = int(f[3])
            evs, stop = events(start, f[5])
            reads.append([chrom, start, stop, [[k, s, e, n] for k, s, e, n in evs if n >= args.min_piece]])
    with open(path + ".tmp", "w") as fh:
        json.dump(reads, fh)
    os.replace(path + ".tmp", path)
    return reads


def cluster(all_events, gap):
    """all_events: (kind, chrom, pos, end, length, sample). Greedy clustering by position and length."""
    sites = []
    for kind, chrom, pos, end, length, sample in sorted(all_events, key=lambda x: (x[0], x[1], x[2], x[4])):
        target = None
        for s in reversed(sites):
            if s["type"] != kind or s["chrom"] != chrom or pos - s["last"] > 10 * gap:
                break
            if pos - s["last"] <= gap and abs(length - s["median"]) <= 0.3 * max(length, s["median"]):
                target = s
                break
        if target is None:
            sites.append({"type": kind, "chrom": chrom, "last": pos, "median": length,
                          "members": [(pos, end, length, sample)]})
            continue
        target["members"].append((pos, end, length, sample))
        target["last"] = max(target["last"], pos)
        lens = sorted(m[2] for m in target["members"])
        target["median"] = lens[len(lens) // 2]
    return sites


def genotype(reads, site, args):
    lo, hi = site["pos"] - args.window, site["end"] + args.window
    need = max(50, site["length"] // 2)
    spanning, sizes = 0, []
    for chrom, start, stop, evs in reads:
        if chrom != site["chrom"] or start > lo - args.margin or stop < hi + args.margin:
            continue
        spanning += 1
        size = sum(n for k, s, e, n in evs if k == site["type"] and e >= lo and s <= hi
                   and n <= 3 * site["length"])
        if need <= size <= 2 * site["length"]:
            sizes.append(size)
    if spanning < args.min_spanning:
        return "not evaluable", spanning, len(sizes)
    if len(sizes) >= args.min_reads and len(sizes) / spanning >= args.min_fraction:
        return "present", spanning, len(sizes)
    return "absent", spanning, len(sizes)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--samples", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--cache", required=True, help="directory for the per-sample read summaries")
    p.add_argument("--out", required=True, help="output prefix (.sites.tsv, .carriers.tsv, .json)")
    p.add_argument("--url", default=URL)
    p.add_argument("--min-length", type=int, default=50)
    p.add_argument("--min-piece", type=int, default=30)
    p.add_argument("--cluster", type=int, default=100)
    p.add_argument("--window", type=int, default=150)
    p.add_argument("--margin", type=int, default=100)
    p.add_argument("--mapq", type=int, default=10)
    p.add_argument("--min-reads", type=int, default=2)
    p.add_argument("--min-fraction", type=float, default=0.10)
    p.add_argument("--min-spanning", type=int, default=4)
    p.add_argument("--threads", type=int, default=8)
    args = p.parse_args(argv)
    os.makedirs(args.cache, exist_ok=True)
    with open(args.samples) as fh:
        samples = [l.split()[0] for l in fh if l.strip() and not l.startswith("#")]
    with ThreadPoolExecutor(args.threads) as pool:
        data = dict(zip(samples, pool.map(lambda s: fetch(s, args), samples)))
    failed = sorted(s for s, r in data.items() if r is None)
    data = {s: r for s, r in data.items() if r is not None}

    # candidate events: seen in at least --min-reads reads of one sample at about the same place
    candidates = []
    for sample, reads in data.items():
        per = {}
        for chrom, start, stop, evs in reads:
            for k, s, e, n in evs:
                if n >= args.min_length:
                    per.setdefault((k, chrom), []).append((s, e, n))
        for (k, chrom), lst in per.items():
            lst.sort()
            used = [False] * len(lst)
            for i, (s, e, n) in enumerate(lst):
                if used[i]:
                    continue
                group = [j for j in range(i, len(lst)) if not used[j] and lst[j][0] - s <= args.cluster
                         and abs(lst[j][2] - n) <= 0.3 * max(n, lst[j][2])]
                if len(group) >= args.min_reads:
                    for j in group:
                        used[j] = True
                    g = sorted(lst[j] for j in group)
                    mid = g[len(g) // 2]
                    candidates.append((k, chrom, mid[0], mid[1], sorted(x[2] for x in g)[len(g) // 2], sample))
    sites = []
    for c in cluster(candidates, args.cluster):
        m = sorted(c["members"])
        mid = m[len(m) // 2]
        gene = next(g for ch, lo, hi, g in GENES if ch == c["chrom"])
        lo, hi = next((lo, hi) for ch, lo, hi, g in GENES if ch == c["chrom"])
        if not lo <= mid[0] <= hi:
            continue          # reads reach beyond the gene locus; events there are not part of the test
        sites.append({"gene": gene, "chrom": c["chrom"], "pos": mid[0], "end": mid[1], "type": c["type"],
                      "length": c["median"], "minLength": min(x[2] for x in m), "maxLength": max(x[2] for x in m)})
    sites.sort(key=lambda s: (s["chrom"], s["pos"], s["type"]))
    rows = []
    for i, site in enumerate(sites, 1):
        site["id"] = f"ONT_{site['gene']}_{site['type']}_{site['pos']}"
        calls = {s: genotype(r, site, args) for s, r in data.items()}
        site["carriers"] = sorted(s for s, c in calls.items() if c[0] == "present")
        site["absent"] = sum(c[0] == "absent" for c in calls.values())
        site["notEvaluable"] = sum(c[0] == "not evaluable" for c in calls.values())
        for s in site["carriers"]:
            rows.append((s, site["id"], calls[s][1], calls[s][2]))
    sites = [s for s in sites if s["carriers"]]
    # a long insertion that some reads align in two pieces gives a second, shorter cluster with the
    # same carriers next to it: keep the longer one
    dup = {id(b) for a in sites for b in sites if a is not b and a["type"] == b["type"] == "INS"
           and a["chrom"] == b["chrom"] and abs(a["pos"] - b["pos"]) <= args.cluster and b["length"] < a["length"]
           and set(b["carriers"]) <= set(a["carriers"])}
    sites = [s for s in sites if id(s) not in dup]
    with open(args.out + ".json", "w") as fh:
        json.dump({"source": "1KG-ONT Vienna, hg38 CRAM (s3://1kg-ont-vienna/hg38)", "samples": len(data),
                   "failedSamples": failed, "sites": sites}, fh, indent=1)
    with open(args.out + ".sites.tsv", "w") as fh:
        fh.write("#id\tgene\tchrom\tpos\tend\ttype\tlength\tmin_length\tmax_length\tcarriers\tabsent\tnot_evaluable\n")
        for s in sites:
            fh.write("\t".join(str(x) for x in (s["id"], s["gene"], s["chrom"], s["pos"], s["end"], s["type"], s["length"],
                                                s["minLength"], s["maxLength"], len(s["carriers"]), s["absent"],
                                                s["notEvaluable"])) + "\n")
    with open(args.out + ".carriers.tsv", "w") as fh:
        fh.write("#sample\tsite\tspanning_reads\tsupporting_reads\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"{len(data)} samples, {len(sites)} sites, {len(failed)} failed", file=sys.stderr)


if __name__ == "__main__":
    main()

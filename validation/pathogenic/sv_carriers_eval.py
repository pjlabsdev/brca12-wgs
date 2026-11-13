#!/usr/bin/env python3
"""Are the known 1000 Genomes SVs in BRCA1/BRCA2 reported for carriers and absent in controls?"""
import argparse
import json
import os
import re
import subprocess
import sys

SLOP = 100
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def read_end_signals(report, workdir, sample):
    """Insertion signatures from clipped read ends: from the report, or computed from the BAM for runs
    made before the detector was part of the pipeline."""
    if "insertionSignals" in report:
        return report["insertionSignals"]
    bams = [os.path.join(workdir, sample, "alignment", f"{sample}.markdup.bam")]
    bam = next((b for b in bams if os.path.exists(b)), None)
    if not bam:
        return None
    regions = os.path.join(REPO, "assets", "brca12_regions.bed")
    out = os.path.join(workdir, sample, "insertion_signals.json")
    view = subprocess.Popen(["samtools", "view", "-L", regions, bam], stdout=subprocess.PIPE)
    subprocess.run([sys.executable, os.path.join(REPO, "bin", "brca12_insertion_signals.py"),
                    "--regions", regions, "--out", out], stdin=view.stdout, check=True)
    view.wait()
    return json.load(open(out))


KEEP = ("gene", "chrom", "pos", "end", "targetSiteDuplication", "clippedReadsLeftJunction",
        "clippedReadsRightJunction", "spanningReads", "clippedFraction", "signature")


def reads_at_site(workdir, sample, chrom, pos, slop=15):
    """Raw read counts at an insertion site: reads over the site, and reads clipped (>= 10 bases) at it."""
    bam = os.path.join(workdir, sample, "alignment", f"{sample}.markdup.bam")
    if not os.path.exists(bam):
        return None
    out = subprocess.run(["samtools", "view", "-F", "0xF04", "-q", "20", bam, f"{chrom}:{pos - slop}-{pos + slop}"],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    clipped = 0
    for line in out:
        f = line.split("\t", 6)
        start, ops = int(f[3]), re.findall(r"(\d+)([MIDNSHP=X])", f[5])
        end = start + sum(int(n) for n, op in ops if op in "MDN=X") - 1
        clipped += ((ops[0][1] == "S" and int(ops[0][0]) >= 10 and abs(start - pos) <= slop)
                    or (ops[-1][1] == "S" and int(ops[-1][0]) >= 10 and abs(end - pos) <= slop))
    return {"reads": len(out), "clippedAtSite": clipped}


def table(path):
    rows = []
    for line in open(path):
        if line.startswith("# "):
            continue
        if line.startswith("#"):
            keys = line[1:].rstrip("\n").split("\t")
        elif line.strip():
            rows.append(dict(zip(keys, line.rstrip("\n").split("\t"))))
    return rows


def look(site, report, calls_vcf):
    """Evidence for the site in the report (Manta) and in the small-variant calls."""
    pos, end, kind, length = int(site["pos"]), int(site["end"]), site["type"], int(site["length"])
    found = {"manta": [], "smallVariantCaller": []}
    for sv in report.get("svs", []):
        if sv["chrom"] != site["chrom"]:
            continue
        sv_end = sv["end"] or sv["pos"]
        overlap = min(sv_end, end) - max(sv["pos"], pos)
        # panel coordinates of merged short-read SV calls are not base-exact: 50 % reciprocal overlap
        if kind == "DEL" and sv["type"] == "DEL" and overlap >= 0.5 * (end - pos) and overlap >= 0.5 * (sv_end - sv["pos"]):
            found["manta"].append({**{k: sv[k] for k in ("type", "pos", "end", "length", "genotype", "filter", "flags")},
                                   "startDifference": sv["pos"] - pos, "endDifference": sv_end - end})
        if kind == "ALU" and sv["type"] in ("INS", "BND") and abs(sv["pos"] - pos) <= SLOP:
            found["manta"].append({k: sv[k] for k in ("type", "pos", "length", "genotype", "filter", "flags")})
    out = subprocess.run(["bcftools", "query", "-r", f"{site['chrom']}:{pos - SLOP}-{end + SLOP}",
                          "-f", "%POS\t%REF\t%ALT[\t%GT]\n", calls_vcf], capture_output=True, text=True, check=True)
    for line in out.stdout.splitlines():
        p, ref, alts, gt = line.split("\t")
        for alt in alts.split(","):
            delta = len(alt) - len(ref)
            if (kind == "DEL" and abs(-delta - length) <= 0.1 * length) or (kind == "ALU" and delta >= 50):
                found["smallVariantCaller"].append({"pos": int(p), "refLength": len(ref), "altLength": len(alt), "genotype": gt})
    return found


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites", required=True)
    p.add_argument("--carriers", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--alu-genotypes", help="panel genotypes at the Alu sites for all samples (sample, one column per site)")
    p.add_argument("--more-workdir", action="append", default=[],
                   help="further pipeline outputs (sample folders) evaluated at the Alu sites only")
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)
    sites = {s["id"]: s for s in table(args.sites)}
    rows = []
    alu = {r["sample"]: r for r in table(args.alu_genotypes)} if args.alu_genotypes else {}
    samples = [(c, args.workdir) for c in table(args.carriers)]
    known = {c["sample"] for c, _ in samples}
    for more in args.more_workdir:
        samples += [({"sample": name, "population": "", "sv": None, "genotype_in_panel": None}, more)
                    for name in sorted(alu) if name not in known and os.path.isdir(os.path.join(more, name))]
    for c, workdir in samples:
        s = c["sample"]
        path = os.path.join(workdir, s, "brca12_results.json")
        if not os.path.exists(path):
            continue
        report = json.load(open(path))
        calls = os.path.join(workdir, s, "calls", f"{s}.small.vcf.gz")
        signals = read_end_signals(report, workdir, s) or []
        for sid, site in sites.items():
            expected = c["genotype_in_panel"] if c["sv"] == sid else ("0|0" if c["sv"] == "CONTROL" else None)
            if expected is None and alu.get(s, {}).get(sid, ".") != ".":
                expected = alu[s][sid]
            if expected is None:
                continue                                  # genotype at the other sites not looked up
            found = look(site, report, calls)
            found["readEndSignal"] = [{k: x[k] for k in KEEP} for x in signals if site["type"] == "ALU"
                                      and x["chrom"] == site["chrom"] and abs(x["pos"] - int(site["pos"])) <= SLOP]
            if site["type"] == "ALU":
                found["readsAtSite"] = reads_at_site(workdir, s, site["chrom"], int(site["pos"]))
            rows.append({"sample": s, "population": c["population"], "site": sid, "type": site["type"],
                         "length": int(site["length"]), "panelGenotype": expected,
                         "reported": bool(found["manta"] or found["smallVariantCaller"]), **found})
    summary = {"rows": rows, "bySite": {}}
    for r in rows:
        carrier = r["panelGenotype"] != "0|0"
        t = summary["bySite"].setdefault(r["site"], {"carriers": 0, "carriersReported": 0, "carriersWithReadEndSignal": 0,
                                                     "controls": 0, "controlsReported": 0, "controlsWithReadEndSignal": 0})
        t["carriers" if carrier else "controls"] += 1
        t["carriersReported" if carrier else "controlsReported"] += r["reported"]
        t["carriersWithReadEndSignal" if carrier else "controlsWithReadEndSignal"] += bool(r["readEndSignal"])
    json.dump(summary, open(args.out + ".json", "w"), indent=2, ensure_ascii=False)
    with open(args.out + ".tsv", "w") as handle:
        handle.write("sample\tsite\ttype\tpanel_genotype\treported\tmanta\tsmall_variant_caller\tread_end_signal\n")
        for r in rows:
            line = "\t".join([r["sample"], r["site"], r["type"], r["panelGenotype"], str(r["reported"]),
                              ";".join(f"{m['type']}@{m['pos']}:{m['genotype']}:{m['filter']}" for m in r["manta"]) or "-",
                              ";".join(f"{m['pos']}:{m['altLength'] - m['refLength']:+d}:{m['genotype']}" for m in r["smallVariantCaller"]) or "-",
                              ";".join(f"{m['pos']}:{m['clippedReadsLeftJunction']}+{m['clippedReadsRightJunction']}:{m['clippedFraction']}"
                                       for m in r["readEndSignal"]) or "-"])
            handle.write(line + "\n")
            if r["panelGenotype"] != "0|0" or r["reported"] or r["readEndSignal"]:
                print(line)
    print(json.dumps(summary["bySite"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Are the known 1000 Genomes SVs in BRCA1/BRCA2 reported for carriers and absent in controls?"""
import argparse
import json
import os
import subprocess
import sys

SLOP = 100


def table(path):
    rows = []
    for line in open(path):
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
        if kind == "DEL" and sv["type"] == "DEL" and abs(sv["pos"] - pos) <= SLOP and abs(sv_end - end) <= SLOP:
            found["manta"].append({k: sv[k] for k in ("type", "pos", "end", "length", "genotype", "filter", "flags")})
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
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)
    sites = {s["id"]: s for s in table(args.sites)}
    rows = []
    for c in table(args.carriers):
        s = c["sample"]
        path = os.path.join(args.workdir, s, "brca12_results.json")
        if not os.path.exists(path):
            continue
        report = json.load(open(path))
        calls = os.path.join(args.workdir, s, "calls", f"{s}.small.vcf.gz")
        for sid, site in sites.items():
            expected = c["genotype_in_panel"] if c["sv"] == sid else ("0|0" if c["sv"] == "CONTROL" else None)
            if expected is None:
                continue                                  # genotype at the other sites not looked up
            found = look(site, report, calls)
            rows.append({"sample": s, "population": c["population"], "site": sid, "type": site["type"],
                         "length": int(site["length"]), "panelGenotype": expected,
                         "reported": bool(found["manta"] or found["smallVariantCaller"]), **found})
    summary = {"rows": rows, "bySite": {}}
    for r in rows:
        carrier = r["panelGenotype"] != "0|0"
        t = summary["bySite"].setdefault(r["site"], {"carriers": 0, "carriersReported": 0, "controls": 0, "controlsReported": 0})
        t["carriers" if carrier else "controls"] += 1
        t["carriersReported" if carrier else "controlsReported"] += r["reported"]
    json.dump(summary, open(args.out + ".json", "w"), indent=2, ensure_ascii=False)
    with open(args.out + ".tsv", "w") as handle:
        handle.write("sample\tsite\ttype\tpanel_genotype\treported\tmanta\tsmall_variant_caller\n")
        for r in rows:
            line = "\t".join([r["sample"], r["site"], r["type"], r["panelGenotype"], str(r["reported"]),
                              ";".join(f"{m['type']}@{m['pos']}:{m['genotype']}:{m['filter']}" for m in r["manta"]) or "-",
                              ";".join(f"{m['pos']}:{m['altLength'] - m['refLength']:+d}:{m['genotype']}" for m in r["smallVariantCaller"]) or "-"])
            handle.write(line + "\n")
            if r["panelGenotype"] != "0|0" or r["reported"]:
                print(line)
    print(json.dumps(summary["bySite"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

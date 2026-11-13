#!/usr/bin/env python3
"""Was the known pathogenic variant of each 1000 Genomes carrier called? Plus concordance totals."""
import argparse
import json
import os
import subprocess
import sys


def normalised(vcf, reference, region):
    cmd = f"bcftools norm -m -any -f {reference} -r {region} {vcf} 2>/dev/null | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER[\\t%GT\\t%DP\\t%GQ\\t%AD]\\n'"
    out = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, check=True)
    return [line.split("\t") for line in out.stdout.splitlines() if line.strip()]


def norm_one(chrom, pos, ref, alt, reference):
    vcf = f"##fileformat=VCFv4.2\n##contig=<ID={chrom}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n"
    out = subprocess.run(["bcftools", "norm", "-f", reference, "-"], input=vcf, capture_output=True, text=True, check=True)
    f = [l for l in out.stdout.splitlines() if not l.startswith("#")][0].split("\t")
    return f[0], int(f[1]), f[3], f[4]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)
    rows, totals = [], {}
    for line in open(args.table):
        if line.startswith("#") or not line.strip():
            continue
        s, run, prefix, pop, gene, chrom, pos, ref, alt, hgvs, csq, clinvar, nygc_filter = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt = norm_one(chrom, pos, ref, alt, args.reference)
        calls = normalised(os.path.join(args.workdir, s, "calls", f"{s}.small.vcf.gz"), args.reference,
                           f"{chrom}:{pos - 60}-{pos + 60}")
        hit = next((c for c in calls if int(c[1]) == pos and c[2] == ref and c[3] == alt), None)
        row = {"sample": s, "population": pop, "gene": gene, "variant": f"{chrom}:{pos} {ref}>{alt}", "hgvs": hgvs,
               "consequence": csq, "clinvarIn1000gAnnotation": clinvar, "nygcFilter": nygc_filter,
               "detected": bool(hit), "genotype": hit[6] if hit else None,
               "depth": int(hit[7]) if hit and hit[7].isdigit() else None,
               "GQ": int(hit[8]) if hit and hit[8].isdigit() else None, "AD": hit[9] if hit else None}
        report = json.load(open(os.path.join(args.workdir, s, "brca12_results.json")))
        row["inReport"] = str(pos) in json.dumps(report.get("genes", {}))
        bench = json.load(open(os.path.join(args.workdir, s, "benchmark", "summary.json")))
        row["concordance"] = bench["results"]
        for r in bench["results"]:
            t = totals.setdefault((r["region"], r["type"]), {"nygc": 0, "TP": 0, "FP": 0, "FN": 0})
            t["nygc"] += r["truthTotal"]; t["TP"] += r["TP"]; t["FP"] += r["FP"]; t["FN"] += r["FN"]
        rows.append(row)
    summary = {"carriers": rows, "detected": sum(r["detected"] for r in rows), "total": len(rows),
               "concordanceWithNygc": [{"region": k[0], "type": k[1], **v} for k, v in totals.items()]}
    json.dump(summary, open(args.out + ".json", "w"), indent=2)
    with open(args.out + ".tsv", "w") as handle:
        keys = ["sample", "population", "gene", "hgvs", "consequence", "variant", "detected", "genotype", "depth", "GQ", "AD", "nygcFilter"]
        handle.write("\t".join(keys) + "\n")
        for r in rows:
            handle.write("\t".join(str(r[k]) for k in keys) + "\n")
            print("\t".join(str(r[k]) for k in keys))
    for t in summary["concordanceWithNygc"]:
        print(t)
    return 0


if __name__ == "__main__":
    sys.exit(main())

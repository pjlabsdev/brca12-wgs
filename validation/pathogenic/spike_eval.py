#!/usr/bin/env python3
"""Evaluate the spike-in runs: small variants by type and size, large events per caller."""
import argparse
import gzip
import json
import os
import subprocess
import sys

SIZE_CLASSES = [(1, 1), (2, 5), (6, 20), (21, 50)]
SLOP = 200


def vcf_positions(path):
    out = set()
    if not os.path.exists(path):
        return out
    for line in gzip.open(path, "rt"):
        if not line.startswith("#"):
            f = line.split("\t")
            out.add((f[0], int(f[1])))
    return out


def size_class(kind, size):
    if kind == "SNV":
        return "SNV"
    lo, hi = next((a, b) for a, b in SIZE_CLASSES if a <= size <= b)
    name = "deletion" if kind == "DEL" else "insertion/duplication"
    return f"{name} {lo} bp" if lo == hi else f"{name} {lo}-{hi} bp"


def small_calls_near(vcf, chrom, start, end):
    out = subprocess.run(["bcftools", "query", "-r", f"{chrom}:{max(1, start - SLOP)}-{end + SLOP}",
                          "-f", "%POS\t%REF\t%ALT[\t%GT]\n", vcf], capture_output=True, text=True, check=True)
    return [l.split("\t") for l in out.stdout.splitlines() if l.strip()]


def evaluate_large(event, report, calls_vcf):
    gene, kind, start, end, size = event["gene"], event["kind"], event["start"], event["end"], event["size"]
    want_sv = {"LDEL": ("DEL",), "LDUP": ("DUP",), "ALU": ("INS", "BND")}[kind]
    want_cnv = {"LDEL": "Deletion", "LDUP": "Duplikation", "ALU": None}[kind]
    res = {"manta": None, "cnvkit": None, "smallVariantCaller": None, "exonDosageHints": [], "notes": []}
    for sv in report.get("svs", []):
        if sv["gene"] != gene or sv["type"] not in want_sv:
            continue
        sv_end = sv["end"] or sv["pos"]
        ok = abs(sv["pos"] - start) <= SLOP and (kind == "ALU" or abs(sv_end - end) <= SLOP)
        if ok:
            res["manta"] = {k: sv[k] for k in ("type", "pos", "end", "length", "genotype", "filter", "pairedReads", "splitReads", "affectedExons")}
            break
    for cnv in report.get("cnvs", []):
        if want_cnv and cnv["gene"] == gene and cnv["type"] == want_cnv and cnv["start"] <= end and cnv["end"] >= start:
            res["cnvkit"] = {k: cnv[k] for k in ("start", "end", "log2", "copyNumber", "type", "probes", "affectedExons")}
            break
    if kind in ("LDEL", "ALU") and size <= 400:
        for pos, ref, alt, gt in small_calls_near(calls_vcf, event["chrom"], start, end):
            for a in alt.split(","):
                delta = len(a) - len(ref)
                if (kind == "LDEL" and abs(-delta - size) <= 0.1 * size) or (kind == "ALU" and delta >= 50):
                    res["smallVariantCaller"] = {"pos": int(pos), "refLength": len(ref), "altLength": len(a), "genotype": gt}
    want_hint = "möglicher Verlust" if kind == "LDEL" else "möglicher Zugewinn" if kind == "LDUP" else None
    res["exonDosageHints"] = [h["exon"] for h in report.get("exonDosageHints", []) if h["gene"] == gene and h["type"] == want_hint]
    res["affectedExonDepthRatios"] = {x["label"]: x["dosageRatio"] for x in report.get("exonCoverage", [])
                                      if x["gene"] == gene and x["start"] <= end and x["end"] >= start}
    res["affectedExonStatus"] = sorted({x["status"] for x in report.get("exonCoverage", [])
                                        if x["gene"] == gene and x["start"] <= end and x["end"] >= start})
    res["needsAttention"] = report["summary"].get("needsAttention")
    res["detected"] = bool(res["manta"] or res["cnvkit"] or res["smallVariantCaller"])
    return res


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--simdir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--small", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)
    small_events = {}
    for line in open(args.small):
        if not line.startswith("#") and line.strip():
            f = line.rstrip("\n").split("\t")
            small_events.setdefault(f[0], []).append({"gene": f[1], "kind": f[2], "chrom": f[3], "pos": int(f[4]), "ref": f[5],
                                                      "alt": f[6], "size": int(f[7]), "context": f[8], "note": f[9]})
    samples = sorted(d for d in os.listdir(args.outdir) if os.path.exists(os.path.join(args.outdir, d, "brca12_results.json")))
    summary = {"samples": [], "smallByClass": {}, "smallMissed": [], "largeEvents": [], "unexpected": []}
    for s in samples:
        report = json.load(open(os.path.join(args.outdir, s, "brca12_results.json")))
        truth = json.load(open(os.path.join(args.simdir, s, f"{s}.truth.json")))
        bench = json.load(open(os.path.join(args.outdir, s, "benchmark", "summary.json")))
        run = os.path.join(args.outdir, s, "benchmark", "vcfeval_gene_loci")
        tp, fn = vcf_positions(f"{run}/tp-baseline.vcf.gz"), vcf_positions(f"{run}/fn.vcf.gz")
        for e in small_events.get(s, []):
            key = (e["chrom"], e["pos"])
            found = key in tp and key not in fn
            cls = summary["smallByClass"].setdefault(size_class(e["kind"], e["size"]), {"spiked": 0, "found": 0})
            cls["spiked"] += 1
            cls["found"] += found
            if not found:
                summary["smallMissed"].append({"sample": s, **e})
        fp = sum(r["FP"] for r in bench["results"] if r["region"] == "gene_loci")
        calls_vcf = os.path.join(args.outdir, s, "calls", f"{s}.small.vcf.gz")
        matched_sv, matched_cnv = set(), set()
        for event in truth["largeEvents"]:
            res = evaluate_large(event, report, calls_vcf)
            summary["largeEvents"].append({**event, **res})
            if res["manta"]:
                matched_sv.add((event["gene"], res["manta"]["pos"]))
            if res["cnvkit"]:
                matched_cnv.add((event["gene"], res["cnvkit"]["start"]))
        truth_genes = {e["gene"] for e in truth["largeEvents"]}
        for sv in report.get("svs", []):
            if (sv["gene"], sv["pos"]) not in matched_sv:
                summary["unexpected"].append({"sample": s, "caller": "Manta", "inGeneWithSpikedEvent": sv["gene"] in truth_genes,
                                              **{k: sv[k] for k in ("gene", "type", "pos", "end", "length", "filter")}})
        for cnv in report.get("cnvs", []):
            if cnv["type"] != "kein Hinweis" and (cnv["gene"], cnv["start"]) not in matched_cnv:
                summary["unexpected"].append({"sample": s, "caller": "CNVkit", "inGeneWithSpikedEvent": cnv["gene"] in truth_genes,
                                              **{k: cnv[k] for k in ("gene", "type", "start", "end", "log2", "copyNumber")}})
        summary["samples"].append({"sample": s, "simulation": truth["stats"], "smallSpiked": truth["smallVariants"],
                                   "smallFalsePositivesGeneLoci": fp, "largeSpiked": len(truth["largeEvents"]),
                                   "svCount": len(report.get("svs", [])),
                                   "cnvCandidates": sum(1 for c in report.get("cnvs", []) if c["type"] != "kein Hinweis"),
                                   "exonDosageHintCount": len(report.get("exonDosageHints", [])),
                                   "uncertainExons": report["summary"]["uncertainExonCount"],
                                   "toolVersions": report.get("toolVersions")})
    json.dump(summary, open(args.out + ".json", "w"), indent=2, ensure_ascii=False)
    print("class\tspiked\tfound")
    for k, v in summary["smallByClass"].items():
        print(f"{k}\t{v['spiked']}\t{v['found']}")
    for e in summary["largeEvents"]:
        print(e["sample"], e["gene"], e["kind"], e["size"], "manta" if e["manta"] else "-", "cnvkit" if e["cnvkit"] else "-",
              "small" if e["smallVariantCaller"] else "-", "hints:" + ",".join(e["exonDosageHints"]), "ratios:" + json.dumps(e["affectedExonDepthRatios"], ensure_ascii=False), sep="\t")
    print("unexpected:", len(summary["unexpected"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())

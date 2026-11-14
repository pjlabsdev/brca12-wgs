#!/usr/bin/env python3
"""Evaluate the mosaic spike-ins (run_mosaic.sh): per allele fraction, how many spiked events are
reported as mosaic candidates (Mutect2), as germline calls (HaplotypeCaller), or by either, and how
many reported variants were not spiked (false positives; the background's own variants are not in
the simulated windows).

A call matches a spiked event on the same chromosome within 10 bp with the same length change
(and the same alternative base for SNVs), so that differently placed indels in repeats still match.
"""
import argparse
import gzip
import json
import os
import statistics


def read_vcf(path):
    opener = gzip.open if path.endswith(".gz") else open
    rows = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            for alt in f[4].split(","):
                if alt in ("*", "<NON_REF>"):
                    continue
                gt = f[9].split(":")[0] if len(f) > 9 else ""
                if gt in ("0/0", "0|0", "./.", "."):
                    continue
                rows.append({"chrom": f[0], "pos": int(f[1]), "ref": f[3], "alt": alt, "filter": f[6]})
    return rows


def same(a, b):
    if a["chrom"] != b["chrom"] or abs(a["pos"] - b["pos"]) > 10:
        return False
    if len(a["ref"]) - len(a["alt"]) != len(b["ref"]) - len(b["alt"]):
        return False
    if len(a["ref"]) == len(a["alt"]) == 1:
        return a["pos"] == b["pos"] and a["alt"] == b["alt"]
    return True


def evaluate(workdir, sample):
    truth_path = os.path.join(workdir, "sim", sample, f"{sample}.truth.vcf")
    truth = read_vcf(truth_path)
    meta = json.load(open(os.path.join(workdir, "sim", sample, f"{sample}.truth.json")))
    result = json.load(open(os.path.join(workdir, "out", sample, "brca12_results.json")))
    germline = [v for g in result["genes"].values() for v in g["smallVariants"]]
    mosaic = result.get("mosaicCandidates", [])
    found = {"mosaic": 0, "germline": 0, "either": 0}
    fractions = []
    for t in truth:
        m = next((v for v in mosaic if same(t, v)), None)
        g = any(same(t, v) for v in germline)
        found["mosaic"] += m is not None
        found["germline"] += g
        found["either"] += bool(m) or g
        if m is not None and m.get("alleleFraction") is not None:
            fractions.append(m["alleleFraction"])
    fp_mosaic = [v for v in mosaic if not any(same(t, v) for t in truth)]
    fp_germline = [v for v in germline if not any(same(t, v) for t in truth)]
    return {
        "sample": sample, "spikedAlleleFraction": None if not truth else meta.get("alleleFraction"),
        "spikedEvents": len(truth),
        "reportedAsMosaicCandidate": found["mosaic"], "calledByGermlineCaller": found["germline"],
        "reportedEither": found["either"],
        "medianObservedFraction": round(statistics.median(fractions), 3) if fractions else None,
        "falsePositiveMosaicCandidates": len(fp_mosaic),
        "falsePositiveGermlineCalls": len(fp_germline),
        "falsePositiveExamples": [f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']} AF {v.get('alleleFraction')}"
                                  for v in fp_mosaic[:5]],
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--workdir", required=True)
    p.add_argument("--out", required=True, help="output prefix (.json, .tsv)")
    args = p.parse_args(argv)
    samples = sorted(d for d in os.listdir(os.path.join(args.workdir, "out"))
                     if os.path.exists(os.path.join(args.workdir, "out", d, "brca12_results.json")))
    rows = [evaluate(args.workdir, s) for s in samples]
    json.dump(rows, open(args.out + ".json", "w"), indent=1)
    keys = ["sample", "spikedAlleleFraction", "spikedEvents", "reportedAsMosaicCandidate", "calledByGermlineCaller",
            "reportedEither", "medianObservedFraction", "falsePositiveMosaicCandidates", "falsePositiveGermlineCalls"]
    with open(args.out + ".tsv", "w") as handle:
        handle.write("\t".join(keys) + "\n")
        for r in rows:
            handle.write("\t".join("-" if r[k] is None else str(r[k]) for k in keys) + "\n")
    print(open(args.out + ".tsv").read())


if __name__ == "__main__":
    main()

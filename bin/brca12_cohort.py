#!/usr/bin/env python3
"""Cohort summary over several brca12-wgs results (one run or batch of samples).

    python3 bin/brca12_cohort.py --results 'runs/*/brca12_results.json' --run-name RUN01 --outdir cohort/

Writes cohort_samples.tsv, cohort_variants.tsv, cohort_contamination.tsv and cohort.xlsx.

Contamination hint (research use only): a small variant called in sample A with a low allele
fraction (default 3 to 25 %) while it is clearly heterozygous or homozygous in another sample B of
the same batch can mean contamination or carry-over from B. The hint lists B; it does not prove
contamination and needs checking (for example identity SNPs or a repeated library).
"""
import argparse
import csv
import glob
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from brca12lib import tables, xlsx  # noqa: E402

RUO = ("Nur Forschung (RUO). Keine medizinische Interpretation ohne klinisch validierte Bestätigung. "
       "Research use only. No medical interpretation without clinically validated confirmation.")


def load(patterns):
    paths = []
    for pattern in patterns:
        found = sorted(glob.glob(pattern)) or ([pattern] if os.path.exists(pattern) else [])
        for path in found:
            if os.path.isdir(path):
                path = os.path.join(path, "brca12_results.json")
            if os.path.exists(path) and path not in paths:
                paths.append(path)
    results = []
    for path in paths:
        with open(path) as handle:
            result = json.load(handle)
        result["_path"] = path
        results.append(result)
    names = [r.get("sample") for r in results]
    if len(set(names)) != len(names):
        raise SystemExit("[brca12] duplicate sample names in the results: " + ", ".join(sorted({n for n in names if names.count(n) > 1})))
    return results


def variant_key(row):
    # SMALL_HEADER: sample gene chrom pos_GRCh38 pos_GRCh37 ref alt ...
    return (row[2], row[3], row[5], row[6])


def contamination(rows, low=(0.03, 0.25), clear_het=(0.30, 0.70), clear_hom=0.85):
    by_key = {}
    for row in rows:
        by_key.setdefault(variant_key(row), []).append(row)
    hints = []
    for key, carriers in by_key.items():
        clear = [r for r in carriers if r[9] is not None and (clear_het[0] <= r[9] <= clear_het[1] or r[9] >= clear_hom)]
        for r in carriers:
            fraction = r[9]
            if fraction is None or not low[0] <= fraction <= low[1]:
                continue
            sources = [c[0] for c in clear if c[0] != r[0]]
            if sources:
                hints.append({"sample": r[0], "gene": r[1], "chrom": key[0], "pos": key[1], "ref": key[2],
                              "alt": key[3], "hgvsc": r[13], "allele_fraction": fraction,
                              "possible_source": ", ".join(sources)})
    return hints


def coverage_summary(result):
    depths = [e.get("meanDepth") for e in result.get("exonCoverage") or [] if e.get("meanDepth") is not None]
    uncertain = [e.get("label") for e in result.get("exonCoverage") or [] if e.get("status") != "ok"]
    return (round(statistics.median(depths), 1) if depths else None), uncertain


def build(results, run_name):
    variant_rows = [row for r in results for row in tables.small_rows(r)]
    counts = {}
    for row in variant_rows:
        counts[variant_key(row)] = counts.get(variant_key(row), 0) + 1
    hints = contamination(variant_rows)
    hint_keys = {(h["sample"], h["chrom"], h["pos"], h["ref"], h["alt"]): h for h in hints}
    n = len(results)
    variants = [["run"] + tables.SMALL_HEADER + ["cohort_carriers", "cohort_fraction", "contamination_hint"]]
    for row in variant_rows:
        key = variant_key(row)
        hint = hint_keys.get((row[0],) + key)
        variants.append([run_name] + row + [counts[key], round(counts[key] / n, 3) if n else None,
                                            f"niedriger Allelanteil, in {hint['possible_source']} klar vorhanden" if hint else None])
    per_sample_hints = {}
    for h in hints:
        per_sample_hints[h["sample"]] = per_sample_hints.get(h["sample"], 0) + 1
    samples = [["run", "sample", "small_variants", "clinvar_pathogenic", "lof_candidates", "cnv_candidates",
                "exon_dosage_hints", "svs", "uncertain_exons", "median_exon_depth", "low_coverage_intervals",
                "contamination_hints", "needs_attention", "result_file"]]
    for r in results:
        s = r.get("summary") or {}
        median, uncertain = coverage_summary(r)
        samples.append([run_name, r.get("sample"), s.get("smallVariantCount"), s.get("clinvarPathogenicCount"),
                        s.get("lossOfFunctionCandidateCount"), s.get("cnvCandidateCount"), s.get("exonDosageHintCount"),
                        s.get("svCount"), "; ".join(uncertain), median, s.get("lowCoverageIntervalCount"),
                        per_sample_hints.get(r.get("sample"), 0), s.get("needsAttention"), r["_path"]])
    contam = [["run", "sample", "gene", "chrom", "pos", "ref", "alt", "hgvsc", "allele_fraction", "possible_source"]]
    contam += [[run_name, h["sample"], h["gene"], h["chrom"], h["pos"], h["ref"], h["alt"], h["hgvsc"],
                h["allele_fraction"], h["possible_source"]] for h in hints]
    return samples, variants, contam


def write_tsv(path, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in rows:
            writer.writerow(["" if v is None else v for v in row])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", nargs="+", required=True, help="brca12_results.json files, directories or glob patterns")
    parser.add_argument("--run-name", default="BRCA12")
    parser.add_argument("--outdir", default="cohort")
    args = parser.parse_args(argv)
    results = load(args.results)
    if not results:
        raise SystemExit("[brca12] no brca12_results.json found")
    os.makedirs(args.outdir, exist_ok=True)
    samples, variants, contam = build(results, args.run_name)
    write_tsv(os.path.join(args.outdir, "cohort_samples.tsv"), samples)
    write_tsv(os.path.join(args.outdir, "cohort_variants.tsv"), variants)
    write_tsv(os.path.join(args.outdir, "cohort_contamination.tsv"), contam)
    xlsx.write(os.path.join(args.outdir, "cohort.xlsx"),
               [("Proben", samples), ("Varianten", variants), ("Kontaminationshinweise", contam),
                ("Hinweis RUO", [["Hinweis"], [RUO]])])
    print(f"[brca12] {len(results)} samples, {len(variants) - 1} variant calls, "
          f"{len(contam) - 1} contamination hints -> {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

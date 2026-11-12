#!/usr/bin/env python3
"""Fasst alle BRCA1/BRCA2-Ergebnisse zu brca12_results.json, brca12_report.md und .html zusammen.

The JSON file (schema genome.brca12.v1) is the machine-readable interface for downstream reports.
"""
import argparse
import datetime
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from brca12lib import coverage, render, structural, variants  # noqa: E402

SCHEMA = "genome.brca12.v1"

GENES = {
    "BRCA1": {"strand": "-", "transcript": "ENST00000357654 / NM_007294.4 (MANE Select)"},
    "BRCA2": {"strand": "+", "transcript": "ENST00000380152 / NM_000059.4 (MANE Select)"},
}

LIMITATIONS = [
    "Nur GRCh38/hg38. Koordinaten BRCA1 chr17:43044295-43170327 (Minus-Strang), BRCA2 chr13:32315077-32400268 (Plus-Strang).",
    "Kurz-Read-WGS erfasst Mosaike unter etwa 15 % Allelanteil, komplexe Umlagerungen und Bereiche mit Pseudogen-Homologie (BRCA1-Exons 1-2 und BRCA1P1) nur eingeschränkt.",
    "CNVs stammen aus einer Einzelprobe mit flacher CNVkit-Referenz und gelten nur als Kandidaten. Einzelexon-Deletionen oder -Duplikationen unter etwa 1 kb sind unsicher; MLPA ist die übliche Bestätigung.",
    "Manta meldet nur SVs mit Bruchpunkt im Suchfenster um BRCA1/BRCA2. Insertionen mobiler Elemente (z. B. Alu in BRCA2-Exon 3) werden nicht sicher erkannt.",
    "ACMG/AMP: nur maschinelle Hinweise auf einzelne Kriterien (ENIGMA-Schwellen). Keine Klassifikation; PS3, PM5, PS1, PP1, PS4 und Spleißanalysen fehlen.",
    "ClinVar-Angaben hängen von der Version der verwendeten ClinVar-VCF ab.",
]


LIMITATIONS_EN = [
    "GRCh38/hg38 only. Coordinates BRCA1 chr17:43044295-43170327 (minus strand), BRCA2 chr13:32315077-32400268 (plus strand).",
    "Short-read WGS has limited sensitivity for mosaics below about 15% allele fraction, complex rearrangements and regions with pseudogene homology (BRCA1 exons 1-2 and BRCA1P1).",
    "CNVs come from a single sample with a flat CNVkit reference and are candidates only. Single-exon deletions or duplications below about 1 kb are uncertain; MLPA is the usual confirmation.",
    "Manta only reports SVs with a breakpoint in the search window around BRCA1/BRCA2. Mobile element insertions (e.g. Alu in BRCA2 exon 3) are not reliably detected.",
    "ACMG/AMP: automated hints on individual criteria only (ENIGMA thresholds). No classification; PS3, PM5, PS1, PP1, PS4 and splicing analyses are missing.",
    "ClinVar entries depend on the version of the ClinVar VCF used.",
]


def read_versions(path):
    versions = {}
    with open(path) as handle:
        for line in handle:
            if "\t" in line:
                key, value = line.rstrip("\n").split("\t", 1)
                versions[key] = value.strip() or None
    return versions


def read_exact_regions(path):
    regions = []
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 4:
                regions.append((fields[0], int(fields[1]), int(fields[2]), fields[3]))
    return regions


def build(args):
    with open(args.reference_check) as handle:
        check = json.load(handle)
    regions = read_exact_regions(args.regions)
    medians = coverage.read_chrom_median(args.chrom_median)
    exons = coverage.exon_table(args.cov_q20, args.cov_thresholds, args.cov_q0, medians,
                                args.min_exon_depth, args.min_exon_pct20, args.max_mapq0)
    uncertain = [e for e in exons if e["status"] != "ok"]
    small = dict(item.split("=", 1) for item in args.small)
    genes = {}
    for gene, meta in GENES.items():
        chrom, start, end, _ = next(r for r in regions if r[3] == gene)
        rows = variants.small_variants(small[gene], gene, uncertain) if gene in small else []
        genes[gene] = {"symbol": gene, "chrom": chrom, "start": start + 1, "end": end,
                       "strand": meta["strand"], "transcript": meta["transcript"],
                       "coverage": coverage.gene_summary(exons, gene), "smallVariants": rows}
    low = coverage.low_coverage_intervals(args.cov_perbase, regions, args.low_depth)
    cnvs = [c for c in structural.cnv_segments(args.cnv, args.cnr, regions, exons)]
    svs = structural.structural_variants(args.sv, regions, exons)
    versions = read_versions(args.versions)
    versions.setdefault("nextflow", args.nextflow_version)
    invocation = " ".join(["nextflow run main.nf", f"--sample {args.sample}",
                           f"--reference {args.reference}",
                           (f"--bam {args.inputs}" if args.input_mode == "bam"
                            else "--reads_r1 {} --reads_r2 {}".format(*args.inputs.split(",", 1))),
                           f"--snv_caller {args.snv_caller}", "-profile <conda|standard>"])
    result = {
        "schema": SCHEMA,
        "sample": args.sample,
        "generatedAt": datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat(),
        "assembly": "GRCh38",
        "researchUseOnly": True,
        "disclaimer": {"de": render.DISCLAIMER_DE, "en": render.DISCLAIMER_EN},
        "clinicalClassification": None,
        "input": {"mode": args.input_mode, "files": [os.path.basename(p) for p in args.inputs.split(",")]},
        "reference": {"path": args.reference, "assembly": check.get("assembly"), "check": check},
        "snvCaller": "GATK HaplotypeCaller" if args.snv_caller == "haplotypecaller" else "DeepVariant",
        "thresholds": {"minExonDepth": args.min_exon_depth, "minExonPct20x": args.min_exon_pct20,
                       "lowDepth": args.low_depth, "maxLowMapqFraction": args.max_mapq0},
        "genes": genes,
        "exonCoverage": exons,
        "lowCoverageIntervals": low,
        "cnvs": [c for c in cnvs if c["type"] != "kein Hinweis"],
        "cnvSegmentsChecked": len(cnvs),
        "exonDosageHints": structural.exon_dosage_hints(exons),
        "svs": svs,
        "toolVersions": versions,
        "invocation": invocation,
        "nextflowSession": args.session_id,
        "limitations": LIMITATIONS,
        "limitationsEN": LIMITATIONS_EN,
    }
    result["summary"] = summarize(result)
    return result


def summarize(result):
    small = [v for g in result["genes"].values() for v in g["smallVariants"]]
    lof = [v for v in small if any(h["code"].startswith("PVS1") for h in v["acmgHints"])]
    clinvar_plp = [v for v in small if v["clinvar"]["significance"]
                   and "athogenic" in v["clinvar"]["significance"]
                   and "onflicting" not in v["clinvar"]["significance"]]
    return {
        "smallVariantCount": len(small),
        "clinvarPathogenicCount": len(clinvar_plp),
        "lossOfFunctionCandidateCount": len(lof),
        "cnvCandidateCount": len(result["cnvs"]),
        "exonDosageHintCount": len(result["exonDosageHints"]),
        "svCount": len(result["svs"]),
        "uncertainExonCount": sum(1 for e in result["exonCoverage"] if e["status"] != "ok"),
        "lowCoverageIntervalCount": len(result["lowCoverageIntervals"]),
        "needsAttention": bool(clinvar_plp or lof or result["cnvs"] or result["svs"]),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("sample", "input-mode", "inputs", "reference", "reference-check", "exons", "regions",
                 "sv", "cnv", "cnr", "cov-q20", "cov-thresholds", "cov-q0", "cov-perbase",
                 "chrom-median", "versions", "commands", "out-json", "out-md", "out-html"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--small", action="append", default=[])
    p.add_argument("--snv-caller", default="haplotypecaller")
    p.add_argument("--min-exon-depth", type=float, default=20)
    p.add_argument("--min-exon-pct20", type=float, default=95)
    p.add_argument("--low-depth", type=float, default=10)
    p.add_argument("--max-mapq0", type=float, default=0.10)
    p.add_argument("--nextflow-version", default=None)
    p.add_argument("--session-id", default=None)
    p.add_argument("--lang", choices=("en", "de"), default="en", help="language of the HTML/Markdown report")
    args = p.parse_args(argv)
    result = build(args)
    with open(args.out_json, "w") as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False)
    with open(args.out_md, "w") as handle:
        handle.write(render.markdown(result, args.lang))
    with open(args.out_html, "w") as handle:
        handle.write(render.html_report(result, args.lang))
    return 0


if __name__ == "__main__":
    sys.exit(main())

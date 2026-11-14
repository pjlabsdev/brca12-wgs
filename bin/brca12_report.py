#!/usr/bin/env python3
"""Fasst alle BRCA1/BRCA2-Ergebnisse zu brca12_results.json, brca12_report.md und .html zusammen.

The JSON file (schema genome.brca12.v1) is the machine-readable interface for downstream reports.
"""
import argparse
import datetime
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from brca12lib import coverage, enigma_tables, evidence, render, structural, variants  # noqa: E402

SCHEMA = "genome.brca12.v1"

GENES = {
    "BRCA1": {"strand": "-", "transcript": "ENST00000357654 / NM_007294.4 (MANE Select)"},
    "BRCA2": {"strand": "+", "transcript": "ENST00000380152 / NM_000059.4 (MANE Select)"},
}

LIMITATIONS = [
    "Nur GRCh38/hg38. Koordinaten BRCA1 chr17:43044295-43170327 (Minus-Strang), BRCA2 chr13:32315077-32400268 (Plus-Strang).",
    "Kurz-Read-WGS erfasst Mosaike unter etwa 15 % Allelanteil, komplexe Umlagerungen und Bereiche mit Pseudogen-Homologie (BRCA1-Exons 1-2 und BRCA1P1) nur eingeschränkt.",
    "CNVs stammen aus einer Einzelprobe mit flacher CNVkit-Referenz und gelten nur als Kandidaten. Einzelexon-Deletionen oder -Duplikationen unter etwa 1 kb sind unsicher; MLPA ist die übliche Bestätigung.",
    "Manta meldet nur SVs mit Bruchpunkt im Suchfenster um BRCA1/BRCA2. Insertionen mobiler Elemente werden nicht sicher erkannt; die zusätzliche Suche nach Insertionssignaturen in abgeschnittenen Leseenden erkennt nur Alu-Sequenz, braucht mindestens 3 Lesungen je Bruchstelle und liefert Hinweise, keine Aufrufe.",
    "ACMG/AMP: nur maschinelle Hinweise auf einzelne Kriterien nach ClinGen ENIGMA v1.2. Keine Klassifikation; PVS1/PM5_PTC (Tabelle 4), PS3/BS3 (Tabelle 9) und PP4/BP5 stammen aus den mitgelieferten ENIGMA-Tabellen, PS1, PM3, BS2, BS4, PP1, PS4 und RNA-Daten der Probe fehlen. Die gnomAD-Lesetiefe wird nicht geprüft.",
    "ClinVar-Angaben hängen von der Version der verwendeten ClinVar-VCF ab.",
]


LIMITATIONS_EN = [
    "GRCh38/hg38 only. Coordinates BRCA1 chr17:43044295-43170327 (minus strand), BRCA2 chr13:32315077-32400268 (plus strand).",
    "Short-read WGS has limited sensitivity for mosaics below about 15% allele fraction, complex rearrangements and regions with pseudogene homology (BRCA1 exons 1-2 and BRCA1P1).",
    "CNVs come from a single sample with a flat CNVkit reference and are candidates only. Single-exon deletions or duplications below about 1 kb are uncertain; MLPA is the usual confirmation.",
    "Manta only reports SVs with a breakpoint in the search window around BRCA1/BRCA2. Mobile element insertions are not reliably detected; the additional search for insertion signatures in clipped read ends recognises Alu sequence only, needs at least 3 reads per junction and gives hints, not calls.",
    "ACMG/AMP: automated hints on individual criteria following ClinGen ENIGMA v1.2. No classification; PVS1/PM5_PTC (Table 4), PS3/BS3 (Table 9) and PP4/BP5 come from the bundled ENIGMA tables, PS1, PM3, BS2, BS4, PP1, PS4 and RNA data of the sample are missing. gnomAD read depth is not checked.",
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


def exon_number(label):
    match = re.match(r"Exon (\d+)", label or "")
    return int(match.group(1)) if match else None


def add_exon_event_codes(event, kind, labels):
    """ENIGMA Table 4 codes for whole-exon deletions/duplications (tandem status unknown)."""
    event["enigmaCodes"] = []
    if kind is None:
        return
    for label in labels:
        number = exon_number(label)
        row = enigma_tables.exon_event_row(event["gene"], number, kind) if number else None
        if row:
            event["enigmaCodes"].append({"exon": label, "event": row["variant"], "code": row["code"],
                                         "note": row.get("warning") or None})
    if event["enigmaCodes"]:
        text = "; ".join(f"{c['event']}: {c['code']}" for c in event["enigmaCodes"])
        event.setdefault("flags", []).append(f"ENIGMA Tabelle 4 (falls bestätigt): {text}")


ASSAY_LIMITATIONS = {
    "exome": ("Exom-/Panel-Daten (Hybrid-Capture): Introns und Promotor sind nur am Rand der Ziele abgedeckt; CNVs aus CNVkit hybrid mit flacher Referenz sind ohne --cnv_reference (Normalproben derselben Methode) unsicher; Manta läuft im Exom-Modus.",
              "Exome/panel data (hybrid capture): introns and the promoter are covered only at the edges of the targets; CNVs from CNVkit hybrid with a flat reference are unreliable without --cnv_reference (normal samples of the same assay); Manta runs in exome mode."),
    "amplicon": ("Amplikon-Daten: Duplikate werden nicht markiert; Varianten unter Primern sind ohne Primer-Abschnitt (--primers_bed) unzuverlässig; Allelausfall durch Varianten in Primerbindestellen ist möglich; große Deletionen und Strukturvarianten sind mit Amplikons kaum erkennbar; CNVs aus CNVkit amplicon brauchen eine Referenz aus Normalproben (--cnv_reference).",
                 "Amplicon data: duplicates are not marked; variants under primers are unreliable without primer clipping (--primers_bed); allele dropout from variants in primer binding sites is possible; large deletions and structural variants are hardly detectable with amplicons; CNVs from CNVkit amplicon need a reference from normal samples (--cnv_reference)."),
}
ASSAY_LIMITATIONS["panel"] = ASSAY_LIMITATIONS["exome"]


def mosaic_candidates(paths, germline, low_exons, sources, min_af, max_af, min_alt=3, min_depth=20,
                      snv_only=False):
    """Mutect2 PASS calls with an allele fraction in [min_af, max_af], at least `min_alt` supporting reads
    and `min_depth` reads in total, that the germline caller did not report. With `snv_only` (Ion Torrent)
    insertions and deletions are left out: homopolymer errors of the platform look like low-fraction indels."""
    known = {(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in germline}
    rows = []
    for gene, path in paths.items():
        for v in variants.small_variants(path, gene, low_exons, sources):
            fraction = v["alleleFraction"]
            if (v["chrom"], v["pos"], v["ref"], v["alt"]) in known or fraction is None:
                continue
            if not min_af <= fraction <= max_af:
                continue
            if (v["altDepth"] or 0) < min_alt or (v["depth"] or 0) < min_depth:
                continue
            if snv_only and len(v["ref"]) != len(v["alt"]):
                continue
            v["genotype"] = "Mosaik-Kandidat"
            v["qualityFlags"] = [f for f in v["qualityFlags"]
                                 if not f.startswith(("Allelverhältnis", "geringe Genotypqualität"))]
            v["qualityFlags"].append("niedriger Allelanteil: Bestätigung mit orthogonaler Methode nötig")
            rows.append(v)
    return rows


def mark_untargeted(exons, untargeted):
    names = set(untargeted or [])
    for exon in exons:
        if exon["exon"] is not None and f"{exon['gene']}_ex{exon['exon']:02d}" in names:
            exon["status"] = "unsicher"
            exon["reasons"].append("nicht im Zielbereich des Panels")


def build(args):
    with open(args.reference_check) as handle:
        check = json.load(handle)
    regions = read_exact_regions(args.regions)
    medians = coverage.read_chrom_median(args.chrom_median)
    exons = coverage.exon_table(args.cov_q20, args.cov_thresholds, args.cov_q0, medians,
                                args.min_exon_depth, args.min_exon_pct20, args.max_mapq0)
    read_qc = None
    if args.read_qc and os.path.exists(args.read_qc) and os.path.getsize(args.read_qc):
        with open(args.read_qc) as handle:
            read_qc = json.load(handle)
    if args.assay != "wgs":
        mark_untargeted(exons, (check.get("targets") or {}).get("exonsWithoutTarget"))
    uncertain = [e for e in exons if e["status"] != "ok"]
    small = dict(item.split("=", 1) for item in args.small)
    sources = evidence.Evidence(args.gnomad_exomes, args.gnomad_genomes, args.spliceai, args.bayesdel, args.tabix)
    genes = {}
    for gene, meta in GENES.items():
        chrom, start, end, _ = next(r for r in regions if r[3] == gene)
        rows = variants.small_variants(small[gene], gene, uncertain, sources) if gene in small else []
        genes[gene] = {"symbol": gene, "chrom": chrom, "start": start + 1, "end": end,
                       "strand": meta["strand"], "transcript": meta["transcript"],
                       "coverage": coverage.gene_summary(exons, gene), "smallVariants": rows}
    low = coverage.low_coverage_intervals(args.cov_perbase, regions, args.low_depth)
    cnvs = [c for c in structural.cnv_segments(args.cnv, args.cnr, regions, exons)]
    for cnv in cnvs:
        add_exon_event_codes(cnv, "DEL" if cnv["type"] == "Deletion" else "DUP" if cnv["type"] == "Duplikation" else None,
                             cnv.get("affectedExons", []))
    dosage_hints = structural.exon_dosage_hints(exons)
    for item in dosage_hints:
        add_exon_event_codes(item, "DEL" if item["type"] == "möglicher Verlust" else "DUP", [item["exon"]])
    svs = structural.structural_variants(args.sv, regions, exons)
    mosaic_paths = dict(item.split("=", 1) for item in args.mosaic)
    germline = [v for g in genes.values() for v in g["smallVariants"]]
    ion = args.platform == "iontorrent"
    mosaic = mosaic_candidates(mosaic_paths, germline, uncertain, sources, args.mosaic_min_af, args.mosaic_max_af,
                               args.mosaic_min_alt_reads, args.mosaic_min_depth, snv_only=ion) if mosaic_paths else []
    if ion:
        for v in germline:
            if len(v["ref"]) != len(v["alt"]):
                v["qualityFlags"].append("Ion Torrent: Insertionen/Deletionen in Homopolymeren unzuverlässig")
    versions = read_versions(args.versions)
    versions.setdefault("nextflow", args.nextflow_version)
    invocation = " ".join(["nextflow run main.nf", f"--sample {args.sample}",
                           f"--reference {args.reference}",
                           (f"--bam {args.inputs}" if args.input_mode == "bam"
                            else " ".join(f"--reads_r{i} {p}" for i, p in enumerate(args.inputs.split(","), 1) if p)),
                           f"--snv_caller {args.snv_caller}", f"--assay {args.assay}",
                           "-profile <conda|docker|standard>"])
    result = {
        "schema": SCHEMA,
        "sample": args.sample,
        "generatedAt": datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat(),
        "assembly": "GRCh38",
        "researchUseOnly": True,
        "disclaimer": {"de": render.DISCLAIMER_DE, "en": render.DISCLAIMER_EN},
        "clinicalClassification": None,
        "input": {"mode": args.input_mode, "files": [os.path.basename(p) for p in args.inputs.split(",") if p],
                  "singleEnd": args.input_mode == "fastq" and len([p for p in args.inputs.split(",") if p]) == 1},
        "reference": {"path": args.reference, "assembly": check.get("assembly"), "check": check},
        "assay": args.assay,
        "snvCaller": "GATK HaplotypeCaller" if args.snv_caller == "haplotypecaller" else "DeepVariant",
        "readQC": read_qc,
        "platform": args.platform,
        "mosaicCalling": {"enabled": bool(mosaic_paths), "caller": "GATK Mutect2 (tumour-only)" if mosaic_paths else None,
                          "minAlleleFraction": args.mosaic_min_af, "maxAlleleFraction": args.mosaic_max_af,
                          "minAltReads": args.mosaic_min_alt_reads, "minDepth": args.mosaic_min_depth,
                          "snvOnly": ion},
        "mosaicCandidates": mosaic,
        "evidenceSources": dict(sources.configured(), specification="ClinGen ENIGMA BRCA1/BRCA2 VCEP v1.2"),
        "thresholds": {"minExonDepth": args.min_exon_depth, "minExonPct20x": args.min_exon_pct20,
                       "lowDepth": args.low_depth, "maxLowMapqFraction": args.max_mapq0},
        "genes": genes,
        "exonCoverage": exons,
        "lowCoverageIntervals": low,
        "cnvs": [c for c in cnvs if c["type"] != "kein Hinweis"],
        "cnvSegmentsChecked": len(cnvs),
        "exonDosageHints": dosage_hints,
        "svs": svs,
        "insertionSignals": structural.insertion_signals(args.insertion_signals, exons),
        "toolVersions": versions,
        "invocation": invocation,
        "nextflowSession": args.session_id,
        "limitations": LIMITATIONS + ([ASSAY_LIMITATIONS[args.assay][0]] if args.assay in ASSAY_LIMITATIONS else []),
        "limitationsEN": LIMITATIONS_EN + ([ASSAY_LIMITATIONS[args.assay][1]] if args.assay in ASSAY_LIMITATIONS else []),
    }
    if ion:
        result["limitations"].append("Ion Torrent: Insertionen und Deletionen in Homopolymeren sind mit GATK HaplotypeCaller unzuverlässig (plattformtypische Lesefehler); Mosaik-Kandidaten werden nur für SNVs gemeldet.")
        result["limitationsEN"].append("Ion Torrent: insertions and deletions in homopolymers are unreliable with GATK HaplotypeCaller (platform-specific read errors); mosaic candidates are reported for SNVs only.")
    if result["input"]["singleEnd"]:
        result["limitations"].append("Einzel-Reads (single-end): Manta wurde übersprungen, strukturelle Varianten werden nicht gesucht; Duplikate lassen sich nur nach Startposition erkennen.")
        result["limitationsEN"].append("Single-end reads: Manta was skipped, structural variants are not searched; duplicates can only be recognised by start position.")
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
        "insertionSignalCount": len(result["insertionSignals"]),
        "mosaicCandidateCount": len(result.get("mosaicCandidates", [])),
        "readQCFlagCount": len((result.get("readQC") or {}).get("flags", [])),
        "uncertainExonCount": sum(1 for e in result["exonCoverage"] if e["status"] != "ok"),
        "lowCoverageIntervalCount": len(result["lowCoverageIntervals"]),
        "needsAttention": bool(clinvar_plp or lof or result["cnvs"] or result["svs"] or result["insertionSignals"]
                               or result.get("mosaicCandidates")),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("sample", "input-mode", "inputs", "reference", "reference-check", "exons", "regions",
                 "sv", "cnv", "cnr", "cov-q20", "cov-thresholds", "cov-q0", "cov-perbase",
                 "chrom-median", "versions", "commands", "out-json", "out-md", "out-html"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--small", action="append", default=[])
    p.add_argument("--insertion-signals", default=None)
    p.add_argument("--snv-caller", default="haplotypecaller")
    p.add_argument("--assay", choices=("wgs", "exome", "panel", "amplicon"), default="wgs")
    p.add_argument("--read-qc", default=None, help="read_qc.json from brca12_read_qc.py")
    p.add_argument("--mosaic", action="append", default=[], help="GENE=mosaic.vep.vcf.gz (Mutect2, VEP)")
    p.add_argument("--mosaic-min-af", type=float, default=0.02)
    p.add_argument("--mosaic-max-af", type=float, default=0.30)
    p.add_argument("--mosaic-min-alt-reads", type=int, default=3)
    p.add_argument("--mosaic-min-depth", type=int, default=20)
    p.add_argument("--platform", choices=("illumina", "iontorrent"), default="illumina")
    p.add_argument("--min-exon-depth", type=float, default=20)
    p.add_argument("--min-exon-pct20", type=float, default=95)
    p.add_argument("--low-depth", type=float, default=10)
    p.add_argument("--max-mapq0", type=float, default=0.10)
    p.add_argument("--nextflow-version", default=None)
    p.add_argument("--session-id", default=None)
    p.add_argument("--gnomad-exomes", default=None, help="gnomAD exomes VCF with non-cancer FAF (bgzip + tabix)")
    p.add_argument("--gnomad-genomes", default=None, help="gnomAD genomes VCF with non-cancer FAF (bgzip + tabix)")
    p.add_argument("--spliceai", default=None, help="SpliceAI precomputed scores VCF (bgzip + tabix)")
    p.add_argument("--bayesdel", default=None, help="TSV chrom,pos,ref,alt,BayesDel_noAF (bgzip + tabix)")
    p.add_argument("--tabix", default="tabix")
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

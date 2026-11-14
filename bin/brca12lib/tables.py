"""Table rows from brca12_results.json, shared by the XLSX export and the cohort summary."""

SMALL_HEADER = ["sample", "gene", "chrom", "pos_GRCh38", "pos_GRCh37", "ref", "alt", "genotype", "depth",
                "allele_fraction", "gq", "filter", "consequence", "hgvsc", "hgvsp", "bic_names", "exon",
                "clinvar", "clinvar_stars", "acmg_hints", "quality_flags"]


def _join(items):
    return "; ".join(str(i) for i in items if i not in (None, ""))


def small_rows(result):
    rows = []
    for gene in ("BRCA1", "BRCA2"):
        for v in result.get("genes", {}).get(gene, {}).get("smallVariants", []):
            clinvar = v.get("clinvar") or {}
            rows.append([
                result.get("sample"), gene, v.get("chrom"), v.get("pos"), v.get("grch37Pos"), v.get("ref"),
                v.get("alt"), v.get("genotype"), v.get("depth"), v.get("alleleFraction"), v.get("gq"),
                v.get("filter"), v.get("consequence"), v.get("hgvsc"), v.get("hgvsp"),
                _join(v.get("bicNames") or []), v.get("exon"), clinvar.get("significance"),
                clinvar.get("stars"), _join(h.get("code") for h in v.get("acmgHints") or []),
                _join(v.get("qualityFlags") or []),
            ])
    return rows


def sample_sheets(result):
    """Sheets for one sample: summary, small variants, CNV, SV, coverage and tool versions."""
    sample = result.get("sample")
    summary = result.get("summary") or {}
    disclaimer = (result.get("disclaimer") or {}).get("de", "")
    overview = [["Feld / field", "Wert / value"],
                ["Probe / sample", sample],
                ["Referenz / reference", result.get("assembly")],
                ["Erstellt / created", result.get("generatedAt")],
                ["SNV-Caller", result.get("snvCaller")],
                ["Nur Forschung (RUO)", disclaimer]]
    overview += [[key, value] for key, value in summary.items()]
    cnv = [["sample", "gene", "type", "chrom", "start", "end", "copy_number", "log2", "affected_exons", "flags"]]
    cnv += [[sample, c.get("gene"), c.get("type"), c.get("chrom"), c.get("start"), c.get("end"),
             c.get("copyNumber"), c.get("log2"), _join(c.get("affectedExons") or []), _join(c.get("flags") or [])]
            for c in result.get("cnvs") or []]
    cnv += [[sample, d.get("gene"), d.get("type"), None, None, None, None, d.get("ratio"), d.get("exon"),
             _join([d.get("note")] + (d.get("flags") or []))] for d in result.get("exonDosageHints") or []]
    sv = [["sample", "gene", "type", "chrom", "pos", "end", "mate", "genotype", "paired_reads", "split_reads",
           "filter", "affected_exons", "flags"]]
    sv += [[sample, s.get("gene"), s.get("type"), s.get("chrom"), s.get("pos"), s.get("end"),
            f"{s.get('mateChrom')}:{s.get('matePos')}" if s.get("mateChrom") else None, s.get("genotype"),
            s.get("pairedReads"), s.get("splitReads"), s.get("filter"), _join(s.get("affectedExons") or []),
            _join(s.get("flags") or [])] for s in result.get("svs") or []]
    exons = [["sample", "gene", "exon", "chrom", "start", "end", "mean_depth", "pct_10x", "pct_20x", "pct_30x",
              "low_mapq_fraction", "dosage_ratio", "status", "reasons"]]
    exons += [[sample, e.get("gene"), e.get("label"), e.get("chrom"), e.get("start"), e.get("end"),
               e.get("meanDepth"), e.get("pct10x"), e.get("pct20x"), e.get("pct30x"), e.get("lowMapqFraction"),
               e.get("dosageRatio"), e.get("status"), _join(e.get("reasons") or [])]
              for e in result.get("exonCoverage") or []]
    low = [["sample", "gene", "chrom", "start", "end", "length", "min_depth"]]
    low += [[sample, i.get("gene"), i.get("chrom"), i.get("start"), i.get("end"), i.get("length"), i.get("minDepth")]
            for i in result.get("lowCoverageIntervals") or []]
    versions = [["tool", "version"]] + [[k, v] for k, v in (result.get("toolVersions") or {}).items()]
    sheets = [("Übersicht", overview), ("SNVs Indels", [SMALL_HEADER] + small_rows(result)),
              ("CNV", cnv), ("SV", sv), ("Exon-Coverage", exons), ("Geringe Tiefe", low)]
    if (result.get("mosaicCalling") or {}).get("enabled"):
        mosaic = [["sample", "gene", "chrom", "pos", "ref", "alt", "allele_fraction", "alt_reads", "depth",
                   "consequence", "hgvsc", "hgvsp", "clinvar", "flags"]]
        mosaic += [[sample, v.get("gene"), v.get("chrom"), v.get("pos"), v.get("ref"), v.get("alt"),
                    v.get("alleleFraction"), v.get("altDepth"), v.get("depth"), v.get("consequence"), v.get("hgvsc"),
                    v.get("hgvsp"), (v.get("clinvar") or {}).get("significance"), _join(v.get("qualityFlags") or [])]
                   for v in result.get("mosaicCandidates") or []]
        sheets.append(("Mosaik-Kandidaten", mosaic))
    qc = result.get("readQC")
    if qc:
        rows = [["metric", "value"], ["assay", result.get("assay")], ["platform", result.get("platform")]]
        rows += [[k, v] for k, v in (qc.get("regionStats") or {}).items()]
        rows += [[f"target {k}", v] for k, v in (qc.get("targetSummary") or {}).items()]
        rows += [[f"primer clipping {k}", v] for k, v in (qc.get("primerClipping") or {}).items()]
        rows += [["warning", f] for f in qc.get("flags") or []]
        sheets.append(("Lese-QC", rows))
        if qc.get("targets"):
            targets = [["sample", "target", "chrom", "start", "end", "mean_depth", "min_depth", "pct_20x", "status"]]
            targets += [[sample, t.get("name"), t.get("chrom"), t.get("start"), t.get("end"), t.get("meanDepth"),
                         t.get("minDepth"), t.get("pct20x"), t.get("status")] for t in qc["targets"]]
            sheets.append(("Ziel-Coverage", targets))
    return sheets + [("Tool-Versionen", versions)]

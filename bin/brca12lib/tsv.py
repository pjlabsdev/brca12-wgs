"""PREFIX.tsv (one line per gene), PREFIX.findings.tsv (one line per finding with its certainty) and
PREFIX.exons.tsv (one line per exon with its certainty), in the form of hla-wgs. Values are English and language-neutral; the full detail stays in PREFIX.json."""
from . import page
from .i18n import tr

GENE_COLUMNS = ("gene", "result", "statement", "attention", "unclear", "small_variants", "copy_number", "structural",
                "insertion_signals", "mosaic", "coverage", "exons", "exons_uncertain", "mean_depth",
                "min_exon_depth", "low_depth_bp", "transcript")

EXON_COLUMNS = ("gene", "exon", "label", "chrom", "start", "end", "mean_depth", "pct10x", "pct20x", "pct30x",
                "low_mapq_fraction", "dosage_ratio", "status", "reasons")

FINDING_COLUMNS = ("gene", "kind", "level", "reason", "chrom", "start", "end", "change", "hgvsc", "hgvsp",
                "exon", "genotype", "reads", "allele_fraction", "clinvar", "clinvar_stars", "acmg", "certainty",
                "flags")

# reason keys of page.notable() and the words used in the TSV
REASONS = {"kind_plp": "clinvar_pathogenic", "kind_lof": "loss_of_function", "kind_vus": "clinvar_vus",
           "kind_mosaic": "mosaic", "kind_cnv": "cnv", "kind_exon": "exon_copy_number", "kind_sv": "sv",
           "kind_ins": "insertion_signal"}

# read counts an insertion signal may carry, in the order they are written
SIGNAL_READS = ("realignedReads", "markerReads", "outwardPairs", "splitReads", "discordantPairs",
                "clippedReadsLeftJunction", "clippedReadsRightJunction")


def _cell(value):
    if value is None or value == "" or value == []:
        return "-"
    if isinstance(value, float):
        return f"{value:g}"
    return str(value).replace("\t", " ").replace("\n", " ")


CODE_WORDS = ((" (Näherung)", " (approximation)"), (" (Kandidat)", " (candidate)"), (" offen", " open"))


def _code(code):
    """ACMG/AMP code in English: the German suffixes of the JSON codes translated."""
    for de, en in CODE_WORDS:
        code = code.replace(de, en)
    return code


def _flags(item):
    return "; ".join(str(tr(f, "en")) for f in (item.get("qualityFlags") or []) + (item.get("flags") or []))


def _gene_coverage(result):
    """{gene: (exons, uncertain exons, mean depth, lowest exon depth)} over the exons with a depth."""
    out = {}
    for gene in page.GENES:
        exons = [e for e in result["exonCoverage"] if e["gene"] == gene and e.get("meanDepth") is not None]
        depths = [e["meanDepth"] for e in exons]
        out[gene] = (len(exons), sum(1 for e in exons if e["status"] != "ok"),
                     round(sum(depths) / len(depths), 1) if depths else None, min(depths) if depths else None)
    return out


STATEMENT = {
    "attention": ("Befund, der Bestätigung braucht ({n})", "finding that needs confirmation ({n})"),
    "unclear1": ("kein Befund, 1 Eintrag mit unklarer Bedeutung", "no finding, 1 entry of uncertain significance"),
    "unclear": ("kein Befund, {n} Einträge mit unklarer Bedeutung", "no finding, {n} entries of uncertain significance"),
    "none": ("kein Befund", "no finding"),
    "depth": ("Exons {depth}x", "exons {depth}x"),
    "no_depth": ("keine Exon-Tiefe", "no exon depth"),
    "uncertain": ("unsicher: {exons}", "uncertain: {exons}"),
    "all": ("alle Exons abgedeckt", "all exons covered"),
}


def statement(row, uncertain_labels, lang="en"):
    """The line for the gene as one sentence: result, then how completely the gene was checked."""
    def t(key, **kw):
        return STATEMENT[key][0 if lang == "de" else 1].format(**kw)
    if row["result"] == "attention":
        what = t("attention", n=row["attention"])
    elif row["result"] == "unclear":
        what = t("unclear1") if row["unclear"] == 1 else t("unclear", n=row["unclear"])
    else:
        what = t("none")
    if row["mean_depth"] is None:
        cover = t("no_depth")
    else:
        rest = t("uncertain", exons=", ".join(str(tr(x, lang)) for x in uncertain_labels)) if uncertain_labels \
            else t("all")
        cover = t("depth", depth=f"{row['mean_depth']:g}") + "; " + rest
    return f"{row['gene']}: {what}; {cover}"


def uncertain_labels(result, gene):
    return [e["label"] for e in result["exonCoverage"] if e["gene"] == gene and e["status"] != "ok"]


def certainty(kind, item):
    """ok: the call carries no quality flag; limited: it does, or it is a mosaic candidate or an SV without
    FILTER PASS; hint: a read signal or a depth ratio, not a call."""
    if kind in ("insertion_signal", "exon_dosage", "exon_copy_number_hint"):
        return "hint"
    if kind == "mosaic":
        return "limited"
    if kind == "sv":
        return "ok" if item.get("filter") == "PASS" else "limited"
    if kind == "small":
        return "limited" if item.get("qualityFlags") or item.get("filter") not in (None, "PASS") else "ok"
    return "limited" if item.get("flags") else "ok"


def gene_rows(result):
    levels = page.gene_levels(result)
    entries = page.notable(result)
    coverage = _gene_coverage(result)
    minimum = result["thresholds"]["minExonDepth"]
    ecn = result.get("exonCopyNumber") or {}
    rows = []
    for gene in page.GENES:
        mine = [e for e in entries if e[1] == gene]
        exons, uncertain, mean, lowest = coverage[gene]
        limited = uncertain or mean is None or mean < minimum
        low_bp = sum(i.get("length") or (i["end"] - i["start"] + 1)
                     for i in result["lowCoverageIntervals"] if i["gene"] == gene)
        rows.append({
            "gene": gene,
            "result": levels[gene],
            "attention": sum(1 for e in mine if e[0] == "attention"),
            "unclear": sum(1 for e in mine if e[0] == "unclear"),
            "small_variants": len(result["genes"][gene]["smallVariants"]),
            "copy_number": sum(1 for c in result["cnvs"] + ecn.get("calls", []) if c["gene"] == gene),
            "structural": sum(1 for s in result["svs"] if s["gene"] == gene),
            "insertion_signals": sum(1 for i in result.get("insertionSignals", []) if i["gene"] == gene),
            "mosaic": sum(1 for v in result.get("mosaicCandidates", []) if v["gene"] == gene),
            "coverage": "limited" if limited else "ok",
            "exons": exons,
            "exons_uncertain": uncertain,
            "mean_depth": mean,
            "min_exon_depth": lowest,
            "low_depth_bp": low_bp,
            "transcript": result["genes"][gene].get("transcript"),
        })
        row = rows[-1]
        row["statement"] = statement(row, uncertain_labels(result, gene))
    return rows


def change(pos, ref, alt):
    """(start, end, change) as in hla-wgs: the shared anchor bases of the VCF record removed, then
    `C>T` for one base, `delT` for a deletion (start and end on the deleted bases), `insT` for an
    insertion (start and end on the base before it), `CT>GA` for anything else."""
    p = 0
    while p < min(len(ref), len(alt)) - 1 and ref[p] == alt[p]:
        p += 1
    s = 0
    while s < min(len(ref), len(alt)) - p - 1 and ref[-1 - s] == alt[-1 - s]:
        s += 1
    r, a = ref[p:len(ref) - s], alt[p:len(alt) - s]
    start = pos + p
    if len(r) == 1 and len(a) > 1 and a[0] == r:          # insertion after the anchor base
        return start, start, "ins" + a[1:]
    if len(a) == 1 and len(r) > 1 and r[0] == a:          # deletion after the anchor base
        return start + 1, start + len(r) - 1, "del" + r[1:]
    return start, start + len(r) - 1, f"{r}>{a}"


def _small(gene, kind, v, level, reason, certainty):
    start, end, what = change(v["pos"], v["ref"], v["alt"])
    stars = (v.get("clinvar") or {}).get("stars")
    return {"gene": gene, "kind": kind, "level": level, "reason": reason, "chrom": v["chrom"],
            "start": start, "end": end, "change": what,
            "hgvsc": v.get("hgvsc"), "hgvsp": v.get("hgvsp"), "exon": v.get("exon") or v.get("intron"),
            "genotype": tr(v.get("genotype"), "en"),
            "reads": f"{v['altDepth']}/{v['depth']}" if v.get("altDepth") is not None and v.get("depth") is not None
            else v.get("depth"),
            "allele_fraction": v.get("alleleFraction"),
            "clinvar": (v.get("clinvar") or {}).get("significance"), "clinvar_stars": stars,
            "acmg": "; ".join(_code(h["code"]) for h in v.get("acmgHints", [])),
            "certainty": certainty, "flags": _flags(v)}


def _region(gene, kind, item, level, reason, certainty, change, reads, genotype=None):
    return {"gene": gene, "kind": kind, "level": level, "reason": reason, "chrom": item["chrom"],
            "start": item.get("start", item.get("pos")), "end": item.get("end"), "change": change,
            "exon": ", ".join(str(tr(e, "en")) for e in item.get("affectedExons", [])),
            "genotype": genotype, "reads": reads, "certainty": certainty, "flags": _flags(item)}


def finding_rows(result):
    """Every reported finding, in the order of the report: small variants, mosaic candidates, copy number,
    structural variants, insertion signals, exon dosage hints.

    level: attention, unclear (as on the result page), common (a known common allele), hint (a lead that
    does not count as a finding) or none (a listed variant without a finding, for example benign).
    certainty: ok when the call carries no quality flag, limited when it does or when it is a mosaic
    candidate or an SV without FILTER PASS; hint for read signals and depth ratios."""
    tier = {id(e[3]): (e[0], REASONS[e[2]]) for e in page.notable(result)}
    rows = []
    for gene in page.GENES:
        for v in result["genes"][gene]["smallVariants"]:
            level, reason = tier.get(id(v), ("none", "ba1" if page._is_ba1(v) else "listed"))
            rows.append(_small(gene, "small", v, level, reason, certainty("small", v)))
    for v in result.get("mosaicCandidates", []):
        level, reason = tier.get(id(v), ("unclear", "mosaic"))
        rows.append(_small(v["gene"], "mosaic", v, level, reason, certainty("mosaic", v)))
    for c in result["cnvs"]:
        level, reason = tier[id(c)]
        rows.append(_region(c["gene"], "cnv", c, level, reason, certainty("cnv", c),
                            f"{tr(c['type'], 'en')} CN {c['copyNumber']}", f"log2 {c.get('log2')}"))
    ecn = result.get("exonCopyNumber") or {}
    for c in ecn.get("calls", []):
        level, reason = tier[id(c)]
        rows.append(_region(c["gene"], "exon_copy_number", c, level, reason, certainty("exon_copy_number", c),
                            f"{tr(c['type'], 'en')} CN {c['copyNumber']}", f"quality {c.get('quality')}"))
    for c in ecn.get("hints", []):
        rows.append(_region(c["gene"], "exon_copy_number", c, "hint", "exon_copy_number_hint", "hint",
                            f"{tr(c['type'], 'en')} CN {c['copyNumber']}", f"quality {c.get('quality')}"))
    for s in result["svs"]:
        level, reason = tier[id(s)]
        rows.append(_region(s["gene"], "sv", s, level, reason, certainty("sv", s),
                            f"{s['type']} {s.get('caller', '')}".strip(),
                            f"PR {s.get('pairedReads')} SR {s.get('splitReads')}", s.get("genotype")))
    for i in result.get("insertionSignals", []):
        level, reason = tier.get(id(i), ("common", "common_allele"))
        reads = "; ".join(f"{k} {i[k]}" for k in SIGNAL_READS if i.get(k) not in (None, "", []))
        row = _region(i["gene"], "insertion_signal", i, level, reason, "hint", tr(i["type"], "en"), reads,
                      tr(i.get("zygosity"), "en"))
        if i.get("label"):
            row["change"] = f"{row['change']}: {tr(i['label'], 'en')}"
        rows.append(row)
    for d in result.get("exonDosageHints", []):
        item = next((e for e in result["exonCoverage"] if e["gene"] == d["gene"] and e["label"] == d["exon"]), {})
        row = _region(d["gene"], "exon_dosage", dict(item, flags=d.get("flags")), "hint", "dosage_ratio", "hint",
                      f"{tr(d['type'], 'en')} ratio {d['ratio']}", None)
        row["exon"] = tr(d["exon"], "en")
        rows.append(row)
    return rows


def _table(columns, rows):
    return "".join("\t".join(_cell(r.get(c)) for c in columns) + "\n" for r in [dict(zip(columns, columns))] + rows)


def genes_tsv(result):
    return _table(GENE_COLUMNS, gene_rows(result))


def findings_tsv(result):
    return _table(FINDING_COLUMNS, finding_rows(result))


def exon_rows(result):
    return [{"gene": e["gene"], "exon": e.get("exon"), "label": tr(e.get("label"), "en"), "chrom": e["chrom"],
             "start": e["start"], "end": e["end"], "mean_depth": e.get("meanDepth"), "pct10x": e.get("pct10x"),
             "pct20x": e.get("pct20x"), "pct30x": e.get("pct30x"), "low_mapq_fraction": e.get("lowMapqFraction"),
             "dosage_ratio": e.get("dosageRatio"), "status": "ok" if e["status"] == "ok" else "uncertain",
             "reasons": "; ".join(str(tr(x, "en")) for x in e.get("reasons", []))}
            for e in result["exonCoverage"]]


def exons_tsv(result):
    return _table(EXON_COLUMNS, exon_rows(result))

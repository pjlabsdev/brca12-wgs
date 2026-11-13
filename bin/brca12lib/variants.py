"""SNV/Indel-Aufbereitung: VEP-Konsequenz, ClinVar, ACMG/AMP-Hinweise und Qualitätsmarken.

Die ACMG/AMP-Hinweise (brca12lib.enigma) sind maschinelle Vorprüfungen einzelner Kriterien. Sie ergeben
bewusst KEINE Klassifikation (pathogen, VUS, benigne). Massgeblich sind die Spezifikationen
des ClinGen ENIGMA BRCA1/BRCA2 Variant Curation Expert Panel und eine Bewertung in einem
akkreditierten Labor.
"""
import re

from .enigma import acmg_hints
from .vcf import VCF, genotype_label, to_float, to_int

CLINVAR_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_single_submitter": 1,
}


def pick_consequence(record, allele, gene):
    entries = [c for c in record["csq"] if c.get("SYMBOL") == gene] or record["csq"]
    if not entries:
        return {}
    for predicate in (lambda c: c.get("PICK") == "1", lambda c: c.get("MANE_SELECT"),
                      lambda c: c.get("CANONICAL") == "YES"):
        for entry in entries:
            if predicate(entry) and entry.get("Allele") in (allele, "-", None, ""):
                return entry
        for entry in entries:
            if predicate(entry):
                return entry
    return entries[0]


def max_af(csq):
    values = [to_float(csq.get(k)) for k in ("MAX_AF", "gnomADe_AF", "gnomADg_AF")]
    values = [v for v in values if v is not None]
    return max(values) if values else None


def clinvar_entry(csq):
    significance = (csq.get("ClinVar_CLNSIG") or "").replace("_", " ").replace("&", ", ")
    review = csq.get("ClinVar_CLNREVSTAT") or ""
    review_key = review.replace("&", ",")
    stars = CLINVAR_STARS.get(review_key, 0 if review else None)
    return {
        "significance": significance or None,
        "reviewStatus": review_key.replace("_", " ") or None,
        "stars": stars,
        "condition": (csq.get("ClinVar_CLNDN") or "").replace("_", " ").replace("&", ", ") or None,
    }


def quality_flags(record, depth, gq, fraction, genotype, low_exons, mq):
    flags = []
    if record["filter"] not in ("PASS", "."):
        flags.append(f"FILTER={record['filter']}")
    if depth is not None and depth < 15:
        flags.append("geringe Tiefe (DP < 15)")
    if gq is not None and gq < 20:
        flags.append("geringe Genotypqualität (GQ < 20)")
    if genotype == "heterozygot" and fraction is not None and not 0.25 <= fraction <= 0.75:
        flags.append("Allelverhältnis auffällig (mögliches Mosaik oder Artefakt)")
    if mq is not None and mq < 40:
        flags.append("geringe Mapping-Qualität (MQ < 40)")
    for exon in low_exons:
        if exon["chrom"] == record["chrom"] and exon["start"] <= record["pos"] <= exon["end"]:
            flags.append(f"liegt in unsicher abgedecktem Bereich {exon['label']}")
    return flags


def small_variants(path, gene, low_exons=(), evidence=None):
    vcf = VCF(path)
    rows = []
    for record in vcf.records:
        sample = record["sample"]
        gt = sample.get("GT", "")
        genotype = genotype_label(gt)
        if genotype == "homozygot Referenz":
            continue
        ad = [to_int(v) for v in sample.get("AD", "").split(",")] if sample.get("AD") else []
        depth = to_int(sample.get("DP")) or (sum(v for v in ad if v) if ad else None)
        alt_depth = ad[1] if len(ad) > 1 else None
        fraction = (alt_depth / depth) if alt_depth is not None and depth else to_float(sample.get("VAF"))
        gq = to_int(sample.get("GQ"))
        for alt in record["alts"]:
            if alt in ("*", "<NON_REF>"):
                continue
            csq = pick_consequence(record, alt, gene)
            af = max_af(csq)
            found = evidence.for_variant(record["chrom"], record["pos"], record["ref"], alt, gene) if evidence else {}
            gnomad = found.get("gnomad")
            rows.append({
                "gene": gene, "chrom": record["chrom"], "pos": record["pos"],
                "ref": record["ref"], "alt": alt,
                "genotype": genotype, "gt": gt, "depth": depth, "altDepth": alt_depth,
                "alleleFraction": round(fraction, 3) if fraction is not None else None,
                "gq": gq, "filter": record["filter"],
                "variantClass": csq.get("VARIANT_CLASS") or None,
                "consequence": (csq.get("Consequence") or "").replace("&", ", ") or None,
                "impact": csq.get("IMPACT") or None,
                "transcript": csq.get("Feature") or None,
                "maneSelect": csq.get("MANE_SELECT") or None,
                "hgvsc": csq.get("HGVSc") or None, "hgvsp": csq.get("HGVSp") or None,
                "exon": csq.get("EXON") or None, "intron": csq.get("INTRON") or None,
                "maxAF": af,
                "gnomadFAF": gnomad["best"] if gnomad else None,
                "spliceAI": found.get("spliceai"), "bayesDel": found.get("bayesdel"),
                "clinvar": clinvar_entry(csq),
                "acmgHints": acmg_hints(gene, csq, af, record["ref"], alt, found),
                "qualityFlags": quality_flags(record, depth, gq, fraction, genotype, low_exons,
                                              to_float(record["info"].get("MQ"))),
            })
    return rows

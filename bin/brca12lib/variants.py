"""SNV/Indel-Aufbereitung: VEP-Konsequenz, ClinVar, ACMG/AMP-Hinweise und Qualitätsmarken.

Die ACMG/AMP-Hinweise sind maschinelle Vorprüfungen einzelner Kriterien. Sie ergeben
bewusst KEINE Klassifikation (pathogen, VUS, benigne). Massgeblich sind die Spezifikationen
des ClinGen ENIGMA BRCA1/BRCA2 Variant Curation Expert Panel und eine Bewertung in einem
akkreditierten Labor.
"""
import re

from .vcf import VCF, genotype_label, to_float, to_int

LOF_CONSEQUENCES = ("stop_gained", "frameshift_variant", "splice_acceptor_variant",
                    "splice_donor_variant", "start_lost", "transcript_ablation")
SPLICE_NEAR = ("splice_region_variant", "splice_donor_5th_base_variant",
               "splice_donor_region_variant", "splice_polypyrimidine_tract_variant")

CLINVAR_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_single_submitter": 1,
}


def hint(code, de, en, direction):
    return {"code": code, "textDE": de, "textEN": en, "direction": direction}


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


def is_last_exon(exon_field):
    match = re.match(r"(\d+)(?:-\d+)?/(\d+)", exon_field or "")
    return bool(match) and match.group(1) == match.group(2)


def acmg_hints(gene, csq, af):
    consequences = (csq.get("Consequence") or "").split("&")
    hints = []
    if any(c in LOF_CONSEQUENCES for c in consequences):
        if is_last_exon(csq.get("EXON")):
            hints.append(hint("PVS1 (eingeschränkt)",
                              "Loss-of-Function-Konsequenz im letzten Exon: NMD-Escape möglich, PVS1-Stärke nach dem ENIGMA-Entscheidungsbaum prüfen.",
                              "Loss-of-function consequence in the last exon: NMD escape possible, check PVS1 strength with the ENIGMA decision tree.",
                              "pathogenic"))
        else:
            hints.append(hint("PVS1 (Kandidat)",
                              "Vorzeitiger Stopp, Leserasterverschiebung oder kanonische Spleißstelle im MANE-Transkript. PVS1 nach ENIGMA-Spezifikation prüfen.",
                              "Premature stop, frameshift or canonical splice site in the MANE transcript. Check PVS1 against the ENIGMA specification.",
                              "pathogenic"))
        if gene == "BRCA2" and "Lys3326Ter" in (csq.get("HGVSp") or ""):
            hints.append(hint("Bekannte Ausnahme",
                              "BRCA2 p.Lys3326Ter (c.9976A>T) gilt trotz Stoppcodon als benigne für hohes HBOC-Risiko.",
                              "BRCA2 p.Lys3326Ter (c.9976A>T) is considered benign for high HBOC risk despite the stop codon.",
                              "benign"))
    if af is None or af == 0:
        hints.append(hint("PM2_Supporting (Kandidat)",
                          "In gnomAD nicht beobachtet.", "Not observed in gnomAD.", "pathogenic"))
    elif af > 0.001:
        hints.append(hint("BA1 (Kandidat)",
                          f"Höchste gnomAD-Populationsfrequenz {af:.4g} liegt über 0,1 % (ENIGMA-Schwelle; Filterfrequenz prüfen).",
                          f"Highest gnomAD population frequency {af:.4g} is above 0.1% (ENIGMA threshold; verify filtering allele frequency).",
                          "benign"))
    elif af > 0.0001:
        hints.append(hint("BS1 (Kandidat)",
                          f"Höchste gnomAD-Populationsfrequenz {af:.4g} liegt über 0,01 %.",
                          f"Highest gnomAD population frequency {af:.4g} is above 0.01%.",
                          "benign"))
    if "missense_variant" in consequences:
        hints.append(hint("PP3/BP4 offen",
                          "Missense: nur mit kalibrierten Vorhersagen (BayesDel, SpliceAI) bewerten. SIFT/PolyPhen sind hier nur informativ.",
                          "Missense: assess only with calibrated predictors (BayesDel, SpliceAI). SIFT/PolyPhen are informative only.",
                          "neutral"))
    if any(c in SPLICE_NEAR for c in consequences):
        hints.append(hint("Spleißeffekt prüfen",
                          "Nahe einer Spleißstelle: Spleißvorhersage (SpliceAI) und gegebenenfalls RNA-Analyse nötig.",
                          "Near a splice site: splice prediction (SpliceAI) and possibly RNA analysis required.",
                          "neutral"))
    if consequences and set(consequences) <= {"synonymous_variant", "intron_variant", "3_prime_UTR_variant",
                                              "5_prime_UTR_variant", "upstream_gene_variant",
                                              "downstream_gene_variant"}:
        hints.append(hint("BP4/BP7 (Kandidat)",
                          "Keine Proteinänderung vorhergesagt; BP7/BP4 nur nach Ausschluss eines Spleißeffekts.",
                          "No protein change predicted; BP7/BP4 only after excluding a splicing effect.",
                          "benign"))
    return hints


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


def small_variants(path, gene, low_exons=()):
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
                "clinvar": clinvar_entry(csq),
                "acmgHints": acmg_hints(gene, csq, af),
                "qualityFlags": quality_flags(record, depth, gq, fraction, genotype, low_exons,
                                              to_float(record["info"].get("MQ"))),
            })
    return rows

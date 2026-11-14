"""ACMG/AMP criterion hints following the ClinGen ENIGMA BRCA1/BRCA2 VCEP specifications v1.2.

Source: ClinGen Criteria Specification Registry, GN092 (BRCA1) and GN097 (BRCA2), version 1.2.
The hints are automated pre-checks of single criteria. They never combine into a classification.

Implemented from the specification:
  BA1, BS1, BS1_Supporting  filtering allele frequency (FAF) in gnomAD non-cancer, non-founder
                            populations, exomes and genomes considered separately
  PM2_Supporting            absent from gnomAD non-cancer; not for insertions, deletions, delins
  PP3, BP4                  BayesDel no-AF inside the functional domains, SpliceAI for splicing
  BP1_Strong                silent, missense or in-frame outside the domains, SpliceAI <= 0.1
  BP7                       silent inside a domain, or intronic at or beyond +7/-21, if BP4 is met

From the bundled ENIGMA tables (enigma_tables.py, assets/enigma):
  PVS1, PM5_PTC             decision trees per exon, splice site and start codon (Table 4)
  PS3, BS3                  calibrated functional assays (Table 9)
  PP4, BP5                  published combined likelihood ratios

Not implemented: PS1, PM3, BS2, BS4 and RNA evidence of the sample itself. The gnomAD read-depth
conditions (>= 20 for BA1/BS1, >= 25 for PM2_Supporting) are not checked and are named in the text.
"""
import re

from . import enigma_tables

LOF_CONSEQUENCES = ("stop_gained", "frameshift_variant", "splice_acceptor_variant",
                    "splice_donor_variant", "start_lost", "transcript_ablation")
SPLICE_NEAR = ("splice_region_variant", "splice_donor_5th_base_variant",
               "splice_donor_region_variant", "splice_polypyrimidine_tract_variant")
INFRAME = ("inframe_insertion", "inframe_deletion", "protein_altering_variant")

# (Potentially) clinically important functional domains, amino-acid positions (specification v1.2).
DOMAINS = {
    "BRCA1": (("RING", 2, 101), ("coiled-coil", 1391, 1424), ("BRCT", 1650, 1857)),
    "BRCA2": (("PALB2 binding", 10, 40), ("DNA binding", 2481, 3186)),
}
# BayesDel no-AF: PP3 at or above the first value, BP4 at or below the second.
BAYESDEL = {"BRCA1": (0.28, 0.15), "BRCA2": (0.30, 0.18)}
SPLICEAI_PP3 = 0.2
SPLICEAI_BENIGN = 0.1
FAF_BA1 = 0.001
FAF_BS1 = 0.0001
FAF_BS1_SUPPORTING = 0.00002


def hint(code, de, en, direction):
    return {"code": code, "textDE": de, "textEN": en, "direction": direction}


def is_last_exon(exon_field):
    match = re.match(r"(\d+)(?:-\d+)?/(\d+)", exon_field or "")
    return bool(match) and match.group(1) == match.group(2)


def protein_position(csq):
    match = re.match(r"(\d+)", csq.get("Protein_position") or "")
    if not match:
        match = re.search(r"p\.\(?[A-Za-z*]+?(\d+)", csq.get("HGVSp") or "")
    return int(match.group(1)) if match else None


def domain_of(gene, position):
    if position is None:
        return None
    for name, start, end in DOMAINS.get(gene, ()):
        if start <= position <= end:
            return name
    return None


def intron_offset(hgvsc):
    """Signed distance from the nearest exon for an intronic HGVSc, e.g. c.212+7A>G gives 7, c.213-21 gives -21."""
    match = re.search(r"c\.[*-]?\d+([+-])(\d+)", hgvsc or "")
    if not match:
        return None
    return int(match.group(2)) * (1 if match.group(1) == "+" else -1)


def _num(value):
    return f"{value:.3g}"


def _de(value):
    return _num(value).replace(".", ",")


def frequency_hints(is_substitution, evidence, proxy_af):
    """BA1/BS1/BS1_Supporting/PM2_Supporting. `evidence` is the gnomAD lookup result or None."""
    hints = []
    depth_de = "Die gnomAD-Lesetiefe an der Position ist nicht geprüft"
    depth_en = "gnomAD read depth at the position was not checked"
    if evidence and evidence["datasets"]:
        best = evidence["best"]
        names = ", ".join(evidence["datasets"])
        if best and best["faf"] is not None and best["faf"] > FAF_BS1_SUPPORTING:
            faf, where = best["faf"], f"{best['dataset']}, {best['population']}"
            subset_de = "" if best["nonCancer"] else " Achtung: keine non-cancer-Werte in der Datei, verwendet wurde die Gesamtkohorte."
            subset_en = "" if best["nonCancer"] else " Note: no non-cancer values in the file, the full cohort was used."
            if faf > FAF_BA1:
                code, de, en = "BA1", "über 0,1 %", "above 0.1%"
            elif faf > FAF_BS1:
                code, de, en = "BS1", "über 0,01 %", "above 0.01%"
            else:
                code, de, en = "BS1_Supporting", "über 0,002 % und höchstens 0,01 %", "above 0.002% and at most 0.01%"
            hints.append(hint(code,
                              f"Filterfrequenz (FAF) {_de(faf)} in {where} liegt {de}. {depth_de} (Regel: nicht anwenden unter 20x). Nicht für bekannte pathogene Gründervarianten.{subset_de}",
                              f"Filtering allele frequency (FAF) {_num(faf)} in {where} is {en}. {depth_en} (rule: do not apply below 20x). Not for well-established pathogenic founder variants.{subset_en}",
                              "benign"))
        elif not evidence["present"]:
            if len(evidence["datasets"]) < 2:
                hints.append(hint("PM2_Supporting offen",
                                  f"Nicht beobachtet in {names}. PM2_Supporting verlangt gnomAD-Exome und -Genome; die zweite Datei fehlt.",
                                  f"Absent from {names}. PM2_Supporting requires gnomAD exomes and genomes; the second file is missing.",
                                  "neutral"))
            elif is_substitution:
                hints.append(hint("PM2_Supporting",
                                  f"In gnomAD nicht beobachtet ({names}). {depth_de} (Regel: mittlere Tiefe mindestens 25x).",
                                  f"Absent from gnomAD ({names}). {depth_en} (rule: average depth of at least 25x).",
                                  "pathogenic"))
            else:
                hints.append(hint("PM2_Supporting nicht anwendbar",
                                  "In gnomAD nicht beobachtet, aber PM2_Supporting gilt nach ENIGMA nicht für Insertionen, Deletionen und Delins.",
                                  "Absent from gnomAD, but ENIGMA does not apply PM2_Supporting to insertions, deletions or delins.",
                                  "neutral"))
        return hints
    # No gnomAD non-cancer file: only an approximation from the VEP population maximum.
    note_de = "Näherung aus der höchsten VEP-Populationsfrequenz, keine Filterfrequenz aus gnomAD non-cancer. Mit --gnomad_exomes_vcf/--gnomad_genomes_vcf wird nach ENIGMA geprüft."
    note_en = "Approximation from the highest VEP population frequency, not a filtering allele frequency from gnomAD non-cancer. Use --gnomad_exomes_vcf/--gnomad_genomes_vcf for the ENIGMA rule."
    if proxy_af is None or proxy_af == 0:
        if is_substitution:
            hints.append(hint("PM2_Supporting (Näherung)", f"In den VEP-Frequenzen nicht beobachtet. {note_de}",
                              f"Not observed in the VEP frequencies. {note_en}", "pathogenic"))
    elif proxy_af > FAF_BA1:
        hints.append(hint("BA1 (Näherung)", f"Höchste Populationsfrequenz {_de(proxy_af)} liegt über 0,1 %. {note_de}",
                          f"Highest population frequency {_num(proxy_af)} is above 0.1%. {note_en}", "benign"))
    elif proxy_af > FAF_BS1:
        hints.append(hint("BS1 (Näherung)", f"Höchste Populationsfrequenz {_de(proxy_af)} liegt über 0,01 %. {note_de}",
                          f"Highest population frequency {_num(proxy_af)} is above 0.01%. {note_en}", "benign"))
    elif proxy_af > FAF_BS1_SUPPORTING:
        hints.append(hint("BS1_Supporting (Näherung)",
                          f"Höchste Populationsfrequenz {_de(proxy_af)} liegt über 0,002 %. {note_de}",
                          f"Highest population frequency {_num(proxy_af)} is above 0.002%. {note_en}", "benign"))
    return hints


def predictor_hints(gene, csq, spliceai, bayesdel):
    """PP3, BP4, BP1_Strong and BP7 from SpliceAI and BayesDel no-AF."""
    consequences = (csq.get("Consequence") or "").split("&")
    missense = "missense_variant" in consequences
    inframe = any(c in INFRAME for c in consequences)
    silent = "synonymous_variant" in consequences
    canonical_splice = any(c in ("splice_acceptor_variant", "splice_donor_variant") for c in consequences)
    intronic = "intron_variant" in consequences and not canonical_splice
    if not (missense or inframe or silent or intronic) or any(c in LOF_CONSEQUENCES for c in consequences):
        return []
    hints = []
    position = protein_position(csq)
    domain = domain_of(gene, position)
    coding = missense or inframe or silent
    if spliceai is None:
        hints.append(hint("PP3/BP4 offen",
                          "Kein SpliceAI-Wert vorhanden (--spliceai_vcf). Ohne Spleißvorhersage sind PP3, BP4, BP1 und BP7 nicht prüfbar.",
                          "No SpliceAI score available (--spliceai_vcf). PP3, BP4, BP1 and BP7 cannot be checked without a splice prediction.",
                          "neutral"))
        return hints
    if spliceai >= SPLICEAI_PP3:
        hints.append(hint("PP3 (Spleißen)",
                          f"SpliceAI {_de(spliceai)} erreicht die Schwelle 0,2: Spleißeffekt vorhergesagt. RNA-Analyse empfohlen.",
                          f"SpliceAI {_num(spliceai)} reaches the 0.2 threshold: splicing impact predicted. RNA analysis recommended.",
                          "pathogenic"))
    no_splice = spliceai <= SPLICEAI_BENIGN
    if (missense or inframe) and domain:
        high, low = BAYESDEL[gene]
        if bayesdel is None:
            hints.append(hint("PP3/BP4 offen",
                              f"In der Domäne {domain}, aber kein BayesDel-Wert vorhanden (--bayesdel_tsv).",
                              f"Inside the {domain} domain, but no BayesDel score available (--bayesdel_tsv).", "neutral"))
        elif bayesdel >= high:
            hints.append(hint("PP3",
                              f"BayesDel no-AF {_de(bayesdel)} erreicht {_de(high)} in der Domäne {domain}.",
                              f"BayesDel no-AF {_num(bayesdel)} reaches {_num(high)} inside the {domain} domain.", "pathogenic"))
        elif bayesdel <= low and no_splice:
            hints.append(hint("BP4",
                              f"BayesDel no-AF {_de(bayesdel)} höchstens {_de(low)} und SpliceAI {_de(spliceai)} höchstens 0,1 in der Domäne {domain}.",
                              f"BayesDel no-AF {_num(bayesdel)} at most {_num(low)} and SpliceAI {_num(spliceai)} at most 0.1 inside the {domain} domain.",
                              "benign"))
    elif coding and not domain and no_splice:
        if position is None:
            hints.append(hint("BP1/BP4 offen",
                              "Proteinposition unbekannt, Lage zu den funktionellen Domänen nicht prüfbar.",
                              "Protein position unknown, location relative to the functional domains cannot be checked.", "neutral"))
        else:
            hints.append(hint("BP1_Strong",
                              f"Außerhalb der klinisch wichtigen Domänen (Position {position}) und SpliceAI {_de(spliceai)} höchstens 0,1.",
                              f"Outside the clinically important domains (position {position}) and SpliceAI {_num(spliceai)} at most 0.1.",
                              "benign"))
    elif silent and domain and no_splice:
        hints.append(hint("BP4", f"Stille Variante in der Domäne {domain}, SpliceAI {_de(spliceai)} höchstens 0,1.",
                          f"Silent variant inside the {domain} domain, SpliceAI {_num(spliceai)} at most 0.1.", "benign"))
        hints.append(hint("BP7", "Stille Variante in einer Domäne und BP4 erfüllt.",
                          "Silent variant inside a domain and BP4 met.", "benign"))
    if intronic and not coding and no_splice:
        hints.append(hint("BP4", f"Intronisch außerhalb der Positionen ±1,2, SpliceAI {_de(spliceai)} höchstens 0,1.",
                          f"Intronic outside the ±1,2 positions, SpliceAI {_num(spliceai)} at most 0.1.", "benign"))
        offset = intron_offset(csq.get("HGVSc"))
        if offset is not None and (offset >= 7 or offset <= -21):
            hints.append(hint("BP7", f"Intronisch bei {offset:+d}, also bei oder jenseits von +7/-21, und BP4 erfüllt.",
                              f"Intronic at {offset:+d}, at or beyond +7/-21, and BP4 met.", "benign"))
    return hints


def decision_tree_hints(gene, csq, consequences):
    """PVS1 and PM5_PTC from ENIGMA Table 4; a candidate hint if no table row applies."""
    hgvsc = csq.get("HGVSc")
    if any(c in ("stop_gained", "frameshift_variant") for c in consequences):
        row, what = enigma_tables.table4_row(gene, hgvsc, "ptc"), ("vorzeitiges Stoppcodon", "premature termination codon")
    elif any(c in ("splice_acceptor_variant", "splice_donor_variant") for c in consequences):
        row, what = enigma_tables.table4_row(gene, hgvsc, "splice"), ("Spleißstelle ±1,2", "splice site ±1,2")
    elif "start_lost" in consequences:
        row, what = enigma_tables.table4_row(gene, hgvsc, "start"), ("Startcodon", "initiation codon")
    else:
        row, what = None, None
    if row:
        return enigma_tables.table4_hints(row, *what)
    if is_last_exon(csq.get("EXON")):
        return [hint("PVS1 (eingeschränkt)",
                     "Loss-of-Function-Konsequenz im letzten Exon ohne passende Zeile in ENIGMA Tabelle 4: NMD-Escape möglich, PVS1-Stärke nach dem Entscheidungsbaum prüfen.",
                     "Loss-of-function consequence in the last exon without a matching row in ENIGMA Table 4: NMD escape possible, check PVS1 strength with the decision tree.",
                     "pathogenic")]
    return [hint("PVS1 (Kandidat)",
                 "Vorzeitiger Stopp, Leserasterverschiebung oder Spleißstelle ohne passende Zeile in ENIGMA Tabelle 4 (z. B. Indel an der Spleißstelle). Stärke nach dem PVS1-Flussdiagramm festlegen.",
                 "Premature stop, frameshift or splice site without a matching row in ENIGMA Table 4 (e.g. an indel at the splice site). Set the strength with the PVS1 flowchart.",
                 "pathogenic")]


def acmg_hints(gene, csq, proxy_af, ref="", alt="", evidence=None):
    """All hints for one variant. `evidence` holds optional gnomAD, SpliceAI and BayesDel lookups."""
    evidence = evidence or {}
    consequences = (csq.get("Consequence") or "").split("&")
    hints = []
    if any(c in LOF_CONSEQUENCES for c in consequences):
        hints += decision_tree_hints(gene, csq, consequences)
        if gene == "BRCA2" and "Lys3326Ter" in (csq.get("HGVSp") or ""):
            hints.append(hint("Bekannte Ausnahme",
                              "BRCA2 p.Lys3326Ter (c.9976A>T) gilt trotz Stoppcodon als benigne für hohes HBOC-Risiko.",
                              "BRCA2 p.Lys3326Ter (c.9976A>T) is considered benign for high HBOC risk despite the stop codon.",
                              "benign"))
    for extra in (enigma_tables.functional_hint(gene, csq.get("HGVSc")),
                  enigma_tables.multifactorial_hint(gene, csq.get("HGVSc"))):
        if extra:
            hints.append(extra)
    is_substitution = len(ref) == 1 and len(alt) == 1
    hints += frequency_hints(is_substitution, evidence.get("gnomad"), proxy_af)
    hints += predictor_hints(gene, csq, evidence.get("spliceai"), evidence.get("bayesdel"))
    if any(c in SPLICE_NEAR for c in consequences) and evidence.get("spliceai") is None:
        hints.append(hint("Spleißeffekt prüfen",
                          "Nahe einer Spleißstelle: Spleißvorhersage (SpliceAI) und gegebenenfalls RNA-Analyse nötig.",
                          "Near a splice site: splice prediction (SpliceAI) and possibly RNA analysis required.",
                          "neutral"))
    return hints

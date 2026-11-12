"""Report language. Result JSON values stay German (stable schema); HTML/Markdown are rendered in `de` or `en`."""
import re

UI = {
    "title": ("BRCA1/BRCA2-Bericht", "BRCA1/BRCA2 report"),
    "important": ("Wichtig", "Important"),
    "reference": ("Referenz", "Reference"),
    "input": ("Eingabe", "Input"),
    "snv_caller": ("SNV/Indel-Caller", "SNV/indel caller"),
    "created": ("Erstellt", "Created"),
    "strand_minus": ("Minus-Strang", "minus strand"),
    "strand_plus": ("Plus-Strang", "plus strand"),
    "small_title": ("{gene} SNVs/Indels ({chrom}, {strand}, {transcript})",
                    "{gene} SNVs/indels ({chrom}, {strand}, {transcript})"),
    "small_headers": (["POS", "REF", "ALT", "Genotyp", "Konsequenz", "Transkript", "HGVSc", "HGVSp",
                       "ClinVar", "ACMG/AMP-Hinweise", "Qualität"],
                      ["POS", "REF", "ALT", "Genotype", "Consequence", "Transcript", "HGVSc", "HGVSp",
                       "ClinVar", "ACMG/AMP hints", "Quality"]),
    "no_small": ("Keine Varianten gegenüber der Referenz gefunden.", "No variants relative to the reference."),
    "cnv_title": ("BRCA1/BRCA2 CNVs (Deletionen/Duplikationen)", "BRCA1/BRCA2 CNVs (deletions/duplications)"),
    "cnv_headers": (["Gen", "Typ", "Bereich", "CN", "log2/Ratio", "Exons", "Hinweise"],
                    ["Gene", "Type", "Region", "CN", "log2/ratio", "Exons", "Notes"]),
    "sv_title": ("BRCA1/BRCA2 strukturelle Varianten", "BRCA1/BRCA2 structural variants"),
    "sv_headers": (["Gen", "Typ", "Position", "Partner", "GT", "PR", "SR", "FILTER", "Hinweise"],
                   ["Gene", "Type", "Position", "Mate", "GT", "PR", "SR", "FILTER", "Notes"]),
    "no_sv": ("Manta hat keine SV mit Bruchpunkt in BRCA1/BRCA2 gemeldet.",
              "Manta reported no SV with a breakpoint in BRCA1/BRCA2."),
    "exon_title": ("Coverage je Exon (MAPQ >= 20)", "Coverage per exon (MAPQ >= 20)"),
    "exon_headers": (["Gen", "Exon", "Bereich", "Mittel", "% >=10x", "% >=20x", "% >=30x", "Anteil MAPQ<20",
                      "Status", "Grund"],
                     ["Gene", "Exon", "Region", "Mean", "% >=10x", "% >=20x", "% >=30x", "Fraction MAPQ<20",
                      "Status", "Reason"]),
    "low_title": ("Unsichere Bereiche (Tiefe < {depth:g}x)", "Uncertain regions (depth < {depth:g}x)"),
    "low_headers": (["Gen", "Bereich", "Länge (bp)", "min. Tiefe"], ["Gene", "Region", "Length (bp)", "Min. depth"]),
    "no_low": ("Keine Bereiche unter der Schwelle.", "No regions below the threshold."),
    "no_entries": ("Keine Einträge.", "No entries."),
    "versions": ("Tool-Versionen", "Tool versions"),
    "version": ("Version", "Version"),
    "unknown": ("unbekannt", "unknown"),
    "commands": ("Reproduzierbare Befehle", "Reproducible commands"),
    "commands_note": ("Alle ausgeführten Befehle stehen vollständig in", "Every executed command is listed verbatim in"),
    "limitations": ("Grenzen", "Limitations"),
    "na": ("n. v.", "n/a"),
}

EXACT = {
    "nicht bestimmt": "not determined", "homozygot Referenz": "homozygous reference",
    "homozygot": "homozygous", "heterozygot": "heterozygous",
    "compound/multiallelisch": "compound/multiallelic",
    "unsicher": "uncertain", "ok": "ok",
    "Deletion": "Deletion", "Duplikation": "Duplication",
    "möglicher Verlust": "possible loss", "möglicher Zugewinn": "possible gain",
    "geringe Tiefe (DP < 15)": "low depth (DP < 15)",
    "geringe Genotypqualität (GQ < 20)": "low genotype quality (GQ < 20)",
    "Allelverhältnis auffällig (mögliches Mosaik oder Artefakt)": "unusual allele balance (possible mosaic or artefact)",
    "geringe Mapping-Qualität (MQ < 40)": "low mapping quality (MQ < 40)",
    "unpräzise Bruchpunkte": "imprecise breakpoints",
    "Bestätigung mit orthogonaler Methode nötig": "confirmation with an orthogonal method required",
    "Kandidat, Bestätigung per MLPA oder validierter CNV-Analyse nötig":
        "candidate, confirmation by MLPA or validated CNV analysis required",
    "Exonweise Tiefe relativ zum Chromosomen-Median; nur ein Hinweis, kein CNV-Call.":
        "Exon depth relative to the chromosome median; a hint only, not a CNV call.",
}

PATTERNS = [
    (re.compile(r"^mittlere Tiefe (.+)x < (.+)x$"), r"mean depth \1x < \2x"),
    (re.compile(r"^nur (.+) % der Basen >= 20x$"), r"only \1% of bases >= 20x"),
    (re.compile(r"^(.+) % der Reads mit MAPQ < 20 \(Pseudogen/Repeat\)$"), r"\1% of reads with MAPQ < 20 (pseudogene/repeat)"),
    (re.compile(r"^liegt in unsicher abgedecktem Bereich (.+)$"), r"located in uncertain region \1"),
    (re.compile(r"^Exon (\d+) \(Legacy (\d+)\)$"), r"Exon \1 (legacy \2)"),
]


def ui(key, lang, **kwargs):
    value = UI[key][0 if lang == "de" else 1]
    return value.format(**kwargs) if kwargs and isinstance(value, str) else value


def tr(text, lang):
    """Translate a German result value for display; unknown values are returned unchanged."""
    if lang == "de" or not isinstance(text, str):
        return text
    if text in EXACT:
        return EXACT[text]
    for pattern, replacement in PATTERNS:
        if pattern.match(text):
            return pattern.sub(replacement, text)
    return text

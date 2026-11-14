"""Lookups in the ClinGen ENIGMA BRCA1/BRCA2 VCEP tables v1.2 (assets/enigma, see SOURCE.json).

    pvs1.tsv            Table 4: PVS1 decision trees with PM5_PTC per exon, splice sites (+-1,2),
                        start codon, whole-exon deletions and duplications
    functional.tsv      Table 9: calibrated functional assays (PS3/BS3)
    multifactorial.tsv  combined likelihood ratios from published clinical data (PP4/BP5)

All lookups use the c. notation of the MANE Select transcripts (NM_007294.4, NM_000059.4), which
VEP reports for ENST00000357654 and ENST00000380152. The tables were built with
bin/brca12_make_enigma.py from the UCSC ENIGMA track hub; check their version against the ClinGen
CSpec registry (GN092, GN097). Research use only.
"""
import csv
import os
import re

ASSETS = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                      "assets", "enigma")

_STRENGTH = {"very strong": "VeryStrong", "verystrong": "VeryStrong", "strong": "Strong",
             "moderate": "Moderate", "supporting": "Supporting"}
_DEFAULT = {"PVS": "VeryStrong", "PS": "Strong", "PM": "Moderate", "PP": "Supporting",
            "BS": "Strong", "BP": "Supporting"}

_CACHE = {}


def _read(name):
    if name not in _CACHE:
        path = os.path.join(ASSETS, name)
        if not os.path.exists(path):
            _CACHE[name] = None
        else:
            with open(path, newline="") as handle:
                _CACHE[name] = list(csv.DictReader(handle, delimiter="\t"))
    return _CACHE[name]


def available():
    return _read("pvs1.tsv") is not None


def c_part(hgvsc):
    hgvsc = (hgvsc or "").strip()
    return hgvsc.split(":", 1)[1] if ":" in hgvsc else hgvsc


def _label(name, strength):
    prefix = re.match(r"[A-Z]+", name).group(0)
    prefix = next((p for p in _DEFAULT if prefix.startswith(p)), prefix)
    return name if strength in (None, _DEFAULT.get(prefix)) else f"{name}_{strength}"


def _direction(name):
    return "benign" if name.startswith("B") else "pathogenic"


def _hint(code, de, en, direction):
    return {"code": code, "textDE": de, "textEN": en, "direction": direction}


def parse_code(text):
    """'PVS1_Strong (RNA)' -> ('PVS1', 'Strong', True); 'PM5_N/A' -> ('PM5_PTC', 'N/A', False)."""
    text = text.strip()
    rna = "(RNA)" in text
    match = re.match(r"([A-Z]+\d)(?:_(.+?))?$", text.replace("(RNA)", "").strip())
    if not match:
        return None
    name, strength = match.group(1), match.group(2)
    if name == "PM5":
        name = "PM5_PTC"          # ENIGMA repurposes PM5 for termination codons, weighted per exon
    if strength and strength.upper() != "N/A":
        strength = _STRENGTH.get(strength.lower(), strength)
    return name, strength, rna


def _coding_position(hgvsc):
    """First coding position of a c. notation: 'c.3000del' -> 3000, 'c.-19+1G>A' -> None (intronic)."""
    match = re.match(r"c\.(-?\d+)(?![+\-\d])", hgvsc or "")
    return int(match.group(1)) if match else None


def _range(position):
    match = re.match(r"c\.(-?\d+)-c\.(-?\d+)$", position)
    return (int(match.group(1)), int(match.group(2))) if match else None


def table4_row(gene, hgvsc, kind):
    """Row of Table 4 for a PTC ('ptc'), a splice-site substitution ('splice') or the start codon ('start')."""
    rows = _read("pvs1.tsv") or []
    hgvsc = c_part(hgvsc)
    if kind == "ptc":
        pos = _coding_position(hgvsc)
        if pos is None:
            return None
        for row in rows:
            span = _range(row["position"])
            if row["gene"] == gene and row["variant"] == "PTC" and span and span[0] <= pos <= span[1]:
                return row
    elif kind == "splice":
        match = re.match(r"(c\.[-*]?\d+[+-][12])[ACGT]>([ACGT])$", hgvsc)
        if not match:
            return None
        position, alt = match.groups()
        for row in rows:
            if row["gene"] == gene and re.sub(r"[ACGT]?>?$", "", row["position"]) == position:
                if alt in {a.strip() for a in row["variant"].split(",")}:
                    return row
    elif kind == "start":
        for row in rows:
            if row["gene"] == gene and row["variant"] == "Start Codon":
                return row
    return None


def table4_hints(row, what_de, what_en):
    """Hints for one Table 4 row, e.g. 'PVS1, PM5_Strong' -> PVS1 and PM5_PTC_Strong."""
    hints = []
    note = "; ".join(v for v in (row.get("observation"), row.get("warning")) if v)
    note_de = f" Hinweis der Tabelle: {note}." if note else ""
    note_en = f" Table note: {note}." if note else ""
    for part in row["code"].split(","):
        parsed = parse_code(part)
        if not parsed:
            continue
        name, strength, rna = parsed
        rna_de = " (gestützt auf RNA-Daten)" if rna else ""
        rna_en = " (based on RNA data)" if rna else ""
        if strength == "N/A":
            hints.append(_hint(f"{name} nicht anwendbar",
                               f"{name} ist nach dem ENIGMA-Entscheidungsbaum für {what_de} in Exon {row['exon']} nicht anwendbar (Tabelle 4).{note_de}",
                               f"{name} does not apply to {what_en} in exon {row['exon']} according to the ENIGMA decision tree (Table 4).{note_en}",
                               "neutral"))
            continue
        hints.append(_hint(_label(name, strength),
                           f"ENIGMA-Entscheidungsbaum, Tabelle 4: {what_de} in Exon {row['exon']}{rna_de}.{note_de}",
                           f"ENIGMA decision tree, Table 4: {what_en} in exon {row['exon']}{rna_en}.{note_en}",
                           _direction(name)))
    return hints


def functional_hint(gene, hgvsc):
    """PS3/BS3 from Table 9 (calibrated functional assays)."""
    for row in _read("functional.tsv") or []:
        if row["gene"] == gene and row["hgvsc"] == c_part(hgvsc):
            if row["code"] not in ("PS3", "BS3") or row["weight"] in ("", "N/A"):
                return _hint("PS3/BS3 nicht erfüllt",
                             f"ENIGMA Tabelle 9: {row['text']}", f"ENIGMA Table 9: {row['text']}", "neutral")
            strength = _STRENGTH.get(row["weight"].lower(), row["weight"])
            return _hint(_label(row["code"], strength), f"ENIGMA Tabelle 9: {row['text']}",
                         f"ENIGMA Table 9: {row['text']}", _direction(row["code"]))
    return None


def multifactorial_hint(gene, hgvsc):
    """PP4/BP5 from the published combined likelihood ratio."""
    for row in _read("multifactorial.tsv") or []:
        if row["gene"] == gene and row["hgvsc"] == c_part(hgvsc):
            match = re.match(r"(PP4|BP5) - \w+ - (.+)$", row["code"])
            if not match:
                return None
            name, strength = match.group(1), _STRENGTH.get(match.group(2).lower())
            lr = row["combined_lr"]
            return _hint(_label(name, strength),
                         f"Kombinierte Likelihood Ratio {lr} aus veröffentlichten klinischen Daten (ENIGMA, multifaktoriell).",
                         f"Combined likelihood ratio {lr} from published clinical data (ENIGMA, multifactorial).",
                         _direction(name))
    return None


def exon_event_row(gene, exon_number, kind, tandem=None):
    """Table 4 row for a deletion ('DEL') or duplication ('DUP') of one exon (MANE numbering)."""
    for row in _read("pvs1.tsv") or []:
        match = re.match(r"(DEL|DUP) E?(\d+)(?:\(\d+\))?\s*(Unknown|Tandem)?$", row["variant"])
        if not match or row["gene"] != gene:
            continue
        if match.group(1) == kind and int(match.group(2)) == exon_number:
            if kind == "DEL" or (match.group(3) or "Unknown") == (tandem or "Unknown"):
                return row
    return None

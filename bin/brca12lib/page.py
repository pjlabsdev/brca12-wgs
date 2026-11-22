"""HTML report page: one file with both languages, two tabs (result, details) and no script.

Language and tab are switched with hidden radio buttons and CSS only, so the file works offline, in
mail previews and with scripts blocked. Printing (and `--report_pdf`) shows the selected language with
both tabs one after the other. The font (IBM Plex Mono, SIL OFL 1.1, see assets/fonts/OFL.txt) is embedded."""
import base64
import html
import os

from . import figures
from .i18n import tr, ui

DOI = "10.5281/zenodo.24173638"
RESEARCH_URL = "https://pjlabs.dev/forschung/brca12-wgs"
FONT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assets", "fonts")

TEXT = {
    "ruo": ("Nur für Forschung (RUO). Kein medizinischer Befund und kein Medizinprodukt.",
            "Research use only (RUO). Not a medical report and not a medical device."),
    "eyebrow": ("BRCA1/BRCA2 aus dem Genom", "BRCA1/BRCA2 from the genome"),
    "tab_r": ("Ergebnis", "Result"),
    "tab_d": ("Details", "Details"),
    "sample": ("Probe", "Sample"),
    "sample_tip": ("Name der Probe (--sample).", "Sample name (--sample)."),
    "input": ("Eingabe", "Input"),
    "input_tip": ("Eingabedateien ohne Ordner. Art: {mode}.", "Input files without folder. Kind: {mode}."),
    "reference": ("Referenz", "Reference"),
    "analysis": ("Analyse", "Analysis"),
    "analysis_tip": ("Zeitpunkt des Berichts (UTC): {ts}", "Time of the report (UTC): {ts}"),
    "program": ("Programm", "Program"),
    "program_tip": ("brca12-wgs {version} · DOI {doi} · Nextflow-Sitzung {session}",
                    "brca12-wgs {version} · DOI {doi} · Nextflow session {session}"),
    "caller": ("SNV/Indel-Rufer", "SNV/indel caller"),
    "depth": ("Tiefe über den Exons", "Depth over the exons"),
    "depth_tip": ("Mittlere Tiefe aller Exons je Gen (Reads mit MAPQ >= 20). Unter {min:g}x gilt ein Exon als unsicher.",
                  "Mean depth of all exons per gene (reads with MAPQ >= 20). Below {min:g}x an exon counts as uncertain."),
    "quality": ("Qualität", "Quality"),
    "quality_tip": ("Bestanden, wenn kein Exon unsicher ist, die mittlere Exon-Tiefe beider Gene mindestens {min:g}x beträgt und die Read-Prüfung keine Warnung meldet. Sonst eingeschränkt, mit Grund.",
                    "Passed if no exon is uncertain, the mean exon depth of both genes is at least {min:g}x and the read check reports no warning. Otherwise limited, with the reason."),
    "passed": ("bestanden", "passed"),
    "limited": ("eingeschränkt", "limited"),
    "q_exons": ("unsichere Exons: {n}", "uncertain exons: {n}"),
    "q_depth": ("{gene}: mittlere Exon-Tiefe {depth:g}x", "{gene}: mean exon depth {depth:g}x"),
    "q_reads": ("Warnungen der Read-Prüfung: {n}", "read-check warnings: {n}"),
    "low_depth": ("Geringe Tiefe: {gene} hat im Mittel {depth:g}x über den Exons (Schwelle {min:g}x). Kleine Varianten, Kopienzahl und strukturelle Varianten sind dort unsicherer; ein leeres Ergebnis sagt wenig.",
                  "Low depth: {gene} has a mean of {depth:g}x over the exons (threshold {min:g}x). Small variants, copy number and structural variants are less certain there; an empty result says little."),
    "card_small": ("SNVs und Indels gegenüber der Referenz", "SNVs and indels relative to the reference"),
    "card_small_tip": ("Alle kleinen Varianten in BRCA1 und BRCA2, auch häufige und harmlose.",
                       "All small variants in BRCA1 and BRCA2, including common and benign ones."),
    "card_plp": ("ClinVar pathogen oder wahrscheinlich pathogen", "ClinVar pathogenic or likely pathogenic"),
    "card_plp_tip": ("Varianten mit ClinVar-Angabe pathogen oder wahrscheinlich pathogen, ohne widersprüchliche Einträge. Dazu {lof} Kandidaten für Funktionsverlust (PVS1-Hinweis).",
                     "Variants with a ClinVar entry pathogenic or likely pathogenic, without conflicting entries. In addition {lof} loss-of-function candidates (PVS1 hint)."),
    "card_struct": ("Kandidaten für CNV und strukturelle Varianten", "Candidates for CNVs and structural variants"),
    "card_struct_tip": ("CNVkit-Kandidaten {cnv}, Exon-Kopienzahl (Kohorte) {exon}, SV-Rufer {sv}, Einschub-Signale {ins}. Schwächere Hinweise stehen in den Details.",
                        "CNVkit candidates {cnv}, exon copy number (cohort) {exon}, SV caller {sv}, insertion signals {ins}. Weaker hints are listed in the details."),
    "card_exons": ("unsichere Exons", "uncertain exons"),
    "card_exons_tip": ("Exons mit zu geringer Tiefe, zu wenig Basen ab 20x oder vielen mehrdeutig kartierten Reads, von {total} Exons.",
                       "Exons with too little depth, too few bases at 20x or many ambiguously mapped reads, out of {total} exons."),
    "notable": ("Zu prüfen", "To review"),
    "notable_headers": (["Gen", "Art", "Angabe", "Genotyp", "ClinVar", "Hinweise"],
                        ["Gene", "Kind", "Finding", "Genotype", "ClinVar", "Notes"]),
    "notable_none": ("Keine Einträge in den geprüften Kategorien (ClinVar pathogen oder wahrscheinlich pathogen, unklare Bedeutung, Funktionsverlust, CNV, strukturelle Varianten, Mosaik). Das schließt eine erbliche Veranlagung nicht aus, siehe Grenzen in den Details.",
                     "No entries in the checked categories (ClinVar pathogenic or likely pathogenic, uncertain significance, loss of function, CNV, structural variants, mosaic). This does not rule out a hereditary predisposition, see the limitations in the details."),
    "notable_note": ("Maschinell zusammengestellt, keine Klassifikation. Jede Angabe braucht eine Bestätigung in einem akkreditierten Labor.",
                     "Compiled automatically, not a classification. Every entry needs confirmation in an accredited laboratory."),
    "kind_plp": ("ClinVar pathogen / wahrscheinlich pathogen", "ClinVar pathogenic / likely pathogenic"),
    "kind_lof": ("Kandidat Funktionsverlust", "Loss-of-function candidate"),
    "kind_vus": ("ClinVar unklare Bedeutung", "ClinVar uncertain significance"),
    "kind_mosaic": ("Mosaik-Kandidat", "Mosaic candidate"),
    "kind_cnv": ("CNV-Kandidat (CNVkit)", "CNV candidate (CNVkit)"),
    "kind_exon": ("Exon-Kopienzahl (Kohorte)", "Exon copy number (cohort)"),
    "kind_sv": ("Strukturelle Variante", "Structural variant"),
    "kind_ins": ("Einschub-Signal", "Insertion signal"),
    "coverage": ("Abdeckung", "Coverage"),
    "cov_headers": (["Gen", "Exons", "mittlere Tiefe", "geringste Exon-Tiefe", "unsichere Exons"],
                    ["Gene", "Exons", "Mean depth", "Lowest exon depth", "Uncertain exons"]),
    "limits": ("Grenzen", "Limitations"),
    "cite": ("Zitieren", "Cite"),
    "research": ("Forschungsbeitrag", "Research article"),
    "footer": ("Erstellt mit brca12-wgs {version}.", "Created with brca12-wgs {version}."),
    "page": ("BRCA1/BRCA2-Bericht", "BRCA1/BRCA2 report"),
}

CSS = """
:root{--bg:#FAF9F5;--surface:#FFFFFF;--fg:#141413;--mut:#3F3E38;--soft-fg:#6E6C63;--line:#D8D5C7;--hair:rgba(216,213,199,.78);--acc:#D97757;--acc-t:#B34D30;--ok:#248A3D;--warn:#B34D30;--soft:#F2F0E8;--mark:#FBEFE6;--shadow:0 1px 2px rgba(20,20,19,.04),0 10px 28px rgba(20,20,19,.055);--r:14px;--ease:cubic-bezier(.22,.61,.36,1)}
@media (prefers-color-scheme:dark){:root{--bg:#1C1C1E;--surface:#2C2C2E;--fg:#FFFFFF;--mut:rgba(235,235,245,.72);--soft-fg:rgba(235,235,245,.55);--line:#38383A;--hair:rgba(255,255,255,.1);--acc-t:#E38B69;--ok:#30D158;--warn:#E38B69;--soft:#242426;--mark:#3A2A20;--shadow:0 1px 2px rgba(0,0,0,.28),0 12px 30px rgba(0,0,0,.3)}}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--fg);font:400 14px/1.6 'IBM Plex Mono',ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;-webkit-font-smoothing:antialiased}
input.sw{position:absolute;opacity:0;pointer-events:none}
.bar{position:sticky;top:0;z-index:5;background:color-mix(in srgb,var(--bg) 82%,transparent);-webkit-backdrop-filter:saturate(180%) blur(18px);backdrop-filter:saturate(180%) blur(18px);border-bottom:1px solid var(--hair)}
.bar-in{max-width:1100px;margin:0 auto;padding:12px 24px;display:flex;align-items:center;justify-content:space-between;gap:16px}
.brand{font-weight:500;letter-spacing:.02em;font-size:13px}.brand i{font-style:normal;color:var(--acc-t)}
.seg{display:inline-flex;background:var(--soft);border:1px solid var(--hair);border-radius:999px;padding:3px;gap:2px}
.seg label{font-size:12px;font-weight:500;border-radius:999px;color:var(--soft-fg);padding:5px 14px;cursor:pointer;transition:background 140ms var(--ease),color 140ms var(--ease)}
.seg label:hover{color:var(--fg)}
#lang-de:checked~header label[for=lang-de],#lang-en:checked~header label[for=lang-en],#tab-r:checked~main label[for=tab-r],#tab-d:checked~main label[for=tab-d]{background:var(--surface);color:var(--fg);box-shadow:0 1px 2px rgba(20,20,19,.12)}
#lang-de:focus-visible~header label[for=lang-de],#lang-en:focus-visible~header label[for=lang-en],#tab-r:focus-visible~main label[for=tab-r],#tab-d:focus-visible~main label[for=tab-d]{outline:2px solid var(--acc);outline-offset:2px}
#lang-de:checked~main.l-en,#lang-en:checked~main.l-de{display:none}
#tab-r:checked~main section.tab-d,#tab-d:checked~main section.tab-r{display:none}
main{max-width:1100px;margin:0 auto;padding:40px 24px 64px}
.ruo{font-size:11px;font-weight:500;letter-spacing:.06em;text-transform:uppercase;color:var(--acc-t);border:1px solid color-mix(in srgb,var(--acc) 40%,transparent);border-radius:999px;display:inline-block;padding:3px 12px;margin:0 0 22px}
.eyebrow{font-size:12px;font-weight:500;letter-spacing:.1em;text-transform:uppercase;color:var(--acc-t);margin-bottom:8px}
h1{font-size:30px;line-height:1.2;font-weight:500;margin:0;letter-spacing:-.01em;overflow-wrap:anywhere}
h2{font-size:12px;font-weight:500;letter-spacing:.1em;text-transform:uppercase;color:var(--soft-fg);margin:40px 0 12px}
nav.tabs{margin:28px 0 4px}
.head{display:grid;grid-template-columns:repeat(4,1fr);margin-top:22px;background:var(--surface);border:1px solid var(--hair);border-radius:var(--r);box-shadow:var(--shadow)}
.head div{padding:12px 16px;border-right:1px solid var(--hair);border-bottom:1px solid var(--hair);min-width:0}.head div:nth-child(4n){border-right:0}.head div:nth-last-child(-n+4){border-bottom:0}
.head span.k{display:block;font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--soft-fg)}.head b{display:block;font-weight:500;font-size:13px;margin-top:3px;overflow-wrap:anywhere}
.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-top:28px}
.card{background:var(--surface);border:1px solid var(--hair);border-radius:var(--r);padding:18px 18px 16px;box-shadow:var(--shadow)}
.card b{display:block;font-size:26px;line-height:1.15;font-weight:500;letter-spacing:-.01em}
.card span.k{display:block;margin-top:6px;color:var(--soft-fg);font-size:12px;line-height:1.4}
.panel{background:var(--surface);border:1px solid var(--hair);border-radius:var(--r);box-shadow:var(--shadow);padding:6px 20px;overflow-x:auto}
table{width:100%;border-collapse:collapse;font-size:12px}
th{font-weight:500;text-align:left;color:var(--soft-fg);font-size:11px;letter-spacing:.06em;text-transform:uppercase;border-bottom:1px solid var(--line);padding:12px 12px 10px 0}
td{border-bottom:1px solid var(--hair);padding:9px 12px 9px 0;vertical-align:top;overflow-wrap:break-word}
tr:last-child td{border-bottom:0}
tr.uncertain td{background:var(--mark);-webkit-print-color-adjust:exact;print-color-adjust:exact}
.note{color:var(--soft-fg);font-size:12px;line-height:1.55;margin:10px 0}
.ok{color:var(--ok)}.warn{color:var(--warn)}
.alert{margin:22px 0 0;padding:12px 16px;border-radius:var(--r);background:color-mix(in srgb,var(--acc) 12%,var(--surface));border:1px solid color-mix(in srgb,var(--acc) 40%,transparent);font-size:12px;line-height:1.55}
.alert.soft{background:var(--surface);border-color:var(--hair)}
.prose{padding:8px 20px}.prose li,.prose p{margin:8px 0;color:var(--mut);font-size:12px;line-height:1.65}.prose ul{padding-left:18px}
pre{background:var(--soft);padding:12px 14px;border-radius:10px;white-space:pre-wrap;overflow-wrap:anywhere;font:400 12px/1.6 'IBM Plex Mono',ui-monospace,monospace;margin:8px 0}
a{color:var(--acc-t);text-decoration:none}a:hover{text-decoration:underline}
.report-figure{margin:0 0 14px;padding:12px 14px;background:#fff;color:#1d1d1f;border:1px solid var(--hair);border-radius:var(--r);break-inside:avoid}
.report-figure svg{max-width:100%;height:auto}
.figure-note{font-size:12px;color:#636366;margin:4px 0 0}
.tip{position:relative;cursor:help;text-decoration:underline dotted color-mix(in srgb,var(--soft-fg) 70%,transparent);text-underline-offset:3px;text-decoration-thickness:1px;outline:none}
.tip:hover::after,.tip:focus::after{content:attr(data-tip);position:absolute;left:0;top:calc(100% + 8px);z-index:50;width:max-content;max-width:min(380px,80vw);background:var(--fg);color:var(--bg);font:400 11.5px/1.55 'IBM Plex Mono',ui-monospace,monospace;letter-spacing:0;text-transform:none;white-space:normal;padding:10px 12px;border-radius:10px;box-shadow:0 8px 28px rgba(20,20,19,.22)}
.head div:nth-child(4n) .tip::after,.head div:nth-child(4n+3) .tip::after,.card:nth-child(n+3) .tip::after{left:auto;right:0}
footer{margin-top:48px;border-top:1px solid var(--line);padding-top:18px;color:var(--soft-fg);font-size:12px}
footer p{margin:6px 0}.cite{color:var(--fg)}
@media (max-width:720px){.cards,.head{grid-template-columns:repeat(2,1fr)}.head div:nth-child(4n){border-right:1px solid var(--hair)}.head div:nth-child(2n){border-right:0}.head div:nth-last-child(-n+4){border-bottom:1px solid var(--hair)}.head div:nth-last-child(-n+2){border-bottom:0}h1{font-size:24px}.panel{padding-left:14px;padding-right:14px}.head div:nth-child(2n) .tip::after,.card:nth-child(2n) .tip::after{left:auto;right:0}}
@page{size:A4;margin:16mm 12mm 18mm;@bottom-right{content:counter(page) ' / ' counter(pages);font:400 9px 'IBM Plex Mono',monospace;color:#6E6C63}}
@media print{:root{--bg:#fff;--surface:#fff;--shadow:none}.bar,nav.tabs{display:none}main{padding:0;max-width:none}body{font-size:10px}table{font-size:9px}a{color:inherit}
#tab-r:checked~main section.tab-d,#tab-d:checked~main section.tab-r{display:block}
tr{break-inside:avoid}thead{display:table-header-group}h2{break-after:avoid}.card,.panel,.head{box-shadow:none}.panel{overflow:visible}.tip{text-decoration:none}.tip::after{display:none!important}}
"""


def _t(key, lang, **kwargs):
    value = TEXT[key][0 if lang == "de" else 1]
    return value.format(**kwargs) if kwargs else value


def _fonts():
    """@font-face rules with the embedded font; empty when the font files are missing (system font is used)."""
    rules = []
    for weight in (400, 500):
        path = os.path.join(FONT_DIR, f"IBMPlexMono-Latin-{weight}.woff2")
        try:
            with open(path, "rb") as handle:
                data = base64.b64encode(handle.read()).decode("ascii")
        except OSError:
            return ""
        rules.append(f"@font-face{{font-family:'IBM Plex Mono';font-weight:{weight};font-style:normal;font-display:swap;"
                     f"src:url(data:font/woff2;base64,{data}) format('woff2')}}")
    return "\n".join(rules)


def _css_text(text):
    """Text that is safe inside a CSS string in a <style> block."""
    return "".join(c if c.isalnum() or c in " ._-" else "_" for c in str(text))


def _tip(label, tip):
    return f'<span class="tip" tabindex="0" data-tip="{html.escape(tip, quote=True)}">{html.escape(label)}</span>'


def _gene_depth(result):
    """{gene: (exon count, mean of exon mean depths, lowest exon mean depth, uncertain exons)}"""
    out = {}
    for gene in ("BRCA1", "BRCA2"):
        exons = [e for e in result["exonCoverage"] if e["gene"] == gene and e.get("meanDepth") is not None]
        if exons:
            depths = [e["meanDepth"] for e in exons]
            out[gene] = (len(exons), round(sum(depths) / len(depths), 1), min(depths),
                         sum(1 for e in exons if e["status"] != "ok"))
    return out


def quality(result):
    """(passed, [reason keys with arguments]) for the report header; language-neutral."""
    minimum = result["thresholds"]["minExonDepth"]
    reasons = []
    uncertain = sum(1 for e in result["exonCoverage"] if e["status"] != "ok")
    if uncertain:
        reasons.append(("q_exons", {"n": uncertain}))
    for gene, (_, mean, _, _) in _gene_depth(result).items():
        if mean < minimum:
            reasons.append(("q_depth", {"gene": gene, "depth": mean}))
    flags = len((result.get("readQC") or {}).get("flags", []))
    if flags:
        reasons.append(("q_reads", {"n": flags}))
    return not reasons, reasons


def _is_plp(v):
    sig = v["clinvar"]["significance"] or ""
    return "athogenic" in sig and "onflicting" not in sig


def _is_vus(v):
    sig = v["clinvar"]["significance"] or ""
    return "ncertain" in sig and "onflicting" not in sig


def notable_rows(result, lang):
    """Rows of the 'to review' table on the result tab. Only entries the pipeline itself reports as
    candidates; weaker hints (exon dosage, cohort hints below the call quality) stay in the details."""
    rows = []
    hint_key = "textDE" if lang == "de" else "textEN"

    def small(v, kind):
        name = " ".join(x for x in (v.get("hgvsc"), v.get("hgvsp")) if x) or f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}"
        notes = [f"{h['code']}: {h[hint_key]}" for h in v.get("acmgHints", [])] + [tr(f, lang) for f in v.get("qualityFlags", [])]
        return [v["gene"], _t(kind, lang), name, tr(v.get("genotype"), lang), v["clinvar"]["significance"], "; ".join(notes)]

    for gene in ("BRCA1", "BRCA2"):
        for v in result["genes"][gene]["smallVariants"]:
            lof = any(h["code"].startswith("PVS1") for h in v.get("acmgHints", []))
            if _is_plp(v):
                rows.append(small(v, "kind_plp"))
            elif lof:
                rows.append(small(v, "kind_lof"))
            elif _is_vus(v):
                rows.append(small(v, "kind_vus"))
    for v in result.get("mosaicCandidates", []):
        row = small(v, "kind_mosaic")
        row[3] = f"AF {v.get('alleleFraction')}"
        rows.append(row)

    def exons(item):
        return ", ".join(str(tr(e, lang)) for e in item.get("affectedExons", []))

    def flags(item):
        return "; ".join(str(tr(f, lang)) for f in item.get("flags", []))

    for c in result["cnvs"]:
        rows.append([c["gene"], _t("kind_cnv", lang),
                     f"{tr(c['type'], lang)} {c['chrom']}:{c['start']}-{c['end']} · CN {c['copyNumber']} · {exons(c)}",
                     None, None, flags(c)])
    for c in (result.get("exonCopyNumber") or {}).get("calls", []):
        rows.append([c["gene"], _t("kind_exon", lang),
                     f"{tr(c['type'], lang)} {c['chrom']}:{c['start']}-{c['end']} · CN {c['copyNumber']} · {exons(c)}",
                     None, None, flags(c)])
    for s in result["svs"]:
        pos = f"{s['chrom']}:{s['pos']}" + (f"-{s['end']}" if s.get("end") else "")
        rows.append([s["gene"], _t("kind_sv", lang), f"{tr(s['type'], lang)} {pos} · {s.get('filter')}",
                     tr(s.get("genotype"), lang), None, flags(s)])
    for i in result.get("insertionSignals", []):
        rows.append([i["gene"], _t("kind_ins", lang), f"{tr(i['type'], lang)} {i['chrom']}:{i['pos']}-{i['end']}",
                     None, None, flags(i)])
    na = ui("na", lang)
    return [[na if c is None or c == "" else str(c) for c in row] for row in rows]


def _table_html(headers, rows, mark=None):
    esc = html.escape
    out = ["<table><thead><tr>" + "".join(f"<th>{esc(str(h))}</th>" for h in headers) + "</tr></thead><tbody>"]
    for row in rows:
        css = " class='uncertain'" if mark and mark in row else ""
        out.append(f"<tr{css}>" + "".join(f"<td>{esc(str(c))}</td>" for c in row) + "</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def _main(result, lang, sections, disclaimer, limitations):
    esc = html.escape
    t = lambda key, **kw: _t(key, lang, **kw)  # noqa: E731
    pipeline = result.get("pipeline") or {}
    version = pipeline.get("version") or ui("unknown", lang)
    doi = pipeline.get("doi") or DOI
    minimum = result["thresholds"]["minExonDepth"]
    depth = _gene_depth(result)
    summary = result.get("summary") or {}
    passed, reasons = quality(result)
    date = (result.get("generatedAt") or "")[:10]
    if lang == "de" and len(date) == 10:
        date = f"{date[8:10]}.{date[5:7]}.{date[0:4]}"
    quality_text = (f"<span class='ok'>{esc(t('passed'))}</span>" if passed else
                    f"<span class='warn'>{esc(t('limited'))}</span>: "
                    + esc("; ".join(t(key, **args) for key, args in reasons)))
    depth_text = " · ".join(f"{gene} {values[1]:g}x" for gene, values in depth.items()) or ui("na", lang)
    cells = [
        (_tip(t("sample"), t("sample_tip")), esc(result["sample"])),
        (_tip(t("input"), t("input_tip", mode=result["input"]["mode"])), esc(", ".join(result["input"]["files"]))),
        (_tip(t("reference"), result["reference"]["path"]), esc(str(result["reference"]["assembly"]))),
        (_tip(t("analysis"), t("analysis_tip", ts=result.get("generatedAt"))), esc(date)),
        (_tip(t("program"), t("program_tip", version=version, doi=doi, session=result.get("nextflowSession"))),
         esc(f"brca12-wgs {version}")),
        (esc(t("caller")), esc(result["snvCaller"])),
        (_tip(t("depth"), t("depth_tip", min=minimum)), esc(depth_text)),
        (_tip(t("quality"), t("quality_tip", min=minimum)), quality_text),
    ]
    parts = [f"<main lang='{lang}' class='l-{lang}'>",
             f"<div class='ruo'>{esc(t('ruo'))}</div>",
             f"<div class='eyebrow'>{esc(t('eyebrow'))}</div><h1>{esc(ui('title', lang))} · {esc(result['sample'])}</h1>",
             "<div class='head'>" + "".join(f"<div><span class='k'>{k}</span><b>{v}</b></div>" for k, v in cells) + "</div>",
             f"<nav class='tabs seg'><label for='tab-r'>{esc(t('tab_r'))}</label>"
             f"<label for='tab-d'>{esc(t('tab_d'))}</label></nav>",
             "<section class='tab-r'>"]
    structural = (summary.get("cnvCandidateCount", len(result["cnvs"])) + summary.get("exonCopyNumberCallCount", 0)
                  + summary.get("svCount", len(result["svs"])) + summary.get("insertionSignalCount", 0))
    cards = [
        (summary.get("smallVariantCount", 0), _tip(t("card_small"), t("card_small_tip"))),
        (summary.get("clinvarPathogenicCount", 0),
         _tip(t("card_plp"), t("card_plp_tip", lof=summary.get("lossOfFunctionCandidateCount", 0)))),
        (structural, _tip(t("card_struct"), t("card_struct_tip", cnv=summary.get("cnvCandidateCount", len(result["cnvs"])),
                                               exon=summary.get("exonCopyNumberCallCount", 0),
                                               sv=summary.get("svCount", len(result["svs"])),
                                               ins=summary.get("insertionSignalCount", 0)))),
        (sum(v[3] for v in depth.values()), _tip(t("card_exons"), t("card_exons_tip", total=sum(v[0] for v in depth.values())))),
    ]
    parts.append("<div class='cards'>" + "".join(f"<div class='card'><b>{n}</b><span class='k'>{label}</span></div>"
                                                 for n, label in cards) + "</div>")
    for gene, values in depth.items():
        if values[1] < minimum:
            parts.append(f"<div class='alert'>{esc(t('low_depth', gene=gene, depth=values[1], min=minimum))}</div>")
    for flag in (result.get("readQC") or {}).get("flags", []):
        parts.append(f"<div class='alert'>{esc(str(tr(flag, lang)))}</div>")
    parts.append(f"<div class='alert soft'><strong>{esc(ui('important', lang))}:</strong> {esc(disclaimer)}</div>")
    parts.append(f"<h2>{esc(t('notable'))}</h2>")
    rows = notable_rows(result, lang)
    if rows:
        parts.append("<div class='panel'>" + _table_html(t("notable_headers"), rows) + "</div>"
                     f"<p class='note'>{esc(t('notable_note'))}</p>")
    else:
        parts.append(f"<div class='panel prose'><p>{esc(t('notable_none'))}</p></div>")
    if depth:
        parts.append(f"<h2>{esc(t('coverage'))}</h2><div class='panel'>"
                     + _table_html(t("cov_headers"), [[gene, v[0], f"{v[1]:g}x", f"{v[2]:g}x", v[3]] for gene, v in depth.items()])
                     + "</div>")
    parts.append("</section><section class='tab-d'>")
    figs = figures.figures_html(result, lang)
    if figs:
        parts.append(figs)
    uncertain = tr("unsicher", lang)
    for title, table, empty in sections:
        parts.append(f"<h2>{esc(title)}</h2>")
        if not table["rows"]:
            parts.append(f"<div class='panel prose'><p>{esc(empty)}</p></div>")
            continue
        parts.append("<div class='panel'>" + _table_html(table["headers"], table["rows"], uncertain) + "</div>")
    parts.append(f"<h2>{esc(ui('versions', lang))}</h2><div class='panel'>"
                 + _table_html(["Tool", ui("version", lang)],
                               [[k, v or ui("unknown", lang)] for k, v in result["toolVersions"].items()]) + "</div>")
    parts.append(f"<h2>{esc(ui('commands', lang))}</h2>"
                 f"<p class='note'>{esc(ui('commands_note', lang))} <code>brca12_commands.sh</code>.</p>"
                 f"<pre>{esc(result['invocation'])}</pre>")
    parts.append(f"<h2>{esc(t('limits'))}</h2><div class='panel prose'><ul>"
                 + "".join(f"<li>{esc(item)}</li>" for item in limitations) + "</ul></div>")
    parts.append("</section>")
    parts.append(f"<footer><p>{esc(t('footer', version=version))} {esc(disclaimer)}</p>"
                 f"<p>{esc(t('cite'))}: <span class='cite'>brca12-wgs. Zenodo. "
                 f"<a href='https://doi.org/{esc(doi)}'>https://doi.org/{esc(doi)}</a></span> · "
                 f"{esc(t('research'))}: <a href='{RESEARCH_URL}'>{RESEARCH_URL.replace('https://', '')}</a></p></footer></main>")
    return "\n".join(parts)


def html_page(result, lang, sections_for, disclaimer_for, limitations_for):
    """Whole page. `lang` is the language shown first; the other one is one click away."""
    esc = html.escape
    first = "de" if lang == "de" else "en"
    checked = {code: " checked" if code == first else "" for code in ("de", "en")}
    parts = [f"<!DOCTYPE html><html lang='{first}'><head><meta charset='utf-8'>"
             "<meta name='viewport' content='width=device-width, initial-scale=1'>"
             f"<title>{esc(ui('title', first))} · {esc(result['sample'])} · RUO</title>",
             f"<style>{_fonts()}{CSS}"
             f"@page{{@bottom-left{{content:'{_css_text(result['sample'])} · brca12-wgs · RUO';"
             "font:400 9px 'IBM Plex Mono',monospace;color:#6E6C63}}</style></head><body>",
             f"<input class='sw' type='radio' name='lang' id='lang-de' aria-label='Deutsch'{checked['de']}>"
             f"<input class='sw' type='radio' name='lang' id='lang-en' aria-label='English'{checked['en']}>"
             "<input class='sw' type='radio' name='tab' id='tab-r' aria-label='Ergebnis / Result' checked>"
             "<input class='sw' type='radio' name='tab' id='tab-d' aria-label='Details'>",
             "<header class='bar'><div class='bar-in'><div class='brand'>PJ Labs <i>·</i> brca12-wgs</div>"
             "<nav class='lang seg'><label for='lang-de'>Deutsch</label><label for='lang-en'>English</label></nav>"
             "</div></header>"]
    for code in (first, "en" if first == "de" else "de"):
        parts.append(_main(result, code, sections_for(result, code), disclaimer_for(code), limitations_for(result, code)))
    parts.append("</body></html>")
    return "\n".join(parts)

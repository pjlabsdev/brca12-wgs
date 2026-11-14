"""Standalone report as Markdown and HTML, in German (`de`) or English (`en`)."""
import html

from .i18n import tr, ui

DISCLAIMER_DE = ("Nur für Forschung und persönliche Information (RUO). Diese Ergebnisse dürfen nicht ohne "
                 "klinisch validierte Bestätigung (akkreditiertes Labor, z. B. Sanger-Sequenzierung oder MLPA) "
                 "medizinisch interpretiert werden. Die ACMG/AMP-Angaben sind maschinelle Hinweise auf einzelne "
                 "Kriterien und keine Klassifikation. Ein unauffälliges Ergebnis schließt eine erbliche "
                 "BRCA1/BRCA2-Veranlagung nicht aus.")
DISCLAIMER_EN = ("Research use only (RUO). These results must not be interpreted medically without clinically "
                 "validated confirmation (accredited laboratory, e.g. Sanger sequencing or MLPA). The ACMG/AMP "
                 "entries are automated hints on individual criteria, not a classification. An unremarkable "
                 "result does not rule out a hereditary BRCA1/BRCA2 predisposition.")


def disclaimer(lang):
    return DISCLAIMER_DE if lang == "de" else DISCLAIMER_EN


def _table(headers, rows, lang):
    na = ui("na", lang)
    return {"headers": headers,
            "rows": [[na if c is None or c == "" else str(tr(c, lang)) for c in row] for row in rows]}


def _join(items, lang):
    return "; ".join(tr(i, lang) for i in items)


def sections(result, lang="en"):
    out = []
    for gene in ("BRCA1", "BRCA2"):
        info = result["genes"][gene]
        hint_key = "textDE" if lang == "de" else "textEN"
        rows = [[v["pos"], v["ref"], v["alt"], v["genotype"], v["consequence"], v["transcript"],
                 v["hgvsc"], v["hgvsp"], v["clinvar"]["significance"],
                 "; ".join(f"{h['code']}: {h[hint_key]}" for h in v["acmgHints"]),
                 _join(v["qualityFlags"], lang) or "ok"] for v in info["smallVariants"]]
        strand = ui("strand_minus" if info["strand"] == "-" else "strand_plus", lang)
        out.append((ui("small_title", lang, gene=gene, chrom=info["chrom"], strand=strand,
                       transcript=info["transcript"]),
                    _table(ui("small_headers", lang), rows, lang), ui("no_small", lang)))
    if result.get("mosaicCalling", {}).get("enabled"):
        mc = result["mosaicCalling"]
        rows = [[v["gene"], v["pos"], v["ref"], v["alt"], v["alleleFraction"], v["altDepth"], v["depth"],
                 v["consequence"], v["hgvsc"], v["hgvsp"], v["clinvar"]["significance"],
                 _join(v["qualityFlags"], lang)] for v in result.get("mosaicCandidates", [])]
        out.append((ui("mosaic_title", lang, lo=100 * mc["minAlleleFraction"], hi=100 * mc["maxAlleleFraction"]),
                    _table(ui("mosaic_headers", lang), rows, lang), ui("no_mosaic", lang)))
    cnv_rows = [[c["gene"], c["type"], f"{c['chrom']}:{c['start']}-{c['end']}", c["copyNumber"], c["log2"],
                 ", ".join(tr(e, lang) for e in c["affectedExons"]), _join(c["flags"], lang)]
                for c in result["cnvs"]]
    cnv_rows += [[d["gene"], d["type"], tr(d["exon"], lang), None, d["ratio"], tr(d["exon"], lang),
                  _join([d["note"]] + d.get("flags", []), lang)]
                 for d in result["exonDosageHints"]]
    out.append((ui("cnv_title", lang), _table(ui("cnv_headers", lang), cnv_rows, lang), ui("no_entries", lang)))
    sv_rows = [[s["gene"], s["type"], f"{s['chrom']}:{s['pos']}" + (f"-{s['end']}" if s["end"] else ""),
                f"{s['mateChrom']}:{s['matePos']}" if s["mateChrom"] else None, s["genotype"],
                s["pairedReads"], s["splitReads"], s["filter"], _join(s["flags"], lang)] for s in result["svs"]]
    sv_rows += [[i["gene"], i["type"], f"{i['chrom']}:{i['pos']}-{i['end']}", None, None, None,
                 f"{i['clippedReadsLeftJunction']}+{i['clippedReadsRightJunction']}", None, _join(i["flags"], lang)]
                for i in result.get("insertionSignals", [])]
    out.append((ui("sv_title", lang), _table(ui("sv_headers", lang), sv_rows, lang), ui("no_sv", lang)))
    exon_rows = [[e["gene"], e["label"], f"{e['chrom']}:{e['start']}-{e['end']}", e["meanDepth"],
                  e["pct10x"], e["pct20x"], e["pct30x"], e["lowMapqFraction"], e["status"],
                  _join(e["reasons"], lang)] for e in result["exonCoverage"]]
    out.append((ui("exon_title", lang), _table(ui("exon_headers", lang), exon_rows, lang), ui("no_entries", lang)))
    qc = result.get("readQC")
    if qc:
        rs = qc.get("regionStats", {})
        items = [("qc_reads", rs.get("readsInRegions")), ("qc_mapped", rs.get("mappedPct")),
                 ("qc_paired", rs.get("properlyPairedPct")), ("qc_dup", rs.get("duplicatePct")),
                 ("qc_mapq0", rs.get("mapq0Pct")), ("qc_insert", rs.get("medianInsertSize")),
                 ("qc_error", None if rs.get("errorRate") is None else round(100 * rs["errorRate"], 3)),
                 ("qc_readlen", rs.get("averageReadLength")), ("qc_basequal", rs.get("averageBaseQuality"))]
        ts = qc.get("targetSummary") or {}
        if ts:
            items += [("qc_targets", ts.get("targetCount")), ("qc_target_mean", ts.get("meanDepth")),
                      ("qc_target_pct20", ts.get("pct20x")), ("qc_uniformity", ts.get("uniformityPct")),
                      ("qc_fold80", ts.get("fold80BasePenalty")), ("qc_ontarget", ts.get("onTargetReadPct")),
                      ("qc_weak", ts.get("weakTargetCount"))]
        clip = qc.get("primerClipping")
        if clip:
            items += [("qc_clipped", clip.get("clippedPct")), ("qc_short", clip.get("readsFilteredTooShort"))]
        rows = [[ui(k, lang), v] for k, v in items]
        rows += [[ui("qc_flag", lang), f] for f in qc.get("flags", [])]
        out.append((ui("qc_title", lang, assay=result.get("assay", "wgs")),
                    _table(ui("qc_headers", lang), rows, lang), ui("no_entries", lang)))
        weak = [t for t in qc.get("targets", []) if t["status"] != "ok"]
        if qc.get("targets"):
            rows = [[t["name"], f"{t['chrom']}:{t['start']}-{t['end']}", t["meanDepth"], t["minDepth"], t["pct20x"]]
                    for t in weak]
            out.append((ui("weak_title", lang, depth=ts.get("minTargetDepth", 30)),
                        _table(ui("weak_headers", lang), rows, lang), ui("no_weak", lang)))
    low_rows = [[i["gene"], f"{i['chrom']}:{i['start']}-{i['end']}", i["length"], i["minDepth"]]
                for i in result["lowCoverageIntervals"]]
    out.append((ui("low_title", lang, depth=result["thresholds"]["lowDepth"]),
                _table(ui("low_headers", lang), low_rows, lang), ui("no_low", lang)))
    return out


def _limitations(result, lang):
    return result["limitations"] if lang == "de" else result.get("limitationsEN", result["limitations"])


def markdown(result, lang="en"):
    lines = [f"# {ui('title', lang)} · {result['sample']}", "",
             f"> **{ui('important', lang)}:** {disclaimer(lang)}", "",
             f"- {ui('reference', lang)}: {result['reference']['path']} ({result['reference']['assembly']})",
             f"- {ui('input', lang)}: {result['input']['mode']} · {', '.join(result['input']['files'])}",
             f"- {ui('snv_caller', lang)}: {result['snvCaller']}",
             f"- {ui('created', lang)}: {result['generatedAt']}", ""]
    for title, table, empty in sections(result, lang):
        lines += [f"## {title}", ""]
        if not table["rows"]:
            lines += [empty, ""]
            continue
        lines.append("| " + " | ".join(table["headers"]) + " |")
        lines.append("|" + "---|" * len(table["headers"]))
        for row in table["rows"]:
            lines.append("| " + " | ".join(c.replace("|", "/") for c in row) + " |")
        lines.append("")
    lines += [f"## {ui('versions', lang)}", "", f"| Tool | {ui('version', lang)} |", "|---|---|"]
    lines += [f"| {k} | {v or ui('unknown', lang)} |" for k, v in result["toolVersions"].items()]
    lines += ["", f"## {ui('commands', lang)}", "",
              f"{ui('commands_note', lang)} `brca12_commands.sh`.", "",
              "```bash", result["invocation"], "```", "", f"## {ui('limitations', lang)}", ""]
    lines += [f"- {item}" for item in _limitations(result, lang)]
    lines += ["", f"_{disclaimer(lang)}_", ""]
    return "\n".join(lines)


CSS = """
body{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px;color:#1d1d1f;margin:32px;max-width:1200px}
h1{font-weight:500;font-size:20px} h2{font-weight:500;font-size:15px;margin-top:28px}
.disclaimer{border:1px solid #c7c7cc;border-left:4px solid #ff9500;background:#fff8ec;padding:10px 14px;border-radius:8px}
table{border-collapse:collapse;width:100%;margin:8px 0} th,td{border-bottom:1px solid #e5e5ea;padding:4px 6px;text-align:left;vertical-align:top}
th{color:#6e6e73;font-weight:500} .uncertain{background:#fff4e5} pre{background:#f5f5f7;padding:10px;border-radius:8px;white-space:pre-wrap}
@media (prefers-color-scheme: dark){body{background:#1c1c1e;color:#f2f2f7} .disclaimer{background:#3a2a10;border-color:#48484a}
th,td{border-color:#38383a} .uncertain{background:#3a2a10} pre{background:#2c2c2e}}
"""


def html_report(result, lang="en"):
    esc = html.escape
    uncertain = tr("unsicher", lang)
    parts = [f"<!DOCTYPE html><html lang='{lang}'><head><meta charset='utf-8'>"
             f"<title>{esc(ui('title', lang))} · {esc(result['sample'])}</title>",
             f"<style>{CSS}</style></head><body>",
             f"<h1>{esc(ui('title', lang))} · {esc(result['sample'])}</h1>",
             f"<p class='disclaimer'><strong>{esc(ui('important', lang))}:</strong> {esc(disclaimer(lang))}</p>",
             "<ul>",
             f"<li>{esc(ui('reference', lang))}: {esc(result['reference']['path'])} ({esc(result['reference']['assembly'])})</li>",
             f"<li>{esc(ui('input', lang))}: {esc(result['input']['mode'])} · {esc(', '.join(result['input']['files']))}</li>",
             f"<li>{esc(ui('snv_caller', lang))}: {esc(result['snvCaller'])}</li>",
             f"<li>{esc(ui('created', lang))}: {esc(result['generatedAt'])}</li></ul>"]
    for title, table, empty in sections(result, lang):
        parts.append(f"<h2>{esc(title)}</h2>")
        if not table["rows"]:
            parts.append(f"<p>{esc(empty)}</p>")
            continue
        parts.append("<table><tr>" + "".join(f"<th>{esc(h)}</th>" for h in table["headers"]) + "</tr>")
        for row in table["rows"]:
            css = " class='uncertain'" if uncertain in row else ""
            parts.append(f"<tr{css}>" + "".join(f"<td>{esc(c)}</td>" for c in row) + "</tr>")
        parts.append("</table>")
    parts.append(f"<h2>{esc(ui('versions', lang))}</h2><table><tr><th>Tool</th><th>{esc(ui('version', lang))}</th></tr>")
    parts += [f"<tr><td>{esc(k)}</td><td>{esc(v or ui('unknown', lang))}</td></tr>"
              for k, v in result["toolVersions"].items()]
    parts.append(f"</table><h2>{esc(ui('commands', lang))}</h2>"
                 f"<p>{esc(ui('commands_note', lang))} <code>brca12_commands.sh</code>.</p>"
                 f"<pre>{esc(result['invocation'])}</pre><h2>{esc(ui('limitations', lang))}</h2><ul>")
    parts += [f"<li>{esc(item)}</li>" for item in _limitations(result, lang)]
    parts.append(f"</ul><p class='disclaimer'>{esc(disclaimer(lang))}</p></body></html>")
    return "\n".join(parts)

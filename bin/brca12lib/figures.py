"""SVG figures for the HTML report: coverage per exon and exon dosage per gene.

Self-contained inline SVG (no scripts, no external fonts), readable in print and on screen.
"""
from html import escape

TEXT = {
    "de": {"coverage": "{gene}: mittlere Tiefe je Exon (MAPQ ≥ 20)", "dosage": "{gene}: Exon-Dosis relativ zum Chromosomen-Median",
           "exon": "Exon (MANE)", "depth": "Tiefe", "ratio": "Verhältnis", "uncertain": "unsicher",
           "band": "Bereich ohne Hinweis", "note_cov": "Linien: {low}x (geringe Tiefe) und {min}x (Mindesttiefe je Exon). Orange: als unsicher markiertes Exon.",
           "note_dos": "Graues Band 0,65 bis 1,35: kein Hinweis. Punkte darunter oder darüber sind Kandidaten für Verlust oder Zugewinn, keine CNV-Aufrufe; Bestätigung z. B. mit MLPA."},
    "en": {"coverage": "{gene}: mean depth per exon (MAPQ ≥ 20)", "dosage": "{gene}: exon dosage relative to the chromosome median",
           "exon": "Exon (MANE)", "depth": "Depth", "ratio": "Ratio", "uncertain": "uncertain",
           "band": "no-hint range", "note_cov": "Lines: {low}x (low depth) and {min}x (minimum depth per exon). Orange: exon flagged as uncertain.",
           "note_dos": "Grey band 0.65 to 1.35: no hint. Points below or above are candidates for loss or gain, not CNV calls; confirm e.g. with MLPA."},
}
INK, MUTED, GRID, WARN, BAND = "#1d1d1f", "#636366", "#d2d2d7", "#b25000", "#f0f0f2"
W, H, LEFT, RIGHT, TOP, BOTTOM = 720, 230, 52, 12, 34, 40


def _n(value):
    """10.0 -> 10 for labels."""
    return int(value) if float(value).is_integer() else value


def _exons(result, gene):
    rows = [e for e in result.get("exonCoverage") or [] if e.get("gene") == gene]
    return sorted(rows, key=lambda e: (e.get("exon") or 0))


def _frame(title, ylabel, xlabel):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="100%" role="img" '
            f'aria-label="{escape(title)}" style="max-width:{W}px;font-family:IBM Plex Mono,ui-monospace,Menlo,monospace">',
            f'<text x="{LEFT}" y="18" font-size="13" fill="{INK}">{escape(title)}</text>',
            f'<text x="12" y="{TOP + (H - TOP - BOTTOM) / 2}" font-size="10" fill="{MUTED}" '
            f'transform="rotate(-90 12 {TOP + (H - TOP - BOTTOM) / 2})" text-anchor="middle">{escape(ylabel)}</text>',
            f'<text x="{LEFT + (W - LEFT - RIGHT) / 2}" y="{H - 6}" font-size="10" fill="{MUTED}" text-anchor="middle">{escape(xlabel)}</text>']


def coverage_svg(result, gene, lang="en"):
    exons = _exons(result, gene)
    if not exons:
        return ""
    t = TEXT["de" if lang == "de" else "en"]
    thr = result.get("thresholds") or {}
    low, minimum = _n(thr.get("lowDepth", 10)), _n(thr.get("minExonDepth", 20))
    depths = [e.get("meanDepth") or 0 for e in exons]
    top = max(max(depths), minimum * 1.5, 1) * 1.1
    plot_h, plot_w = H - TOP - BOTTOM, W - LEFT - RIGHT
    y = lambda v: TOP + plot_h - (v / top) * plot_h
    step = plot_w / len(exons)
    out = _frame(t["coverage"].format(gene=gene), t["depth"], t["exon"])
    for tick in range(10, int(top) + 1, 10):
        if tick not in (low, minimum):
            out.append(f'<line x1="{LEFT}" x2="{W - RIGHT}" y1="{y(tick):.1f}" y2="{y(tick):.1f}" stroke="{GRID}" stroke-width="0.6"/>')
            out.append(f'<text x="{LEFT - 6}" y="{y(tick) + 3:.1f}" font-size="9" fill="{MUTED}" text-anchor="end">{tick}x</text>')
    for value, dash in ((low, "3 3"), (minimum, "5 3")):
        out.append(f'<line x1="{LEFT}" x2="{W - RIGHT}" y1="{y(value):.1f}" y2="{y(value):.1f}" stroke="{MUTED}" stroke-dasharray="{dash}"/>')
        out.append(f'<text x="{LEFT - 6}" y="{y(value) + 3:.1f}" font-size="9" fill="{MUTED}" text-anchor="end">{value}x</text>')
    out.append(f'<line x1="{LEFT}" x2="{W - RIGHT}" y1="{y(0):.1f}" y2="{y(0):.1f}" stroke="{GRID}"/>')
    for i, e in enumerate(exons):
        value = e.get("meanDepth") or 0
        uncertain = e.get("status") != "ok"
        x = LEFT + i * step + step * 0.18
        width = step * 0.64
        colour = WARN if uncertain else MUTED
        label = f"{e.get('label')}: {value}x" + (f" ({t['uncertain']})" if uncertain else "")
        out.append(f'<rect x="{x:.1f}" y="{y(value):.1f}" width="{width:.1f}" height="{y(0) - y(value):.1f}" fill="{colour}"><title>{escape(label)}</title></rect>')
        if len(exons) <= 30 or i % 2 == 0:
            out.append(f'<text x="{x + width / 2:.1f}" y="{y(0) + 12:.1f}" font-size="9" fill="{WARN if uncertain else MUTED}" text-anchor="middle">{e.get("exon")}</text>')
    out.append("</svg>")
    return "".join(out) + f'<p class="figure-note">{escape(t["note_cov"].format(low=low, min=minimum))}</p>'


def dosage_svg(result, gene, lang="en", low=0.65, high=1.35):
    exons = [e for e in _exons(result, gene) if e.get("dosageRatio") is not None]
    if not exons:
        return ""
    t = TEXT["de" if lang == "de" else "en"]
    top = max(2.0, max(e["dosageRatio"] for e in exons) * 1.1)
    plot_h, plot_w = H - TOP - BOTTOM, W - LEFT - RIGHT
    y = lambda v: TOP + plot_h - (v / top) * plot_h
    step = plot_w / len(exons)
    out = _frame(t["dosage"].format(gene=gene), t["ratio"], t["exon"])
    out.append(f'<rect x="{LEFT}" y="{y(high):.1f}" width="{plot_w}" height="{y(low) - y(high):.1f}" fill="{BAND}"><title>{escape(t["band"])}</title></rect>')
    for value in (0.5, 1.0, 1.5):
        out.append(f'<line x1="{LEFT}" x2="{W - RIGHT}" y1="{y(value):.1f}" y2="{y(value):.1f}" stroke="{GRID if value != 1.0 else MUTED}"/>')
        out.append(f'<text x="{LEFT - 6}" y="{y(value) + 3:.1f}" font-size="9" fill="{MUTED}" text-anchor="end">{value:.1f}</text>')
    for i, e in enumerate(exons):
        ratio = e["dosageRatio"]
        flagged = ratio < low or ratio > high
        cx = LEFT + i * step + step / 2
        out.append(f'<circle cx="{cx:.1f}" cy="{y(ratio):.1f}" r="{4 if flagged else 3}" fill="{WARN if flagged else INK}">'
                   f'<title>{escape(str(e.get("label")))}: {ratio}</title></circle>')
        if len(exons) <= 30 or i % 2 == 0:
            out.append(f'<text x="{cx:.1f}" y="{H - BOTTOM + 12}" font-size="9" fill="{WARN if flagged else MUTED}" text-anchor="middle">{e.get("exon")}</text>')
    out.append("</svg>")
    return "".join(out) + f'<p class="figure-note">{escape(t["note_dos"])}</p>'


def figures_html(result, lang="en"):
    title = "Abbildungen" if lang == "de" else "Figures"
    blocks = []
    for gene in ("BRCA1", "BRCA2"):
        for svg in (coverage_svg(result, gene, lang), dosage_svg(result, gene, lang)):
            if svg:
                blocks.append(f'<figure class="report-figure">{svg}</figure>')
    return f"<h2>{title}</h2>" + "".join(blocks) if blocks else ""

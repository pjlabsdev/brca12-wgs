"""Coverage je Exon (mosdepth) und unsichere Bereiche."""
import gzip
import re


def _rows(path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            yield line.rstrip("\n").split("\t")


def exon_label(name):
    """BRCA1_ex05 -> ('BRCA1', 5, 'Exon 5 (Legacy 6)'). BRCA1 zählt historisch ohne Exon 4."""
    match = re.match(r"(BRCA[12])\D*(\d+)", name)
    if not match:
        return name.split("_")[0], None, name
    gene, number = match.group(1), int(match.group(2))
    label = f"Exon {number}"
    if gene == "BRCA1" and number >= 4:
        label += f" (Legacy {number + 1})"
    return gene, number, label


def read_chrom_median(path):
    medians = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) == 2:
                medians[fields[0]] = float(fields[1])
    return medians


def exon_table(q20_path, thresholds_path, q0_path, chrom_median, min_depth, min_pct20, max_mapq0):
    q0 = {(r[0], int(r[1]), int(r[2])): float(r[-1]) for r in _rows(q0_path)}
    thresholds = {}
    for r in _rows(thresholds_path):
        length = int(r[2]) - int(r[1])
        counts = [int(v) for v in r[4:7]]
        thresholds[(r[0], int(r[1]), int(r[2]))] = [100.0 * c / length if length else 0.0 for c in counts]
    exons = []
    for r in _rows(q20_path):
        chrom, start, end, name, mean = r[0], int(r[1]), int(r[2]), r[3], float(r[-1])
        gene, number, label = exon_label(name)
        pct10, pct20, pct30 = thresholds.get((chrom, start, end), [None, None, None])
        all_mean = q0.get((chrom, start, end))
        mapq0 = (1.0 - mean / all_mean) if all_mean else None
        median = chrom_median.get(chrom)
        dosage = (mean / median) if median else None
        reasons = []
        if mean < min_depth:
            reasons.append(f"mittlere Tiefe {mean:.1f}x < {min_depth:g}x")
        if pct20 is not None and pct20 < min_pct20:
            reasons.append(f"nur {pct20:.1f} % der Basen >= 20x")
        if mapq0 is not None and mapq0 > max_mapq0:
            reasons.append(f"{100 * mapq0:.0f} % der Reads mit MAPQ < 20 (Pseudogen/Repeat)")
        exons.append({
            "gene": gene, "exon": number, "label": label, "chrom": chrom,
            "start": start + 1, "end": end, "length": end - start,
            "meanDepth": round(mean, 1),
            "pct10x": _round(pct10), "pct20x": _round(pct20), "pct30x": _round(pct30),
            "lowMapqFraction": _round(mapq0, 3),
            "dosageRatio": _round(dosage, 2),
            "status": "unsicher" if reasons else "ok",
            "reasons": reasons,
        })
    exons.sort(key=lambda e: (e["gene"], e["exon"] or 0))
    return exons


def _round(value, digits=1):
    return round(value, digits) if value is not None else None


def low_coverage_intervals(perbase_path, regions, threshold, min_length=1):
    """Zusammenhängende Intervalle unter `threshold` innerhalb der Genbereiche (1-basiert)."""
    spans = {g: (c, s, e) for c, s, e, g in regions}
    intervals = []
    current = None
    for r in _rows(perbase_path):
        chrom, start, end, depth = r[0], int(r[1]), int(r[2]), float(r[3])
        gene = next((g for g, (c, s, e) in spans.items() if c == chrom and start < e and end > s), None)
        if gene is None or depth >= threshold:
            if current:
                intervals.append(current)
                current = None
            continue
        c, s, e = spans[gene]
        start, end = max(start, s), min(end, e)
        if current and current["gene"] == gene and current["end"] == start:
            current["end"] = end
            current["minDepth"] = min(current["minDepth"], depth)
        else:
            if current:
                intervals.append(current)
            current = {"gene": gene, "chrom": chrom, "start": start, "end": end, "minDepth": depth}
    if current:
        intervals.append(current)
    result = []
    for i in intervals:
        if i["end"] - i["start"] >= min_length:
            result.append({**i, "start": i["start"] + 1, "length": i["end"] - i["start"]})
    return result


def gene_summary(exons, gene):
    rows = [e for e in exons if e["gene"] == gene]
    if not rows:
        return {}
    total = sum(e["length"] for e in rows)
    mean = sum(e["meanDepth"] * e["length"] for e in rows) / total
    pct20 = sum((e["pct20x"] or 0) * e["length"] for e in rows) / total
    return {"exonCount": len(rows), "codingMeanDepth": round(mean, 1),
            "pct20x": round(pct20, 1),
            "uncertainExons": [e["label"] for e in rows if e["status"] != "ok"]}

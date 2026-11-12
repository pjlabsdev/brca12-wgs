"""CNV-Kandidaten (CNVkit) und strukturelle Varianten (Manta) in BRCA1/BRCA2."""
import re

from .vcf import VCF, to_int


def _overlaps(chrom, start, end, region):
    c, s, e, _ = region
    return chrom == c and start <= e and end > s


def _gene_for(chrom, start, end, regions):
    return [r[3] for r in regions if _overlaps(chrom, start, end, r)]


def _affected_exons(chrom, start, end, exons, gene):
    return [e["label"] for e in exons
            if e["gene"] == gene and e["chrom"] == chrom and e["start"] <= end and e["end"] >= start]


def read_table(path):
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            if line.strip():
                yield dict(zip(header, line.rstrip("\n").split("\t")))


def cnv_segments(cns_path, cnr_path, regions, exons):
    """CNVkit-Segmente, die BRCA1/BRCA2 überlappen. cn != 2 gilt als Kandidat."""
    bins = list(read_table(cnr_path))
    results = []
    for seg in read_table(cns_path):
        chrom, start, end = seg["chromosome"], int(seg["start"]), int(seg["end"])
        for gene in _gene_for(chrom, start, end, regions):
            cn = to_int(seg.get("cn"))
            log2 = float(seg["log2"])
            region = next(r for r in regions if r[3] == gene)
            gene_bins = [float(b["log2"]) for b in bins
                         if b["chromosome"] == chrom and _overlaps(chrom, int(b["start"]), int(b["end"]), region)]
            kind = "Deletion" if cn is not None and cn < 2 else "Duplikation" if cn is not None and cn > 2 else "kein Hinweis"
            overlap_start, overlap_end = max(start, region[1]), min(end, region[2])
            results.append({
                "gene": gene, "chrom": chrom, "start": start + 1, "end": end,
                "overlapStart": overlap_start + 1, "overlapEnd": overlap_end,
                "log2": round(log2, 3), "copyNumber": cn, "type": kind,
                "probes": to_int(seg.get("probes")),
                "geneBinMeanLog2": round(sum(gene_bins) / len(gene_bins), 3) if gene_bins else None,
                "affectedExons": _affected_exons(chrom, overlap_start + 1, overlap_end, exons, gene),
                "caller": "CNVkit (WGS, flache Referenz)",
                "flags": ["Kandidat, Bestätigung per MLPA oder validierter CNV-Analyse nötig"] if kind != "kein Hinweis" else [],
            })
    return results


def exon_dosage_hints(exons, low=0.65, high=1.35):
    hints = []
    for exon in exons:
        ratio = exon.get("dosageRatio")
        if ratio is None or exon["status"] != "ok":
            continue
        if ratio < low or ratio > high:
            hints.append({
                "gene": exon["gene"], "exon": exon["label"], "ratio": ratio,
                "type": "möglicher Verlust" if ratio < low else "möglicher Zugewinn",
                "note": "Exonweise Tiefe relativ zum Chromosomen-Median; nur ein Hinweis, kein CNV-Call.",
            })
    return hints


_BND_MATE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def structural_variants(path, regions, exons):
    vcf = VCF(path)
    results = []
    for record in vcf.records:
        info = record["info"]
        svtype = str(info.get("SVTYPE", "SV"))
        chrom, pos = record["chrom"], record["pos"]
        end = to_int(info.get("END")) or pos
        mate_chrom = mate_pos = None
        alt = record["alts"][0]
        if svtype == "BND":
            match = _BND_MATE.search(alt)
            if match:
                mate_chrom, mate_pos = match.group(1), int(match.group(2))
        genes = set(_gene_for(chrom, pos - 1, end, regions))
        if mate_chrom:
            genes |= set(_gene_for(mate_chrom, mate_pos - 1, mate_pos, regions))
        if not genes:
            continue
        sample = record["sample"]
        pr = sample.get("PR", "")
        sr = sample.get("SR", "")
        for gene in sorted(genes):
            flags = []
            if record["filter"] not in ("PASS", "."):
                flags.append(f"FILTER={record['filter']}")
            if "IMPRECISE" in info:
                flags.append("unpräzise Bruchpunkte")
            flags.append("Bestätigung mit orthogonaler Methode nötig")
            results.append({
                "gene": gene, "id": record["id"], "type": svtype, "chrom": chrom, "pos": pos,
                "end": end if svtype != "BND" else None,
                "mateChrom": mate_chrom, "matePos": mate_pos,
                "length": to_int(str(info.get("SVLEN", "")).lstrip("-")),
                "genotype": sample.get("GT"), "pairedReads": pr or None, "splitReads": sr or None,
                "filter": record["filter"],
                "affectedExons": _affected_exons(chrom, pos, end, exons, gene),
                "caller": "Manta", "flags": flags,
            })
    return results

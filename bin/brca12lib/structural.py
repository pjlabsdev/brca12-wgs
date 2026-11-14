"""CNV-Kandidaten (CNVkit) und strukturelle Varianten (Manta) in BRCA1/BRCA2."""
import json
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
    """Exon depth relative to the chromosome median. A heterozygous deletion halves the depth and
    thereby also makes the exon "unsicher" (low depth), so such exons must not be skipped."""
    hints = []
    for exon in exons:
        ratio = exon.get("dosageRatio")
        if ratio is None:
            continue
        if ratio < low or ratio > high:
            note = "Exonweise Tiefe relativ zum Chromosomen-Median; nur ein Hinweis, kein CNV-Call."
            if exon["status"] != "ok":
                note += " Exon ist zugleich als unsicher markiert: Verlust oder Abdeckungsproblem, nicht unterscheidbar."
            hints.append({
                "gene": exon["gene"], "exon": exon["label"], "ratio": ratio,
                "type": "möglicher Verlust" if ratio < low else "möglicher Zugewinn",
                "note": note,
            })
    return hints


# start of the Alu consensus (left monomer) and its reverse complement
_ALU_MOTIFS = ("GGCCGGGCGCGGTGGCTCA", "TGAGCCACCGCGCCCGGCC")
_BND_MATE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def structural_variants(path, regions, exons):
    vcf = VCF(path)
    results = []
    ids = {record["id"] for record in vcf.records}
    reported_bnd = set()
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
        mate_id = str(info.get("MATEID", "")) or None
        if svtype == "BND":
            # Manta writes one BND record per junction end, linked by MATEID. Report each
            # junction once, from its in-gene end; the outside end only stands in when its
            # in-gene mate record is missing from the VCF.
            if mate_id and mate_id in reported_bnd:
                continue
            if not genes and mate_chrom:
                mate_genes = set(_gene_for(mate_chrom, mate_pos - 1, mate_pos, regions))
                if not mate_genes or (mate_id and mate_id in ids):
                    continue
                genes = mate_genes
            reported_bnd.add(record["id"])
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
            inserted = " ".join([alt, str(info.get("LEFT_SVINSSEQ", "")), str(info.get("RIGHT_SVINSSEQ", ""))]).upper()
            if svtype == "INS" and any(motif in inserted for motif in _ALU_MOTIFS):
                flags.append("eingefügte Sequenz mit Alu-Signatur: mögliche Alu-Insertion")
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
    # Two breakends a few bases apart whose mates lie elsewhere are the signature of an insertion
    # of foreign sequence (e.g. a mobile element such as Alu); Manta cannot name it as INS.
    bnds = [r for r in results if r["type"] == "BND"]
    for a in bnds:
        if any(b is not a and b["gene"] == a["gene"] and b["chrom"] == a["chrom"] and abs(b["pos"] - a["pos"]) <= 50
               for b in bnds):
            a["flags"].insert(0, "Bruchpunktpaar an derselben Stelle: mögliche Insertion (z. B. mobiles Element wie Alu)")
    return results


def insertion_signals(path, exons):
    """Insertionssignaturen aus abgeschnittenen Leseenden (brca12_insertion_signals.py)."""
    if not path:
        return []
    with open(path) as handle:
        rows = json.load(handle)
    results = []
    for row in rows:
        flags = []
        kind = "Insertionssignatur"
        if "tandem copy" in row.get("signature", []):
            kind = "Tandemkopie-Signatur"
            if row.get("clippedReadsLeftJunction") or row.get("clippedReadsRightJunction"):
                flags.append("abgeschnittene Leseenden an beiden Grenzen einer Wiederholungseinheit: mögliche zusätzliche Tandemkopie")
            if "known allele" in row.get("signature", []):
                flags.append(row.get("label") or "bekannte zusätzliche Tandemkopie")
                flags.append("Markerbase der zusätzlichen Kopie in den Lesungen vorhanden")
        elif "Alu" in row.get("signature", []):
            flags.append("abgeschnittene Leseenden mit Alu-Sequenz: mögliche Alu-Insertion")
        else:
            flags.append("abgeschnittene Leseenden von beiden Seiten: mögliche Insertion")
        flags.append("Hinweis aus Leseenden, kein Variantenaufruf; Bestätigung mit orthogonaler Methode nötig")
        results.append(dict(row, type=kind,
                            affectedExons=_affected_exons(row["chrom"], row["pos"], row["end"], exons, row["gene"]),
                            caller="brca12_insertion_signals (Leseenden)", flags=flags))
    return results

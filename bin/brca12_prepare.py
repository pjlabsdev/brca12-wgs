#!/usr/bin/env python3
"""Eingabeprüfung und BED-Aufbereitung für die BRCA1/BRCA2-Pipeline (nur GRCh38).

Prüft anhand der .fai, dass die Referenz GRCh38 ist, gleicht die Contig-Schreibweise
(chr17 bzw. 17) an und erzeugt die Regionsdateien für Calling, SV-Suche, Extraktion und
Exon-Coverage. Zusätzlich zwei kleine Helfer für CNV-Ziele und den Chromosomen-Median.
"""
import argparse
import gzip
import json
import statistics
import sys

GRCH38_LENGTHS = {"13": 114364328, "17": 83257441}
GRCH37_LENGTHS = {"13": 115169878, "17": 81195210}

# Verbindliche Genbereiche (GRCh38, 1-basiert, inklusive) laut Auftrag.
GENE_LOCI = {
    "BRCA1": {"chrom": "17", "start": 43044295, "end": 43170327, "strand": "-"},
    "BRCA2": {"chrom": "13", "start": 32315077, "end": 32400268, "strand": "+"},
}


class InputError(Exception):
    pass


def read_fai(path):
    lengths = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def bare(contig):
    return contig[3:] if contig.lower().startswith("chr") else contig


def check_reference(lengths):
    """Gibt (Contig-Präfix, Bericht) zurück oder wirft InputError für nicht GRCh38."""
    report = {"assembly": None, "ok": False, "contigs": {}, "messages": []}
    style = None
    for number, expected in GRCH38_LENGTHS.items():
        candidates = [c for c in (f"chr{number}", number) if c in lengths]
        if not candidates:
            raise InputError(f"Contig chr{number}/{number} missing from the reference .fai.")
        name = candidates[0]
        style = "chr" if name.startswith("chr") else ""
        observed = lengths[name]
        report["contigs"][name] = observed
        if observed == GRCH37_LENGTHS[number]:
            raise InputError(f"{name} has the GRCh37 length {observed}. This pipeline only accepts GRCh38/hg38.")
        if observed != expected:
            raise InputError(f"{name} has length {observed}, expected {expected} for GRCh38.")
    report.update({"assembly": "GRCh38", "ok": True, "contigStyle": style or "numeric"})
    report["messages"].append("Reference recognised as GRCh38 from the contig lengths of chr13 and chr17.")
    return style, report


def read_bed(path):
    rows = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise InputError(f"BED line with fewer than 3 columns in {path}: {line.strip()}")
            rows.append(fields)
    return rows


def gene_of(name):
    upper = name.upper()
    for gene in GENE_LOCI:
        if upper.startswith(gene):
            return gene
    return None


def normalize_regions(rows, style, lengths):
    """Ordnet jeder Zeile BRCA1/BRCA2 zu und prüft, dass sie im verbindlichen Locus liegt."""
    out = []
    for fields in rows:
        chrom = style + bare(fields[0])
        start, end = int(fields[1]), int(fields[2])
        name = fields[3] if len(fields) > 3 else ""
        gene = gene_of(name)
        if gene is None:
            gene = next((g for g, l in GENE_LOCI.items()
                         if bare(chrom) == l["chrom"] and start < l["end"] and end >= l["start"]), None)
        if gene is None:
            raise InputError(f"Region {fields[0]}:{start}-{end} lies in neither BRCA1 nor BRCA2.")
        locus = GENE_LOCI[gene]
        if bare(chrom) != locus["chrom"] or end < locus["start"] - 5000 or start > locus["end"] + 5000:
            raise InputError(f"Region {name or gene} {chrom}:{start}-{end} lies outside the GRCh38 locus of {gene}.")
        if chrom not in lengths:
            raise InputError(f"Contig {chrom} missing from the reference.")
        out.append({"chrom": chrom, "start": start, "end": end, "gene": gene,
                    "name": name or gene, "strand": fields[5] if len(fields) > 5 else locus["strand"]})
    genes = {r["gene"] for r in out}
    missing = set(GENE_LOCI) - genes
    if missing:
        raise InputError(f"BED has no region for {', '.join(sorted(missing))}.")
    return out


def exact_gene_regions(regions, style):
    """Eine Zeile je Gen: Vereinigung aus BED und verbindlichem Locus (0-basiert)."""
    merged = {}
    for gene, locus in GENE_LOCI.items():
        merged[gene] = [locus["start"] - 1, locus["end"]]
    for r in regions:
        span = merged[r["gene"]]
        span[0] = min(span[0], r["start"])
        span[1] = max(span[1], r["end"])
    return [(style + GENE_LOCI[g]["chrom"], s, e, g) for g, (s, e) in merged.items()]


def padded(rows, padding, lengths):
    return [(c, max(0, s - padding), min(lengths[c], e + padding), g) for c, s, e, g in rows]


def write_bed(path, rows):
    with open(path, "w") as handle:
        for row in sorted(rows, key=lambda r: (r[0], r[1])):
            handle.write("\t".join(str(v) for v in row) + "\n")


def gene_targets(path, style, windows, out_path, exons=None, keep_columns=False):
    """Keep the target (or primer) intervals that overlap BRCA1/BRCA2 and harmonise contig names.

    For targets the exons that no target covers are listed: their coverage cannot be judged.
    """
    rows = []
    for fields in read_bed(path):
        chrom = style + bare(fields[0])
        start, end = int(fields[1]), int(fields[2])
        if end <= start:
            raise InputError(f"BED line with end <= start in {path}: {' '.join(fields)}")
        if any(c == chrom and start < e and end > s for c, s, e, _ in windows):
            name = fields[3] if len(fields) > 3 and fields[3] else f"{chrom}:{start + 1}-{end}"
            extra = fields[4:6] if keep_columns and len(fields) >= 6 else []
            rows.append((chrom, start, end, name, *extra))
    if not rows:
        raise InputError(f"{path} has no interval inside BRCA1/BRCA2 (GRCh38).")
    write_bed(out_path, rows)
    info = {"file": path, "intervalsInGenes": len(rows)}
    if exons is not None:
        info["exonsWithoutTarget"] = [e["name"] for e in exons
                                      if not any(r[0] == e["chrom"] and r[1] < e["end"] and r[2] > e["start"]
                                                 for r in rows)]
    return info


def prepare(args):
    lengths = read_fai(args.fai)
    style, report = check_reference(lengths)
    regions = normalize_regions(read_bed(args.regions), style, lengths)
    exons = normalize_regions(read_bed(args.exons), style, lengths)
    exact = exact_gene_regions(regions, style)
    write_bed(args.out_exact, exact)
    write_bed(args.out_calling, padded(exact, args.calling_padding, lengths))
    write_bed(args.out_sv, padded(exact, args.sv_padding, lengths))
    write_bed(args.out_exons, [(e["chrom"], e["start"], e["end"], e["name"], 0, e["strand"]) for e in exons])
    report["regions"] = [{"chrom": c, "start": s + 1, "end": e, "gene": g} for c, s, e, g in exact]
    if args.targets and args.out_targets:
        report["targets"] = gene_targets(args.targets, style, padded(exact, args.calling_padding, lengths),
                                         args.out_targets, exons)
    if args.primers and args.out_primers:
        report["primers"] = gene_targets(args.primers, style, padded(exact, args.calling_padding, lengths),
                                         args.out_primers, None, keep_columns=True)
    report["exonCount"] = {g: sum(1 for e in exons if e["gene"] == g) for g in GENE_LOCI}
    report["chromLengths"] = {style + n: lengths[style + n] for n in GRCH38_LENGTHS}
    with open(args.out_check, "w") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)


def cnv_targets(check_path, out_path, access_path=None):
    """CNVkit WGS mode needs one BED for --targets and --access: the accessible part of chr13/chr17."""
    with open(check_path) as handle:
        check = json.load(handle)
    lengths = check["chromLengths"]
    if not access_path:
        write_bed(out_path, [(c, 0, length, "chrom") for c, length in lengths.items()])
        return
    rows = []
    with open(access_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] in lengths:
                rows.append((fields[0], int(fields[1]), int(fields[2]), "access"))
    if not rows:
        raise InputError("CNVkit access BED has no intervals on chr13/chr17")
    write_bed(out_path, rows)


def chrom_median(paths, out_path):
    with open(out_path, "w") as out:
        for path in paths:
            depths = {}
            with gzip.open(path, "rt") as handle:
                for line in handle:
                    fields = line.split("\t")
                    value = float(fields[-1])
                    if value > 0:
                        depths.setdefault(fields[0], []).append(value)
            for chrom, values in depths.items():
                out.write(f"{chrom}\t{statistics.median(values):.3f}\n")


def exon_baseline(path, out_path):
    """Targeted assays: the median of the per-exon mean depths (both genes) for every chromosome."""
    means, chroms = [], set()
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            chroms.add(fields[0])
            value = float(fields[-1])
            if value > 0:
                means.append(value)
    baseline = statistics.median(means) if means else 0.0
    with open(out_path, "w") as out:
        for chrom in sorted(chroms):
            out.write(f"{chrom}\t{baseline:.3f}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fai")
    parser.add_argument("--regions")
    parser.add_argument("--exons")
    parser.add_argument("--calling-padding", type=int, default=200)
    parser.add_argument("--sv-padding", type=int, default=1_000_000)
    parser.add_argument("--out-check")
    parser.add_argument("--out-calling")
    parser.add_argument("--out-sv")
    parser.add_argument("--out-exact")
    parser.add_argument("--out-exons")
    parser.add_argument("--targets", help="target BED of an exome, panel or amplicon assay")
    parser.add_argument("--out-targets")
    parser.add_argument("--primers", help="primer BED of an amplicon assay (samtools ampliconclip)")
    parser.add_argument("--out-primers")
    parser.add_argument("--cnv-targets")
    parser.add_argument("--out-cnv-targets")
    parser.add_argument("--cnv-access")
    parser.add_argument("--chrom-median", nargs="+")
    parser.add_argument("--out-chrom-median")
    parser.add_argument("--exon-baseline", help="mosdepth per-exon means (targeted assays)")
    args = parser.parse_args(argv)
    try:
        if args.cnv_targets:
            cnv_targets(args.cnv_targets, args.out_cnv_targets, args.cnv_access)
        elif args.exon_baseline:
            exon_baseline(args.exon_baseline, args.out_chrom_median)
        elif args.chrom_median:
            chrom_median(args.chrom_median, args.out_chrom_median)
        else:
            prepare(args)
    except InputError as error:
        print(f"[brca12] input error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

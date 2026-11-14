#!/usr/bin/env python3
"""Read and target quality control for brca12-wgs (all assays).

Summarises samtools stats over the BRCA1/BRCA2 calling regions (reads, mapping, duplicates,
insert size, error rate), and for exome, panel and amplicon data the coverage of every target:
mean and minimum depth, share of bases >= 20x, weakly covered targets, uniformity (share of
target bases with at least 20 % of the mean depth), fold-80 base penalty and on-target rate.
For amplicon data with a primer BED it also reports how many reads had a primer clipped
(samtools ampliconclip). Writes read_qc.json.
"""
import argparse
import gzip
import json
import statistics
import sys


def read_samtools_stats(path):
    sn, insert = {}, []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "SN":
                key = fields[1].rstrip(":")
                try:
                    sn[key] = float(fields[2])
                except ValueError:
                    pass
            elif fields[0] == "IS":
                insert.append((int(fields[1]), int(fields[2])))
    total = sn.get("raw total sequences", 0)
    mapped = sn.get("reads mapped", 0)
    pairs = sum(n for _, n in insert)
    median_insert = None
    if pairs:
        seen = 0
        for size, n in insert:
            seen += n
            if seen >= pairs / 2:
                median_insert = size
                break

    def pct(a, b):
        return round(100.0 * a / b, 2) if b else None

    return {
        "readsInRegions": int(total),
        "mappedPct": pct(mapped, total),
        "properlyPairedPct": pct(sn.get("reads properly paired", 0), mapped),
        "duplicatePct": pct(sn.get("reads duplicated", 0), mapped),
        "mapq0Pct": pct(sn.get("reads MQ0", 0), mapped),
        "medianInsertSize": median_insert,
        "meanInsertSize": sn.get("insert size average"),
        "errorRate": round(sn["error rate"], 5) if "error rate" in sn else None,
        "averageReadLength": sn.get("average length"),
        "averageBaseQuality": sn.get("average quality"),
    }


def read_idxstats(path):
    mapped = unmapped = 0
    with open(path) as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) >= 4:
                mapped += int(fields[2])
                unmapped += int(fields[3])
    return mapped, unmapped


def read_targets(path):
    targets = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            name = f[3] if len(f) > 3 and f[3] else f"{f[0]}:{int(f[1]) + 1}-{f[2]}"
            targets.append({"name": name, "chrom": f[0], "start": int(f[1]), "end": int(f[2])})
    return targets


def target_coverage(targets, depth_path, min_target_depth):
    """depth_path: samtools depth -a output (chrom, 1-based pos, depth) over the targets, gzip."""
    index = {}
    for i, t in enumerate(targets):
        index.setdefault(t["chrom"], []).append((t["start"], t["end"], i))
    depths = [dict() for _ in targets]
    opener = gzip.open if depth_path.endswith(".gz") else open
    with opener(depth_path, "rt") as handle:
        for line in handle:
            chrom, pos, depth = line.split("\t")[:3]
            pos, depth = int(pos), int(depth)
            for start, end, i in index.get(chrom, ()):
                if start < pos <= end:
                    depths[i][pos] = depth
    all_bases = []
    rows = []
    for t, d in zip(targets, depths):
        length = t["end"] - t["start"]
        values = [d.get(p, 0) for p in range(t["start"] + 1, t["end"] + 1)]
        all_bases.extend(values)
        mean = sum(values) / length if length else 0.0
        rows.append({"name": t["name"], "chrom": t["chrom"], "start": t["start"] + 1, "end": t["end"],
                     "length": length, "meanDepth": round(mean, 1), "minDepth": min(values) if values else 0,
                     "pct20x": round(100.0 * sum(v >= 20 for v in values) / length, 1) if length else None,
                     "status": "schwach abgedeckt" if mean < min_target_depth else "ok"})
    summary = {}
    if all_bases:
        mean = sum(all_bases) / len(all_bases)
        ordered = sorted(all_bases)
        p20 = ordered[int(0.2 * (len(ordered) - 1))]
        summary = {
            "targetCount": len(rows),
            "targetBases": len(all_bases),
            "meanDepth": round(mean, 1),
            "medianTargetDepth": round(statistics.median(r["meanDepth"] for r in rows), 1),
            "pct20x": round(100.0 * sum(v >= 20 for v in all_bases) / len(all_bases), 1),
            "uniformityPct": round(100.0 * sum(v >= 0.2 * mean for v in all_bases) / len(all_bases), 1)
            if mean else None,
            "fold80BasePenalty": round(mean / p20, 2) if p20 else None,
            "weakTargetCount": sum(r["status"] != "ok" for r in rows),
            "minTargetDepth": min_target_depth,
        }
    return rows, summary


def read_clip_stats(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            if ":" in line and not line.startswith("COMMAND"):
                key, value = line.split(":", 1)
                try:
                    values[key.strip()] = int(value.strip())
                except ValueError:
                    pass
    total = values.get("TOTAL READS", 0) - values.get("EXCLUDED", 0)
    clipped = values.get("TOTAL CLIPPED", 0)
    return {"readsConsidered": total, "readsWithPrimerClipped": clipped,
            "clippedPct": round(100.0 * clipped / total, 1) if total else None,
            "readsFilteredTooShort": values.get("FILTERED", 0)}


def build(args):
    result = {"assay": args.assay, "regionStats": read_samtools_stats(args.stats)}
    flags = []
    rs = result["regionStats"]
    if rs["duplicatePct"] is not None and rs["duplicatePct"] > 30 and args.assay != "amplicon":
        flags.append(f"hoher Duplikatanteil ({rs['duplicatePct']:.1f} %)")
    if rs["mappedPct"] is not None and rs["mappedPct"] < 95:
        flags.append(f"geringer Anteil kartierter Reads ({rs['mappedPct']:.1f} %)")
    if rs["errorRate"] is not None and rs["errorRate"] > 0.01:
        flags.append(f"hohe Fehlerrate ({100 * rs['errorRate']:.2f} %)")
    if args.idxstats:
        mapped, unmapped = read_idxstats(args.idxstats)
        result["totalMappedReads"] = mapped
        result["totalUnmappedReads"] = unmapped
    if args.targets and args.target_depth:
        rows, summary = target_coverage(read_targets(args.targets), args.target_depth, args.min_target_depth)
        if args.on_target_reads and args.idxstats and result.get("totalMappedReads"):
            with open(args.on_target_reads) as handle:
                on_target = int(handle.read().split()[0])
            summary["onTargetReadPct"] = round(100.0 * on_target / result["totalMappedReads"], 1)
        result["targets"] = rows
        result["targetSummary"] = summary
        if summary.get("weakTargetCount"):
            flags.append(f"{summary['weakTargetCount']} Ziel(e) mit mittlerer Tiefe < {args.min_target_depth:g}x")
        if summary.get("uniformityPct") is not None and summary["uniformityPct"] < 90:
            flags.append(f"ungleichmäßige Abdeckung (nur {summary['uniformityPct']:.1f} % der Zielbasen >= 0,2 x Mittel)")
        if args.out_baseline and rows:
            # Per-exon dosage baseline for targeted data: median depth over all targets (both genes),
            # written for every chromosome so that a whole-gene loss remains visible.
            baseline = statistics.median(r["meanDepth"] for r in rows if r["meanDepth"] > 0) if any(
                r["meanDepth"] > 0 for r in rows) else 0
            with open(args.out_baseline, "w") as out:
                for chrom in sorted({r["chrom"] for r in rows}):
                    out.write(f"{chrom}\t{baseline:.3f}\n")
    if args.clip_stats:
        result["primerClipping"] = read_clip_stats(args.clip_stats)
        pct = result["primerClipping"]["clippedPct"]
        if pct is not None and pct < 80:
            flags.append(f"nur {pct:.1f} % der Reads mit abgeschnittenem Primer")
    result["flags"] = flags
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--assay", default="wgs", choices=("wgs", "exome", "panel", "amplicon"))
    p.add_argument("--stats", required=True, help="samtools stats output over the calling regions")
    p.add_argument("--idxstats", default=None)
    p.add_argument("--targets", default=None, help="target BED (exome/panel/amplicon)")
    p.add_argument("--target-depth", default=None, help="samtools depth -a over the targets (gzip)")
    p.add_argument("--on-target-reads", default=None, help="file with the number of reads overlapping the targets")
    p.add_argument("--clip-stats", default=None, help="samtools ampliconclip statistics")
    p.add_argument("--min-target-depth", type=float, default=30)
    p.add_argument("--out-baseline", default=None, help="write the dosage baseline for targeted data here")
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)
    with open(args.out, "w") as handle:
        json.dump(build(args), handle, indent=1, ensure_ascii=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())

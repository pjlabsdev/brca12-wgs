#!/usr/bin/env python3
"""Pipeline calls against assembly-based truth (HPRC v1.1 pangenome VCF) inside BRCA1/BRCA2.

Small variants (length difference below 50 bp) are compared with RTG vcfeval per sample, as in
giab_benchmark.py. There is no confident-region file for this truth; positions where the sample has
a missing allele in the pangenome VCF and the surroundings of its larger variants are excluded.
Variants of 50 bp or more are grouped into sites and looked up in the Manta calls, the
small-variant calls and the read-end insertion signatures.
"""
import argparse
import gzip
import json
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "validation"))
import giab_benchmark  # noqa: E402

GENES = (("chr17", 43044294, 43170327, "BRCA1"), ("chr13", 32315076, 32400268, "BRCA2"))
SV_MIN = 50
SLOP = 150


def records(path, sample):
    """(chrom, pos, ref, alt, gt) of one sample for every record of the region VCF."""
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                col = f.index(sample)
                continue
            yield f[0], int(f[1]), f[3], f[4], f[col].split(":")[0]


def sites(truth_dir):
    """Positions with a variant of SV_MIN or more, merged into sites when closer than SLOP."""
    found = []
    for chrom in ("chr13", "chr17"):
        with gzip.open(os.path.join(truth_dir, f"hprc.{chrom}.vcf.gz"), "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 6)
                if abs(len(f[4]) - len(f[3])) >= SV_MIN:
                    found.append((f[0], int(f[1]), int(f[1]) + len(f[3]) - 1, len(f[4]) - len(f[3])))
    merged = []
    for chrom, start, end, delta in sorted(found):
        if merged and merged[-1]["chrom"] == chrom and start - merged[-1]["end"] <= SLOP:
            merged[-1]["end"] = max(merged[-1]["end"], end)
            merged[-1]["lengths"].append(delta)
        else:
            merged.append({"chrom": chrom, "start": start, "end": end, "lengths": [delta]})
    for m in merged:
        m["id"] = f"{m['chrom']}:{m['start']}"
        m["gene"] = next(g for c, s, e, g in GENES if c == m["chrom"] and s < m["start"] <= e)
        m["lengths"] = sorted(set(m["lengths"]))
    return merged


def prepare_truth(truth_dir, sample, outdir):
    """Per-sample truth VCF with small variants, and a BED of the gene regions minus unusable parts."""
    excluded, lines, header = [], [], []
    for chrom in ("chr13", "chr17"):
        path = os.path.join(truth_dir, f"hprc.{chrom}.vcf.gz")
        if not header:
            with gzip.open(path, "rt") as handle:
                for line in handle:
                    if line.startswith("##"):
                        header.append(line)
                    else:
                        break
        for c, pos, ref, alt, gt in records(path, sample):
            alleles = gt.replace("/", "|").split("|")
            if "." in alleles or len(alleles) != 2:
                excluded.append((c, pos - 11, pos + len(ref) + 10))
                continue
            if alleles == ["0", "0"]:
                continue
            if abs(len(alt) - len(ref)) >= SV_MIN:
                excluded.append((c, pos - SLOP, pos + len(ref) + SLOP))
                continue
            lines.append("\t".join([c, str(pos), ".", ref, alt, ".", "PASS", ".", "GT", "|".join(alleles)]) + "\n")
    vcf = os.path.join(outdir, f"{sample}.truth.vcf")
    with open(vcf, "w") as handle:
        handle.writelines(h for h in header if h.startswith(("##fileformat", "##contig")))
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")
        handle.writelines(lines)
    subprocess.run(["bgzip", "-f", vcf], check=True)
    subprocess.run(["tabix", "-f", vcf + ".gz"], check=True)
    keep = []
    for chrom, start, end, _ in GENES:
        cuts = sorted((max(s, start), min(e, end)) for c, s, e in excluded if c == chrom and e > start and s < end)
        at = start
        for s, e in cuts:
            if s > at:
                keep.append((chrom, at, s))
            at = max(at, e)
        if at < end:
            keep.append((chrom, at, end))
    bed = os.path.join(outdir, f"{sample}.truth.bed")
    with open(bed, "w") as handle:
        handle.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in keep)
    return vcf + ".gz", bed, len(lines)


def site_genotype(truth_dir, sample, site):
    """carrier / reference / unknown and the allele lengths of the sample at a site."""
    path = os.path.join(truth_dir, f"hprc.{site['chrom']}.vcf.gz")
    alt, missing, lengths = 0, False, []
    for c, pos, ref, a, gt in records(path, sample):
        if not (site["start"] <= pos <= site["end"]) or abs(len(a) - len(ref)) < SV_MIN:
            continue
        alleles = gt.replace("/", "|").split("|")
        if "." in alleles:
            missing = True
        n = sum(x not in ("0", ".") for x in alleles)
        if n:
            alt += n
            lengths.append(len(a) - len(ref))
    if alt:
        return ("homozygous" if alt >= 2 and len(set(lengths)) == 1 else "carrier"), lengths
    return ("unknown" if missing else "reference"), []


def evidence(site, report, calls_vcf):
    lo, hi = site["start"] - SLOP, site["end"] + SLOP
    found = {"manta": [], "smallVariantCaller": [], "readEndSignal": []}
    for sv in report.get("svs", []):
        if sv["chrom"] == site["chrom"] and lo <= sv["pos"] <= hi:
            found["manta"].append({k: sv[k] for k in ("type", "pos", "end", "length", "genotype", "filter")})
    out = subprocess.run(["bcftools", "query", "-r", f"{site['chrom']}:{lo}-{hi}", "-f", "%POS\t%REF\t%ALT[\t%GT]\n", calls_vcf],
                         capture_output=True, text=True, check=True)
    for line in out.stdout.splitlines():
        pos, ref, alts, gt = line.split("\t")
        for alt in alts.split(","):
            if abs(len(alt) - len(ref)) >= 30:
                found["smallVariantCaller"].append({"pos": int(pos), "lengthDifference": len(alt) - len(ref), "genotype": gt})
    for x in report.get("insertionSignals", []):
        if x["chrom"] == site["chrom"] and lo <= x["pos"] <= hi:
            found["readEndSignal"].append({k: x[k] for k in ("pos", "end", "targetSiteDuplication", "clippedReadsLeftJunction",
                                                              "clippedReadsRightJunction", "clippedFraction", "signature")})
    return found


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", required=True)
    p.add_argument("--truth-dir", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--sdf", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args(argv)
    names = [line.split("\t")[0] for line in open(args.samples) if line.strip() and not line.startswith("#")]
    all_sites = sites(args.truth_dir)
    small, sv_rows, done = {}, [], []
    for s in names:
        run = os.path.join(args.workdir, s)
        if not os.path.exists(os.path.join(run, "brca12_results.json")):
            continue
        done.append(s)
        report = json.load(open(os.path.join(run, "brca12_results.json")))
        calls = os.path.join(run, "calls", f"{s}.small.vcf.gz")
        bench = os.path.join(run, "benchmark_hprc")
        os.makedirs(bench, exist_ok=True)
        truth, bed, n = prepare_truth(args.truth_dir, s, bench)
        giab_benchmark.main(["--sample", s, "--calls", calls, "--truth", truth, "--truth-bed", bed,
                             "--regions", os.path.join(REPO, "assets", "brca12_regions.bed"),
                             "--exons", os.path.join(REPO, "assets", "brca12_exons.bed"),
                             "--sdf", args.sdf, "--outdir", bench, "--threads", str(args.threads)])
        summary = json.load(open(os.path.join(bench, "summary.json")))
        for r in summary["results"]:
            t = small.setdefault((r["region"], r["type"]), {"truthTotal": 0, "TP": 0, "FP": 0, "FN": 0, "perSample": {}})
            for k in ("truthTotal", "TP", "FP", "FN"):
                t[k] += r[k]
            t["perSample"][s] = {k: r[k] for k in ("truthTotal", "TP", "FP", "FN")}
        for site in all_sites:
            state, lengths = site_genotype(args.truth_dir, s, site)
            found = evidence(site, report, calls)
            sv_rows.append({"sample": s, "site": site["id"], "gene": site["gene"], "truth": state, "truthLengths": lengths, **found})
    by_site = {}
    for r in sv_rows:
        t = by_site.setdefault(r["site"], {"carriers": 0, "carriersCalled": 0, "carriersWithReadEndSignal": 0,
                                           "reference": 0, "referenceCalled": 0, "referenceWithReadEndSignal": 0, "unknown": 0})
        if r["truth"] == "unknown":
            t["unknown"] += 1
            continue
        key = "reference" if r["truth"] == "reference" else "carriers"
        t[key] += 1
        t[key + "Called"] += bool(r["manta"] or r["smallVariantCaller"])
        t[key + "WithReadEndSignal"] += bool(r["readEndSignal"])
    result = {"samples": done, "sites": all_sites, "bySite": by_site,
              "small": [{"region": k[0], "type": k[1], **v} for k, v in sorted(small.items())], "svRows": sv_rows}
    json.dump(result, open(args.out + ".json", "w"), indent=1)
    with open(args.out + ".tsv", "w") as handle:
        handle.write("region\ttype\ttruth\tTP\tFP\tFN\n")
        for k, v in sorted(small.items()):
            handle.write(f"{k[0]}\t{k[1]}\t{v['truthTotal']}\t{v['TP']}\t{v['FP']}\t{v['FN']}\n")
            print(k, {x: v[x] for x in ("truthTotal", "TP", "FP", "FN")})
    for site in all_sites:
        print(site["id"], site["gene"], site["lengths"], by_site.get(site["id"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())

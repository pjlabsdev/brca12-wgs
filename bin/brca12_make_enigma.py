#!/usr/bin/env python3
"""Build assets/enigma/ from the UCSC track hub of the ClinGen ENIGMA BRCA1/BRCA2 VCEP specifications.

Source: https://hgdownload.soe.ucsc.edu/hubs/enigma/ (hg38 bigBed tracks, specifications v1.2).
The hub states that the data are for research purposes only and asks users to check the version
against the ClinGen CSpec registry (GN092 BRCA1, GN097 BRCA2) before use.

Needs pyBigWig (only for this build step, not for pipeline runs):
    python3 bin/brca12_make_enigma.py --out assets/enigma [--bigbed-dir DIR]
"""
import argparse
import csv
import datetime
import hashlib
import json
import pathlib
import urllib.request

HUB = "https://hgdownload.soe.ucsc.edu/hubs/enigma"
TRACKS = {
    "BRCAclinDomains": "domains.tsv",
    "BRCAexonWeights": "exon_weights.tsv",
    "BRCAsplicing": "pvs1.tsv",
    "BRCAfunctionalAssays": "functional.tsv",
    "BRCAmfa": "multifactorial.tsv",
}


def fetch(name, directory):
    path = directory / f"{name}.bb"
    if not path.exists():
        urllib.request.urlretrieve(f"{HUB}/hg38/{name}.bb", path)
    return path


def entries(path):
    import pyBigWig
    bb = pyBigWig.open(str(path))
    for chrom, length in sorted(bb.chroms().items()):
        for start, end, rest in bb.entries(chrom, 0, length) or []:
            yield chrom, start, end, rest.split("\t")


def clean(value):
    return (value or "").replace("\t", " ").replace("\n", " ").strip()


def c_part(hgvs):
    """'NM_007294.4:c.68_69del' or 'c.68_69del' -> 'c.68_69del'."""
    hgvs = clean(hgvs)
    return hgvs.split(":", 1)[1] if ":" in hgvs else hgvs


def build(out, bigbed_dir):
    out.mkdir(parents=True, exist_ok=True)
    bigbed_dir.mkdir(parents=True, exist_ok=True)
    sources = {}
    for track, filename in TRACKS.items():
        path = fetch(track, bigbed_dir)
        sources[track] = {"url": f"{HUB}/hg38/{track}.bb",
                          "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        rows = []
        for chrom, start, end, f in entries(path):
            if track == "BRCAclinDomains":
                # name score strand thickStart thickEnd rgb gene NM aa mouseOver
                rows.append([chrom, start, end, f[6], clean(f[0]), f[7], f[8]])
            elif track == "BRCAexonWeights":
                # name score strand thickS thickE NM exon c.start c.end p.start p.end nmdBoundary pm5Strength points ...
                rows.append([chrom, start, end, clean(f[0]).split(" ")[0], f[2], clean(f[6]),
                             clean(f[7]), clean(f[8]), clean(f[11]), clean(f[12]), clean(f[13])])
            elif track == "BRCAsplicing":
                # name score strand thickS thickE rgb gene NM exon position variant code observations warnings
                rows.append([chrom, start, end, f[6], f[2], clean(f[8]), clean(f[9]), clean(f[10]),
                             clean(f[11]), clean(f[12]), clean(f[13])])
            elif track == "BRCAfunctionalAssays":
                # name score strand thickS thickE rgb NM gene p. code weight text splicePublished spliceAI splice pubs r1..r4
                rows.append([chrom, start, end, f[7], c_part(f[0]), clean(f[8]), clean(f[9]),
                             clean(f[10]), clean(f[11]), clean(f[13])])
            elif track == "BRCAmfa":
                # name score strand thickS thickE rgb LLR code ...
                rows.append([chrom, start, end, "BRCA2" if chrom == "chr13" else "BRCA1",
                             c_part(f[0]), clean(f[6]), clean(f[7])])
        header = {
            "BRCAclinDomains": ["chrom", "start0", "end", "gene", "domain", "transcript", "aa_range"],
            "BRCAexonWeights": ["chrom", "start0", "end", "gene", "strand", "exon", "c_start", "c_end",
                                "nmd_boundary_p", "pm5_ptc_strength", "pm5_points"],
            "BRCAsplicing": ["chrom", "start0", "end", "gene", "strand", "exon", "position", "variant",
                             "code", "observation", "warning"],
            "BRCAfunctionalAssays": ["chrom", "start0", "end", "gene", "hgvsc", "hgvsp", "code", "weight",
                                     "text", "spliceai"],
            "BRCAmfa": ["chrom", "start0", "end", "gene", "hgvsc", "combined_lr", "code"],
        }[track]
        with open(out / filename, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)
        sources[track]["rows"] = len(rows)
        sources[track]["file"] = filename
    meta = {
        "name": "ClinGen ENIGMA BRCA1/BRCA2 VCEP specifications",
        "version": "1.2",
        "cspec": ["https://cspec.genome.network/cspec/ui/svi/doc/GN092",
                  "https://cspec.genome.network/cspec/ui/svi/doc/GN097"],
        "hub": f"{HUB}/enigma.html",
        "note": "Research use only. Check the version against the ClinGen CSpec registry before use.",
        "built": datetime.date.today().isoformat(),
        "tracks": sources,
    }
    (out / "SOURCE.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="assets/enigma")
    parser.add_argument("--bigbed-dir", default="build/enigma-bigbed")
    args = parser.parse_args()
    meta = build(pathlib.Path(args.out), pathlib.Path(args.bigbed_dir))
    for track, info in meta["tracks"].items():
        print(f"{track}: {info['rows']} rows -> {info['file']}")


if __name__ == "__main__":
    main()

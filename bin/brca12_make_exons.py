#!/usr/bin/env python3
"""Builds assets/brca12_exons.bed from Ensembl (GRCh38, MANE Select transcripts).

Run when maintaining the pipeline, not during an analysis:
    python3 bin/brca12_make_exons.py > assets/brca12_exons.bed
Expected: BRCA1 ENST00000357654 with 23 exons, BRCA2 ENST00000380152 with 27 exons.
"""
import json
import sys
import urllib.request

TRANSCRIPTS = {"BRCA1": ("ENST00000357654", 23), "BRCA2": ("ENST00000380152", 27)}
SERVER = "https://rest.ensembl.org"


def fetch(path):
    request = urllib.request.Request(SERVER + path, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def main():
    lines = []
    for gene, (transcript, expected) in TRANSCRIPTS.items():
        data = fetch(f"/overlap/id/{transcript}?feature=exon")
        exons = [e for e in data if e.get("Parent") == transcript]
        if len(exons) != expected:
            sys.exit(f"{gene}: {len(exons)} exons instead of {expected}")
        if any(e.get("assembly_name") not in (None, "GRCh38") for e in exons):
            sys.exit(f"{gene}: Ensembl did not return GRCh38")
        strand = "-" if exons[0]["strand"] == -1 else "+"
        for exon in sorted(exons, key=lambda e: e["rank"]):
            lines.append((f"chr{exon['seq_region_name']}", exon["start"] - 1, exon["end"],
                          f"{gene}_ex{exon['rank']:02d}", exon["exon_id"], strand))
    print("# BRCA1/BRCA2 MANE Select exons, GRCh38 (0-based BED), source: Ensembl REST " + SERVER)
    for row in sorted(lines, key=lambda r: (r[0], r[1])):
        print("\t".join(str(v) for v in row))


if __name__ == "__main__":
    main()

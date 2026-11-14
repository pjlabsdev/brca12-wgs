#!/usr/bin/env python3
"""Build assets/exon_cnv_panel_1kg30x.json from the per-sample read intervals of panel.tsv
(written by fetch.py), with the parameters frozen in frozen_rules.json. Research use only.

usage: build_panel.py EXONS.bed DATADIR panel.tsv frozen_rules.json OUT.json"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bin"))
sys.path.insert(0, HERE)
from brca12lib import exon_cnv, coverage  # noqa: E402
import cnvmodel  # noqa: E402
import fetch_controls  # noqa: E402

bed, data, panel_list, frozen, out = sys.argv[1:6]
rules = json.load(open(frozen))["parameters"]
pad = rules["pad_bp"]
exons = {}
for line in open(bed):
    if not line.strip() or line.startswith("#"):
        continue
    chrom, s, e, name = line.split("\t")[:4]
    gene, number, label = coverage.exon_label(name)
    exons[name] = {"chrom": chrom, "start": int(s), "end": int(e), "gene": gene, "label": label}
samples = []
for sid in (l.split("\t")[0] for l in open(panel_list) if l.strip()):
    d = json.load(open(os.path.join(data, sid + ".json")))
    depths = {n: cnvmodel.exon_depth(d["exons"][n], e["start"] - pad, e["end"] + pad) for n, e in exons.items()}
    depths.update(d["controls"])
    samples.append(depths)
panel = exon_cnv.build_panel(samples, exons, fetch_controls.CONTROLS, pad, rules["p_change"], rules["min_q"],
                             rules["min_sd_log2"],
                             description="80 unrelated 1000 Genomes 30x samples (NYGC, GRCh38, BWA-MEM); "
                                         "valid for 30x PCR-free Illumina WGS aligned the same way")
json.dump(panel, open(out, "w"), indent=1)
print(len(samples), "samples,", len(panel["exons"]), "exons ->", out)

#!/usr/bin/env python3
"""Build the application note (HTML with inline SVG figures) from the validation result files.

Every number and figure in the note is read from the files under validation/ so that the text
cannot drift from the data. Usage:

    python3 docs/paper/build_paper.py [--fonts DIR] > docs/paper/brca12-wgs-application-note.html

--fonts points to a directory with the IBM Plex Mono and Instrument Sans woff2 files (fontsource
layout); they are embedded so that the PDF does not depend on installed fonts.
"""
import argparse
import base64
import glob
import gzip
import html
import json
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
VAL = os.path.join(ROOT, "validation")
PATHO = os.path.join(VAL, "pathogenic")
RES = os.path.join(PATHO, "results")

INK, SUB, LINE, PAPER = "#141413", "#56544C", "#D8D5CC", "#FAF9F5"
GOLD = "#A57A27"                       # brand mark only, never a chart colour
OK, BAD, WARN, BLUE, GREY = "#2A7350", "#B3261E", "#8F5B00", "#2F5D8A", "#9A978E"
MONO = "'IBM Plex Mono', ui-monospace, Menlo, monospace"
SANS = "'Instrument Sans', 'Helvetica Neue', Arial, sans-serif"

LOGO = """<svg xmlns="http://www.w3.org/2000/svg" height="26" viewBox="0 0 600 95" role="img" aria-label="PJ Labs">
<g><path d="M 60,18.5 L 70,23.1 L 80,29.4 L 90,36.5 L 100,44.2 L 110,51.8 L 120,59.2 L 130,65.8 L 140,71.4 L 150,75.6 L 160,78.2 L 170,79.1 L 180,78.2 L 190,75.6 L 200,71.4 L 210,65.8 L 220,59.2 L 230,51.8 L 240,44.2 L 250,36.5" fill="none" stroke="#A57A27" stroke-width="13" stroke-linecap="round" stroke-linejoin="round"/>
<path d="M 80,70.6 L 90,63.5 L 100,55.8 L 110,48.2 L 120,40.8 L 130,34.2 L 140,28.6 L 150,24.4 L 160,21.8 L 170,20.9 L 180,21.8 L 190,24.4 L 200,28.6 L 210,34.2 L 220,40.8 L 230,48.2 L 240,55.8 L 250,63.5 L 260,70.6 L 270,76.9 L 280,81.5" fill="none" stroke="#A57A27" stroke-width="13" stroke-linecap="round" stroke-linejoin="round" opacity="0.35"/>
<line x1="144" y1="40.92" x2="144" y2="59.08" stroke="#A57A27" stroke-width="5" stroke-linecap="round"/>
<line x1="170" y1="34.90" x2="170" y2="65.10" stroke="#A57A27" stroke-width="5" stroke-linecap="round"/>
<line x1="196" y1="40.92" x2="196" y2="59.08" stroke="#A57A27" stroke-width="5" stroke-linecap="round"/></g>
<text x="315" y="72" font-family="'IBM Plex Mono', ui-monospace, Menlo, monospace" font-weight="500" font-size="46" fill="#141413">PJ Labs</text></svg>"""


REFS = {
    'bwa': 'Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997; 2013.',
    'samtools': 'Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008.',
    'gatk': 'Poplin R, Ruano-Rubio V, DePristo MA, et al. Scaling accurate genetic variant discovery to tens of thousands of samples. bioRxiv 201178; 2018.',
    'deepvariant': 'Poplin R, Chang PC, Alexander D, et al. A universal SNP and small-indel variant caller using deep neural networks. Nat Biotechnol. 2018;36:983–987.',
    'manta': 'Chen X, Schulz-Trieglaff O, Shaw R, et al. Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics. 2016;32(8):1220–1222.',
    'cnvkit': 'Talevich E, Shain AH, Botton T, Bastian BC. CNVkit: genome-wide copy number detection and visualization from targeted DNA sequencing. PLoS Comput Biol. 2016;12(4):e1004873.',
    'vep': 'McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biol. 2016;17:122.',
    'mosdepth': 'Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics. 2018;34(5):867–868.',
    'nextflow': 'Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017;35:316–319.',
    'acmg': 'Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants. Genet Med. 2015;17(5):405–424.',
    'clinvar': 'Landrum MJ, Lee JM, Benson M, et al. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Res. 2018;46(D1):D1062–D1067.',
    'mane': 'Morales J, Pujar S, Loveland JE, et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. Nature. 2022;604:310–315.',
    'gnomad': 'Karczewski KJ, Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature. 2020;581:434–443.',
    'brcaanalyzer': 'Kechin A. BRCA-analyzer: automatic workflow for analysis of BRCA1/2 NGS data. https://github.com/aakechin/BRCA-analyzer (accessed 12 November 2026).',
    'enigma': 'Parsons MT, de la Hoya M, Richardson ME, et al. Evidence-based recommendations for gene-specific ACMG/AMP variant classification from the ClinGen ENIGMA BRCA1 and BRCA2 Variant Curation Expert Panel. Am J Hum Genet. 2024;111(9):2044–2058.',
    'pvs1': 'Abou Tayoun AN, Pesaran T, DiStefano MT, et al. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. Hum Mutat. 2018;39(11):1517–1524.',
    'giab1': 'Zook JM, Chapman B, Wang J, et al. Integrating human sequence data sets provides a resource of benchmark SNP and indel genotype calls. Nat Biotechnol. 2014;32:246–251.',
    'giab2': 'Wagner J, Olson ND, Harris L, et al. Benchmarking challenging small variants with linked and long reads. Cell Genomics. 2022;2(5):100128.',
    'vcfeval': 'Cleary JG, Braithwaite R, Gaastra K, et al. Comparing variant call files for performance benchmarking of next-generation sequencing variant calling pipelines. bioRxiv 023754; 2015.',
    'kgp': 'Byrska-Bishop M, Evani US, Zhao X, et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. Cell. 2022;185(18):3426–3440.',
    'getrm': 'Barker SD, Bale S, Booker J, et al. Development and characterization of reference materials for MTHFR, SERPINA1, RET, BRCA1, and BRCA2 genetic testing. J Mol Diagn. 2009;11(6):553–561.',
    'kuchenbaecker': 'Kuchenbaecker KB, Hopper JL, Barnes DR, et al. Risks of breast, ovarian, and contralateral breast cancer for BRCA1 and BRCA2 mutation carriers. JAMA. 2017;317(23):2402–2416.',
    'sluiter': 'Sluiter MD, van Rensburg EJ. Large genomic rearrangements of the BRCA1 and BRCA2 genes: review of the literature and report of a novel BRCA1 mutation. Breast Cancer Res Treat. 2011;125(2):325–349.',
    'teugels': 'Teugels E, De Brakeleer S, Goelen G, Lissens W, Sermijn E, De Grève J. De novo Alu element insertions targeted to a sequence common to the BRCA1 and BRCA2 genes. Hum Mutat. 2005;26(3):284.',
    'machado': 'Machado PM, Brandão RD, Cavaco BM, et al. Screening for a BRCA2 rearrangement in high-risk breast/ovarian cancer families: evidence for a founder effect and analysis of the associated phenotypes. J Clin Oncol. 2007;25(15):2027–2034.',
    'sarek': 'Garcia M, Juhos S, Larsson M, et al. Sarek: a portable workflow for whole-genome sequencing analysis of germline and somatic variants. F1000Research. 2020;9:63.',
    'spliceai': 'Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. Predicting splicing from primary sequence with deep learning. Cell. 2019;176(3):535–548.',
    'bayesdel': 'Feng BJ. PERCH: a unified framework for disease gene prioritization. Hum Mutat. 2017;38(3):243–251.',
    'mlpa': 'Schouten JP, McElgunn CJ, Waaijer R, Zwijnenburg D, Diepvens F, Pals G. Relative quantification of 40 nucleic acid sequences by multiplex ligation-dependent probe amplification. Nucleic Acids Res. 2002;30(12):e57.',
    'melt': 'Gardner EJ, Lam VK, Harris DN, et al. The Mobile Element Locator Tool (MELT): population-scale mobile element discovery and biology. Genome Res. 2017;27(11):1916–1929.',
    'gatksv': 'Collins RL, Brand H, Karczewski KJ, et al. A structural variation reference for medical and population genetics. Nature. 2020;581:444–451.',
    'mutect2': 'Benjamin D, Sato T, Cibulskis K, Getz G, Stewart C, Lichtenstein L. Calling somatic SNVs and indels with Mutect2. bioRxiv 861054; 2019.',
    'giabdata': 'Zook JM, Catoe D, McDaniel J, et al. Extensive sequencing of seven human genomes to characterize benchmark reference materials. Sci Data. 2016;3:160025.',
    'hprc': 'Liao WW, Asri M, Ebler J, et al. A draft human pangenome reference. Nature. 2023;617:312–324.',
    'mgcactus': 'Hickey G, Monlong J, Ebler J, et al. Pangenome graph construction from genome alignments with Minigraph-Cactus. Nat Biotechnol. 2024;42:663–673.',
    'ont1kg': 'Schloissnig S, Pani S, Ebler J, et al. Structural variation in 1,019 diverse humans based on long-read sequencing. Nature. 2025;644(8076):442–452.',
    'wilson': 'Wilson EB. Probable inference, the law of succession, and statistical inference. J Am Stat Assoc. 1927;22(158):209–212.',
}


# ------------------------------------------------------------------ data
def load(path):
    with open(path) as handle:
        return json.load(handle)


def table(path):
    rows, keys = [], None
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                keys = line[1:].rstrip("\n").split("\t")
            elif line.strip():
                rows.append(dict(zip(keys, line.rstrip("\n").split("\t"))))
    return rows


def tsv(path):
    with open(path) as handle:
        lines = [l.rstrip("\n").split("\t") for l in handle if l.strip()]
    return [dict(zip(lines[0], l)) for l in lines[1:]]


TARGETED = (("HG002_exome_agilent", "HG002, Agilent SureSelect v5 exome"),
            ("HG001_exome_truseq", "HG001, Illumina TruSeq exome"),
            ("HG002_ampliseq_ion", "HG002, Ion AmpliSeq exome (single-end)"))


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def exons():
    out = {"BRCA1": [], "BRCA2": []}
    for line in open(os.path.join(ROOT, "assets", "brca12_exons.bed")):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and "_ex" in f[3]:
            gene, number = f[3].split("_ex")
            out[gene].append((int(number), int(f[1]) + 1, int(f[2])))
    return out


class Data:
    def __init__(self):
        self.giab = [r for s in ("HG001", "HG002") for r in load(os.path.join(VAL, "results", f"{s}.summary.json"))["results"]]
        self.carriers = load(os.path.join(RES, "carriers_summary.json"))
        self.spike = load(os.path.join(RES, "spike_summary.json"))
        self.spike_before = load(os.path.join(RES, "spike_summary.before_fix.json"))
        self.sv = load(os.path.join(RES, "sv_carriers_summary.json"))
        self.small_events = table(os.path.join(PATHO, "small_events.tsv"))
        self.large_events = table(os.path.join(PATHO, "large_events.tsv"))
        self.carrier_table = table(os.path.join(PATHO, "carriers_1000g.tsv"))
        self.sv_sites = table(os.path.join(PATHO, "sv_sites_1000g.tsv"))
        self.depth = table(os.path.join(RES, "depth_bins.tsv"))
        self.exons = exons()
        self.targeted = {n: {"bench": tsv(os.path.join(VAL, "targeted", "results", f"{n}.benchmark.tsv")),
                             "qc": load(os.path.join(VAL, "targeted", "results", f"{n}.qc.json"))} for n, _ in TARGETED}
        self.mosaic = load(os.path.join(RES, "mosaic_summary.json"))
        self.hprc = load(os.path.join(VAL, "hprc", "results", "hprc_summary.json"))
        self.ont = tsv(os.path.join(VAL, "longread", "results", "ours_vs_shortread.tsv"))
        self.tandem_hprc = tsv(os.path.join(VAL, "hprc", "results", "tandem_copy.tsv"))
        self.tandem_ont = tsv(os.path.join(VAL, "longread", "results", "tandem_copy.tsv"))

        disc = os.path.join(VAL, "longread", "results", "all_disc.json")
        self.ont_all = load(disc) if os.path.exists(disc) else None

    @staticmethod
    def tandem(rows, truth):
        sub = [r for r in rows if r["truth"] == truth]
        return sum(r["signature"] == "yes" for r in sub), len(sub)

    def target_row(self, name, kind):
        r = next(x for x in self.targeted[name]["bench"] if x["region"] == "targets" and x["type"] == kind)
        return int(r["TP"]), int(r["truthTotal"]), int(r["FP"])

    def hprc_small(self, region, kind):
        r = next(x for x in self.hprc["small"] if x["region"] == region and x["type"] == kind)
        return r["TP"], r["truthTotal"], r["FP"], r["FN"]

    def giab_total(self, region, kind):
        rows = [r for r in self.giab if r["region"] == region and r["type"] == kind]
        return sum(r["TP"] for r in rows), sum(r["truthTotal"] for r in rows), sum(r["FP"] for r in rows)

    def concordance(self, region, kind):
        r = next(x for x in self.carriers["concordanceWithNygc"] if x["region"] == region and x["type"] == kind)
        return r["TP"], r["nygc"], r["FP"]

    def large(self, which="spike"):
        return (self.spike if which == "spike" else self.spike_before)["largeEvents"]

    def count_large(self, key, which="spike"):
        return sum(1 for e in self.large(which) if e.get(key))

    def sv_site(self, sid):
        return self.sv["bySite"][sid]


# ------------------------------------------------------------------ svg helpers
def esc(text):
    return html.escape(str(text), quote=True)


def text(x, y, s, size=8, anchor="start", fill=INK, mono=False, weight=400, italic=False, rotate=None):
    style = f"font-family:{MONO if mono else SANS};font-size:{size}px;font-weight:{weight};" + ("font-style:italic;" if italic else "")
    tr = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    return f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" fill="{fill}" style="{style}"{tr}>{esc(s)}</text>'


def svg(width, height, body):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" '
            f'role="img">{"".join(body)}</svg>')


def fmt_bp(n):
    if n >= 1000:
        v = n / 1000
        return f"{v:.0f} kb" if v >= 10 or abs(v - round(v)) < 0.05 else f"{v:.1f} kb"
    return f"{n} bp"


# ------------------------------------------------------------------ figure 1: workflow
def fig_workflow():
    w, h = 720, 250
    b = []

    def box(x, y, bw, bh, title, sub, fill="#FFFFFF", stroke=LINE):
        b.append(f'<rect x="{x}" y="{y}" width="{bw}" height="{bh}" rx="5" fill="{fill}" stroke="{stroke}"/>')
        b.append(text(x + bw / 2, y + 15, title, 8.5, "middle", mono=True, weight=500))
        for i, line in enumerate(sub):
            b.append(text(x + bw / 2, y + 27 + 10 * i, line, 7.2, "middle", SUB))

    def arrow(x1, y1, x2, y2):
        b.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{SUB}" stroke-width="1" marker-end="url(#ah)"/>')

    b.append(f'<defs><marker id="ah" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="6" markerHeight="6" orient="auto">'
             f'<path d="M0,0 L8,4 L0,8 z" fill="{SUB}"/></marker></defs>')
    box(8, 96, 100, 48, "Input", ["FASTQ or BAM/CRAM", "WGS, exome, panel"])
    box(132, 96, 104, 48, "Check + align", ["GRCh38 only", "BWA-MEM, samtools"])
    arrow(108, 120, 131, 120)
    rows = [("Small variants", ["GATK HaplotypeCaller", "or DeepVariant"], 8),
            ("Structural variants", ["Manta, genes ±1 Mb"], 68),
            ("Copy number", ["CNVkit +", "exon depth ratios"], 128),
            ("Coverage", ["mosdepth per exon,", "MAPQ-aware"], 188)]
    for title, sub, y in rows:
        box(268, y, 128, 48, title, sub)
        arrow(236, 120, 267, y + 24)
    box(428, 8, 128, 48, "Annotation", ["VEP 113, MANE Select,", "ClinVar"])
    arrow(396, 32, 427, 32)
    box(588, 8, 124, 48, "ACMG/AMP hints", ["ClinGen ENIGMA v1.2,", "single criteria only"])
    arrow(556, 32, 587, 32)
    box(588, 128, 124, 48, "Report", ["HTML, Markdown, JSON", "no classification"], fill="#F3F1EA")
    arrow(650, 56, 650, 127)
    for y, ty in ((92, 140), (152, 152), (212, 164)):
        arrow(396, y, 587, ty)
    b.append(f'<rect x="8" y="176" width="228" height="62" rx="5" fill="none" stroke="{WARN}" stroke-dasharray="3 2"/>')
    b.append(text(18, 192, "Research use only", 8.5, mono=True, weight=500, fill=WARN))
    b.append(text(18, 205, "Findings are candidates. Confirmation in an", 7.2, fill=SUB))
    b.append(text(18, 215, "accredited laboratory (Sanger, MLPA) is required", 7.2, fill=SUB))
    b.append(text(18, 225, "before any medical interpretation.", 7.2, fill=SUB))
    return svg(w, h, b)


# ------------------------------------------------------------------ figure 2: gene maps
def fig_gene_maps(d):
    w, h = 720, 430
    b = []
    left, right = 62, 706
    panels = [("BRCA1", "chr17", 43044295, 43170327, "minus strand, NM_007294.4", 20),
              ("BRCA2", "chr13", 32315077, 32400268, "plus strand, NM_000059.4", 232)]
    for gene, chrom, lo, hi, label, top in panels:
        def X(g):
            return left + (g - lo) / (hi - lo) * (right - left)
        axis = top + 78
        b.append(text(8, top + 8, gene, 10, mono=True, weight=500, italic=True))
        b.append(text(left, top + 8, f"{chrom}:{lo:,}–{hi:,} · {label}", 7.2, fill=SUB))
        # real carriers (lollipops)
        seen = {}
        for c in d.carrier_table:
            if c["gene"] == gene:
                seen.setdefault((int(c["pos"]), c["hgvs_c"]), 0)
                seen[(int(c["pos"]), c["hgvs_c"])] += 1
        placed = []
        for i, ((pos, name), n) in enumerate(sorted(seen.items())):
            x = X(pos)
            level = 0
            while any(abs(x - px) < 62 and level == pl for px, pl in placed):
                level += 1
            placed.append((x, level))
            y = axis - 26 - 14 * level
            b.append(f'<line x1="{x:.1f}" y1="{axis - 7}" x2="{x:.1f}" y2="{y + 3:.1f}" stroke="{BAD}" stroke-width="0.9"/>')
            b.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.2" fill="{BAD}"/>')
            anchor = "end" if x > right - 90 else "start"
            b.append(text(x + (5 if anchor == "start" else -5), y + 2.6, name + (f" (×{n})" if n > 1 else ""), 6.8, anchor, mono=True))
        # gene body and exons
        b.append(f'<line x1="{left}" y1="{axis}" x2="{right}" y2="{axis}" stroke="{INK}" stroke-width="1"/>')
        for number, start, end in d.exons[gene]:
            x1, x2 = X(start), X(end)
            b.append(f'<rect x="{x1:.1f}" y="{axis - 6}" width="{max(1.3, x2 - x1):.1f}" height="12" fill="{INK}"/>')
        first = min(d.exons[gene], key=lambda e: e[0])
        last = max(d.exons[gene], key=lambda e: e[0])
        for number, start, end in (first, last):
            b.append(text(X((start + end) / 2), axis + 16, f"ex {number}", 6.5, "middle", SUB, mono=True))
        # spiked small variants
        y0 = axis + 24
        for e in d.small_events:
            if e["gene"] == gene:
                x = X(int(e["pos"]))
                colour = BLUE if e["kind"] == "SNV" else (OK if e["kind"] == "INS" else WARN)
                b.append(f'<line x1="{x:.1f}" y1="{y0}" x2="{x:.1f}" y2="{y0 + 9}" stroke="{colour}" stroke-width="0.8"/>')
        b.append(text(8, y0 + 8, "spike-in", 6.8, fill=SUB))
        # large spiked events
        y1 = y0 + 18
        rows = [e for e in d.large_events if e["gene"] == gene]
        rows.sort(key=lambda e: -int(e["size"]))
        for i, e in enumerate(rows):
            y = y1 + 9 * i
            start, end = int(e["start"]), int(e["end"])
            x1, x2 = X(max(lo - (hi - lo) * 0.03, start)), X(min(hi + (hi - lo) * 0.03, end))
            colour = OK if e["kind"] == "LDUP" else (BAD if e["kind"] == "LDEL" else INK)
            if e["kind"] == "ALU":
                b.append(f'<path d="M{x1 - 3:.1f},{y:.1f} L{x1 + 3:.1f},{y:.1f} L{x1:.1f},{y + 6:.1f} z" fill="{INK}"/>')
                lab = "Alu insertion, exon 3"
            else:
                b.append(f'<rect x="{x1:.1f}" y="{y:.1f}" width="{max(2.2, x2 - x1):.1f}" height="5" fill="{colour}" opacity="0.8"/>')
                lab = ("dup " if e["kind"] == "LDUP" else "del ") + fmt_bp(int(e["size"]))
            lx = x2 + 4 if x2 < right - 80 else x1 - 4
            b.append(text(lx, y + 5, lab, 6.5, "start" if x2 < right - 80 else "end", SUB, mono=True))
        # real SVs
        for s in d.sv_sites:
            if s["gene"] == gene:
                x = X(int(s["pos"]))
                b.append(f'<path d="M{x:.1f},{axis + 8} l-3.5,7 l7,0 z" fill="none" stroke="{INK}" stroke-width="0.9"/>')
    ly = h - 14
    items = [("circle", BAD, "pathogenic variant, real carrier (1000 Genomes)"), ("tick", BLUE, "spiked SNV"),
             ("tick", WARN, "spiked deletion ≤ 40 bp"), ("tick", OK, "spiked insertion ≤ 25 bp"),
             ("bar", BAD, "spiked large deletion"), ("bar", OK, "spiked tandem duplication"),
             ("tri", INK, "real SV (1000 Genomes)")]
    x = 8
    for kind, colour, label in items:
        if kind == "circle":
            b.append(f'<circle cx="{x + 3}" cy="{ly - 3}" r="3" fill="{colour}"/>')
        elif kind == "tick":
            b.append(f'<line x1="{x + 3}" y1="{ly - 7}" x2="{x + 3}" y2="{ly + 1}" stroke="{colour}" stroke-width="1.2"/>')
        elif kind == "bar":
            b.append(f'<rect x="{x}" y="{ly - 5}" width="8" height="4" fill="{colour}" opacity="0.8"/>')
        else:
            b.append(f'<path d="M{x + 4},{ly - 7} l-3.5,7 l7,0 z" fill="none" stroke="{colour}" stroke-width="0.9"/>')
        b.append(text(x + 11, ly, label, 6.6, fill=SUB))
        x += 17 + len(label) * 3.25
    return svg(w, h, b)


# ------------------------------------------------------------------ figure 3: recall by class
def recall_rows(d):
    rows = []
    k, n, _ = d.giab_total("gene_loci", "SNV")
    rows.append(("GIAB truth sets (HG001, HG002)", "SNVs, gene loci", k, n))
    k, n, _ = d.giab_total("gene_loci", "INDEL")
    rows.append(("", "Indels, gene loci", k, n))
    ks, ns, _ = d.giab_total("exons_pm20", "SNV")
    ki, ni, _ = d.giab_total("exons_pm20", "INDEL")
    rows.append(("", "All variants, exons ±20 bp", ks + ki, ns + ni))
    rows.append(("Real carriers (1000 Genomes)", "Pathogenic SNVs/indels", d.carriers["detected"], d.carriers["total"]))
    order = ["SNV", "deletion 1 bp", "deletion 2-5 bp", "deletion 6-20 bp", "deletion 21-50 bp",
             "insertion/duplication 1 bp", "insertion/duplication 2-5 bp", "insertion/duplication 6-20 bp",
             "insertion/duplication 21-50 bp"]
    for i, key in enumerate(order):
        v = d.spike["smallByClass"][key]
        rows.append(("Spike-in, small variants" if i == 0 else "", key[0].upper() + key[1:].replace("-", "–") + ("s" if key == "SNV" else ""), v["found"], v["spiked"]))
    rows.append(("Spike-in, large events", "Deletions, duplications, Alu (Manta)", d.count_large("manta"), len(d.large())))
    rows.append(("", "Same events, CNVkit segments", d.count_large("cnvkit"), len(d.large())))
    for sid, label in (("HGSV_216998", "Deletion 1,081 bp (BRCA1)"), ("HGSV_216996", "Deletion 271 bp (BRCA1)"),
                       ("HGSV_181048", "Alu insertion, common (BRCA2)"), ("HGSV_181052", "Alu insertion, rare (BRCA2)")):
        s = d.sv_site(sid)
        alu = "Alu" in label
        rows.append(("Real SVs (1000 Genomes)" if sid == "HGSV_216998" else "", label + (", called" if alu else ""), s["carriersReported"], s["carriers"]))
        if alu:
            rows.append(("", label + ", read-end signature", s["carriersWithReadEndSignal"], s["carriers"]))
    for i, (name, label) in enumerate(TARGETED):
        k, n, _ = d.target_row(name, "SNV")
        rows.append(("Exome and amplicon (GIAB)" if i == 0 else "", label.split(", ", 1)[1].replace(" (single-end)", "") + ", SNVs in targets", k, n))
    for i, (region, kind, label) in enumerate((("exons_pm20", "SNV", "SNVs, exons ±20 bp"), ("exons_pm20", "INDEL", "Indels, exons ±20 bp"),
                                               ("gene_loci", "SNV", "SNVs, gene loci"), ("gene_loci", "INDEL", "Indels, gene loci"))):
        k, n, _, _ = d.hprc_small(region, kind)
        rows.append(("Assembly truth (39 HPRC samples)" if i == 0 else "", label, k, n))
    for sid, label, key in (("chr13:32348295", "Alu insertion, read-end signature", "carriersWithReadEndSignal"),
                            ("chr17:43100595", "VNTR length change, called", "carriersCalled"),
                            ("chr17:43128264", "Tandem-repeat insertion 126 bp, called", "carriersCalled")):
        s = d.hprc["bySite"][sid]
        rows.append(("", label, s[key], s["carriers"]))
    rows.append(("", "Same insertion, tandem-copy signature", *d.tandem(d.tandem_hprc, "carrier")))
    rows.append(("Long reads (1KG-ONT)", "Same insertion, tandem-copy signature", *d.tandem(d.tandem_ont, "carrier")))
    for i, m in enumerate(r for r in d.mosaic if r["spikedEvents"]):
        rows.append(("Spike-in, low allele fraction" if i == 0 else "", f"{m['spikedAlleleFraction'] * 100:.0f} % of reads, any call", m["reportedEither"], m["spikedEvents"]))
    return rows


def fig_recall(d):
    rows = recall_rows(d)
    w, rh, top = 720, 15.5, 26
    h = top + rh * len(rows) + 30
    b = []
    x0, x1 = 400, 630

    def X(p):
        return x0 + p * (x1 - x0)

    for p in (0, 0.25, 0.5, 0.75, 1.0):
        b.append(f'<line x1="{X(p):.1f}" y1="{top - 6}" x2="{X(p):.1f}" y2="{top + rh * len(rows)}" stroke="{LINE}" stroke-width="0.6"/>')
        b.append(text(X(p), top + rh * len(rows) + 12, f"{p:.2f}".rstrip("0").rstrip(".") if p not in (0, 1) else str(int(p)), 7, "middle", SUB, mono=True))
    b.append(text((x0 + x1) / 2, top + rh * len(rows) + 25, "Fraction detected (95 % Wilson interval)", 7.5, "middle", SUB))
    b.append(text(8, top - 12, "Data set", 7.5, fill=SUB))
    b.append(text(196, top - 12, "Variant class", 7.5, fill=SUB))
    b.append(text(706, top - 12, "found / total", 7.5, "end", SUB))
    for i, (group, label, k, n) in enumerate(rows):
        y = top + rh * i + rh / 2
        if group:
            if i:
                b.append(f'<line x1="8" y1="{y - rh / 2:.1f}" x2="706" y2="{y - rh / 2:.1f}" stroke="{LINE}" stroke-width="0.6"/>')
            b.append(text(8, y + 2.6, group, 7.6, weight=500))
        b.append(text(196, y + 2.6, label, 7.6))
        lo, hi = wilson(k, n)
        colour = OK if k == n else (WARN if k / n >= 0.5 else BAD)
        b.append(f'<line x1="{X(lo):.1f}" y1="{y:.1f}" x2="{X(hi):.1f}" y2="{y:.1f}" stroke="{colour}" stroke-width="1.4"/>')
        b.append(f'<circle cx="{X(k / n):.1f}" cy="{y:.1f}" r="3" fill="{colour}"/>')
        b.append(text(706, y + 2.6, f"{k} / {n}", 7.6, "end", mono=True))
    return svg(w, h, b)


# ------------------------------------------------------------------ figure 7: low allele fractions
def fig_mosaic(d):
    mos = [m for m in d.mosaic if m["spikedEvents"]]
    w, top, bh, gap = 720, 18, 15, 7
    x0, x1 = 90, 600
    h = top + len(mos) * (bh + gap) + 34
    b = []
    for p in (0, 0.25, 0.5, 0.75, 1.0):
        x = x0 + p * (x1 - x0)
        b.append(f'<line x1="{x:.1f}" y1="{top - 4}" x2="{x:.1f}" y2="{top + len(mos) * (bh + gap)}" stroke="{LINE}" stroke-width="0.6"/>')
        b.append(text(x, top + len(mos) * (bh + gap) + 11, f"{p * 100:.0f} %", 7, "middle", SUB, mono=True))
    b.append(text((x0 + x1) / 2, top + len(mos) * (bh + gap) + 24, "Share of spiked variants", 7.5, "middle", SUB))
    b.append(text(8, top - 6, "Allele fraction", 7.5, fill=SUB))
    for i, m in enumerate(mos):
        y = top + i * (bh + gap)
        n = m["spikedEvents"]
        parts = ((m["reportedAsMosaicCandidate"], BLUE), (m["calledByGermlineCaller"], OK), (n - m["reportedEither"], GREY))
        x = x0
        for k, colour in parts:
            wd = k / n * (x1 - x0)
            if wd > 0:
                b.append(f'<rect x="{x:.1f}" y="{y}" width="{wd:.1f}" height="{bh}" fill="{colour}"/>')
                if wd > 18:
                    b.append(text(x + wd / 2, y + bh / 2 + 2.6, str(k), 7.2, "middle", "#FFFFFF", mono=True))
            x += wd
        b.append(text(8, y + bh / 2 + 2.6, f"{m['spikedAlleleFraction'] * 100:.0f} %", 7.6, mono=True))
        b.append(text(x1 + 8, y + bh / 2 + 2.6, f"{m['reportedEither']} / {n} reported", 7.4, mono=True))
    lx = x0
    for label, colour in (("mosaic candidate", BLUE), ("germline caller", OK), ("not reported", GREY)):
        b.append(f'<rect x="{lx}" y="{h - 9}" width="8" height="8" fill="{colour}"/>')
        b.append(text(lx + 12, h - 2, label, 7.2))
        lx += 110
    return svg(w, h + 2, b)


# ------------------------------------------------------------------ figure 4: large events per caller
def fig_large(d):
    events = sorted(d.large(), key=lambda e: e["size"])
    before = {(e["sample"], e["gene"]): e for e in d.large("before")}
    cols = [("Manta", "manta"), ("CNVkit 1.0.2", "cnv_before"), ("CNVkit 1.1.0", "cnvkit"),
            ("Exon depth hint", "hint"), ("Small-variant caller", "smallVariantCaller")]
    w, rh, top = 720, 16, 44
    h = top + rh * len(events) + 34
    b = []
    cx0, cw = 392, 64
    b.append(text(8, top - 10, "Spiked event (heterozygous)", 7.5, fill=SUB))
    b.append(text(300, top - 10, "Size", 7.5, "end", SUB))
    for j, (name, _) in enumerate(cols):
        words = name.split(" ")
        if len(words) > 2:
            words = [" ".join(words[:-1]), words[-1]]
        for i, wd in enumerate(words if len(name) > 11 else [name]):
            b.append(text(cx0 + cw * j, top - 20 + 10 * i + (10 if len(name) <= 11 else 0), wd, 7.2, "middle", SUB))
    for i, e in enumerate(events):
        y = top + rh * i + rh / 2
        if i % 2 == 0:
            b.append(f'<rect x="4" y="{y - rh / 2:.1f}" width="{w - 8}" height="{rh}" fill="#F3F1EA"/>')
        kind = {"LDEL": "Deletion", "LDUP": "Tandem duplication", "ALU": "Alu insertion"}[e["kind"]]
        note = e["note"].split(",")[0].split(";")[0]
        note = re.sub(r" \(legacy exon \d+\)", "", note)
        note = note.replace("deletion of ", "").replace("tandem duplication of ", "").replace("Alu insertion in ", "")
        b.append(text(8, y + 2.6, f"{e['gene']}  {kind}", 7.6, mono=True))
        b.append(text(150, y + 2.6, note[:44], 7.2, fill=SUB))
        b.append(text(300, y + 2.6, fmt_bp(e["size"]), 7.6, "end", mono=True))
        # size bar on log scale
        bx0, bx1 = 306, 356
        frac = (math.log10(e["size"]) - 2) / (math.log10(135000) - 2)
        b.append(f'<rect x="{bx0}" y="{y - 2:.1f}" width="{max(1.5, frac * (bx1 - bx0)):.1f}" height="4" fill="{GREY}"/>')
        ratios = [v for v in e.get("affectedExonDepthRatios", {}).values() if v is not None]
        deviates = [v for v in ratios if (v < 0.65 if e["kind"] == "LDEL" else v > 1.35)]
        state = {"manta": bool(e["manta"]), "cnvkit": bool(e["cnvkit"]),
                 "cnv_before": bool(before.get((e["sample"], e["gene"]), {}).get("cnvkit")),
                 "hint": bool(e["exonDosageHints"]), "smallVariantCaller": bool(e["smallVariantCaller"])}
        applicable = {"manta": True, "cnvkit": e["kind"] != "ALU", "cnv_before": e["kind"] != "ALU",
                      "hint": bool(ratios) and e["kind"] != "ALU", "smallVariantCaller": e["size"] <= 400}
        for j, (_, key) in enumerate(cols):
            x = cx0 + cw * j
            if not applicable[key]:
                b.append(text(x, y + 2.6, "n/a", 6.8, "middle", GREY))
            elif state[key]:
                b.append(f'<circle cx="{x}" cy="{y:.1f}" r="4.2" fill="{OK}"/>')
                b.append(f'<path d="M{x - 2.2},{y:.1f} l1.6,1.8 l2.9,-3.4" fill="none" stroke="#fff" stroke-width="1.1"/>')
            else:
                b.append(f'<circle cx="{x}" cy="{y:.1f}" r="4.2" fill="none" stroke="{BAD}" stroke-width="1.1"/>')
    yb = top + rh * len(events) + 14
    for j, (_, key) in enumerate(cols):
        x = cx0 + cw * j
        if key == "cnv_before":
            k = sum(1 for e in events if before.get((e["sample"], e["gene"]), {}).get("cnvkit"))
        elif key == "hint":
            k = sum(1 for e in events if e["exonDosageHints"])
        else:
            k = sum(1 for e in events if e.get(key))
        b.append(text(x, yb, f"{k} / {len(events)}", 7.6, "middle", mono=True, weight=500))
    b.append(text(8, yb, "Detected", 7.6, weight=500))
    b.append(text(8, yb + 13, "Filled: reported with the correct type at the spiked position. Open: not reported. Size bar on a log scale.", 6.8, fill=SUB))
    return svg(w, h, b)


# ------------------------------------------------------------------ figure 5: depth profiles
def fig_depth(d):
    panels = []
    for key in sorted({r["panel"] for r in d.depth}):
        rows = [r for r in d.depth if r["panel"] == key]
        panels.append((key, rows))
    w, ph = 720, 118
    h = ph * len(panels) + 22
    b = []
    left, right = 46, 706
    for i, (key, rows) in enumerate(panels):
        top = 12 + ph * i
        lo, hi = min(int(r["start"]) for r in rows), max(int(r["end"]) for r in rows)
        y0, y1 = top + 14, top + ph - 26

        def X(g):
            return left + (g - lo) / (hi - lo) * (right - left)

        def Y(v):
            v = max(-2.0, min(1.2, v))
            return y1 - (v + 2.0) / 3.2 * (y1 - y0)

        title = rows[0]["title"]
        b.append(text(left, top + 6, title, 7.6, weight=500))
        for v in (-2, -1, 0, 1):
            b.append(f'<line x1="{left}" y1="{Y(v):.1f}" x2="{right}" y2="{Y(v):.1f}" stroke="{LINE}" stroke-width="0.6"/>')
            b.append(text(left - 5, Y(v) + 2.5, str(v), 6.8, "end", SUB, mono=True))
        b.append(f'<line x1="{left}" y1="{Y(-1):.1f}" x2="{right}" y2="{Y(-1):.1f}" stroke="{BAD}" stroke-width="0.6" stroke-dasharray="3 2"/>')
        ev_s, ev_e = int(rows[0]["event_start"]), int(rows[0]["event_end"])
        if ev_e > ev_s:
            b.append(f'<rect x="{X(ev_s):.1f}" y="{y0}" width="{X(ev_e) - X(ev_s):.1f}" height="{y1 - y0}" fill="{BAD}" opacity="0.08"/>')
        gs, ge = int(rows[0]["gene_start"]), int(rows[0]["gene_end"])
        b.append(f'<line x1="{X(gs):.1f}" y1="{y1 + 6}" x2="{X(ge):.1f}" y2="{y1 + 6}" stroke="{INK}" stroke-width="3"/>')
        b.append(text(X((gs + ge) / 2), y1 + 17, rows[0]["gene"], 7, "middle", italic=True))
        for r in rows:
            v = float(r["log2"])
            inside = ev_s <= (int(r["start"]) + int(r["end"])) / 2 <= ev_e
            b.append(f'<circle cx="{X((int(r["start"]) + int(r["end"])) / 2):.1f}" cy="{Y(v):.1f}" r="1.1" fill="{BAD if inside else GREY}"/>')
        b.append(text(right, y1 + 17, f"{rows[0]['chrom']}:{lo / 1e6:.2f}–{hi / 1e6:.2f} Mb", 6.8, "end", SUB, mono=True))
    b.append(text(10, 12 + ph * len(panels) / 2, "log2 depth ratio (1 kb bins)", 7.2, "middle", SUB, rotate=-90))
    return svg(w, h, b)


# ------------------------------------------------------------------ figure 6: read ends at the Alu sites
def fig_alu(d):
    rows = [r for r in d.sv["rows"] if r["type"] == "ALU" and r.get("readsAtSite") and r["readsAtSite"]["reads"]]
    groups = [("0|0", "Panel: no insertion (0|0)"), ("het", "Panel: heterozygous"), ("hom", "Panel: homozygous")]

    def group(r):
        g = r["panelGenotype"]
        return "0|0" if g == "0|0" else "hom" if g == "1|1" else "het"

    w, h = 720, 196
    left, right, top, bottom = 58, 560, 16, 150
    b = []

    def Y(v):
        return bottom - v * (bottom - top)

    for v in (0, 0.25, 0.5, 0.75, 1.0):
        b.append(f'<line x1="{left}" y1="{Y(v):.1f}" x2="{right}" y2="{Y(v):.1f}" stroke="{LINE}" stroke-width="0.6"/>')
        b.append(text(left - 6, Y(v) + 2.5, f"{v:.2f}", 6.8, "end", SUB, mono=True))
    b.append(text(14, (top + bottom) / 2, "fraction of reads clipped at the site", 7.2, "middle", SUB, rotate=-90))
    step = (right - left) / len(groups)
    for i, (key, label) in enumerate(groups):
        members = sorted((r for r in rows if group(r) == key), key=lambda r: (r["site"], r["sample"]))
        cx = left + step * (i + 0.5)
        b.append(text(cx, bottom + 14, label, 7.4, "middle"))
        found = sum(1 for r in members if r["readEndSignal"])
        b.append(text(cx, bottom + 25, f"{len(members)} evaluations, signature in {found}", 6.8, "middle", SUB, mono=True))
        n = len(members)
        for j, r in enumerate(members):
            x = cx + ((j + 0.5) / n - 0.5) * step * 0.72
            y = Y(r["readsAtSite"]["clippedAtSite"] / r["readsAtSite"]["reads"])
            rare = r["site"] == "HGSV_181052"
            fill = BLUE if r["readEndSignal"] else "#FFFFFF"
            stroke = BLUE if r["readEndSignal"] else (BAD if key != "0|0" else GREY)
            if rare:
                b.append(f'<rect x="{x - 2.6:.1f}" y="{y - 2.6:.1f}" width="5.2" height="5.2" fill="{fill}" stroke="{stroke}" stroke-width="1"/>')
            else:
                b.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.8" fill="{fill}" stroke="{stroke}" stroke-width="1"/>')
            if key != "0|0" and not r["readEndSignal"]:
                b.append(text(x + 5, y + 2.5, r["sample"], 6.4, "start", BAD, mono=True))
    lx, ly = right + 14, top + 40
    for k, (shape, fill, stroke, lab) in enumerate((("c", BLUE, BLUE, "signature reported"), ("c", "#FFFFFF", GREY, "no signature, panel 0|0"),
                                                    ("c", "#FFFFFF", BAD, "no signature, panel carrier"),
                                                    ("c", "#FFFFFF", SUB, "circle: chr13:32,348,302"), ("r", "#FFFFFF", SUB, "square: chr13:32,392,272"))):
        y = ly + 12 * k
        if shape == "c":
            b.append(f'<circle cx="{lx + 3}" cy="{y}" r="2.8" fill="{fill}" stroke="{stroke}" stroke-width="1"/>')
        else:
            b.append(f'<rect x="{lx + 0.4}" y="{y - 2.6}" width="5.2" height="5.2" fill="{fill}" stroke="{stroke}" stroke-width="1"/>')
        b.append(text(lx + 11, y + 2.5, lab, 6.8, "start", SUB))
    return svg(w, h, b)


# ------------------------------------------------------------------ html
def fonts_css(directory):
    if not directory:
        return ""
    css = []
    wanted = [("IBM Plex Mono", "ibm-plex-mono", ("400-normal", "500-normal", "400-italic")),
              ("Instrument Sans", "instrument-sans", ("400-normal", "500-normal", "400-italic", "500-italic"))]
    for family, slug, variants in wanted:
        for variant in variants:
            weight, style = variant.split("-")
            for subset in ("latin", "latin-ext"):
                hits = glob.glob(os.path.join(directory, "**", f"{slug}-{subset}-{variant}.woff2"), recursive=True)
                if not hits:
                    continue
                data = base64.b64encode(open(hits[0], "rb").read()).decode()
                rng = "U+0000-00FF,U+2013-2014,U+2018-201D,U+2022,U+2026,U+2032,U+2212,U+2264-2265" if subset == "latin" else "U+0100-024F,U+1E00-1EFF"
                css.append(f"@font-face{{font-family:'{family}';font-weight:{weight};font-style:{style};"
                           f"src:url(data:font/woff2;base64,{data}) format('woff2');unicode-range:{rng};}}")
    return "\n".join(css)


def pct(k, n):
    return f"{100 * k / n:.1f} %"


def build(d, fonts, version, date, doi):
    c = d.carriers
    snv_k, snv_n, _ = d.giab_total("gene_loci", "SNV")
    ind_k, ind_n, _ = d.giab_total("gene_loci", "INDEL")
    ex_k = sum(d.giab_total("exons_pm20", t)[0] for t in ("SNV", "INDEL"))
    ex_n = sum(d.giab_total("exons_pm20", t)[1] for t in ("SNV", "INDEL"))
    small_found = sum(v["found"] for v in d.spike["smallByClass"].values())
    small_total = sum(v["spiked"] for v in d.spike["smallByClass"].values())
    n_large = len(d.large())
    manta = d.count_large("manta")
    cnv_after = d.count_large("cnvkit")
    cnv_before = d.count_large("cnvkit", "before")
    hints = sum(1 for e in d.large() if e["exonDosageHints"])
    cnv_found = sorted((e for e in d.large() if e["cnvkit"]), key=lambda e: e["size"])
    cnv_min = fmt_bp(cnv_found[0]["size"]) if cnv_found else "n/a"
    del_big, del_small = d.sv_site("HGSV_216998"), d.sv_site("HGSV_216996")
    alu_c, alu_r = d.sv_site("HGSV_181048"), d.sv_site("HGSV_181052")
    alu_k, alu_n = alu_c["carriersReported"] + alu_r["carriersReported"], alu_c["carriers"] + alu_r["carriers"]
    alu_sig = alu_c["carriersWithReadEndSignal"] + alu_r["carriersWithReadEndSignal"]
    alu_ctrl = alu_c["controls"] + alu_r["controls"]
    alu_ctrl_sig = alu_c["controlsWithReadEndSignal"] + alu_r["controlsWithReadEndSignal"]
    controls = max(s["controls"] for s in d.sv["bySite"].values())
    control_hits = sum(s["controlsReported"] for s in d.sv["bySite"].values())
    cs_k, cs_n, cs_fp = d.concordance("exons_pm20", "SNV")
    ci_k, ci_n, ci_fp = d.concordance("exons_pm20", "INDEL")
    unexpected = [u for u in d.spike["unexpected"] if not (u["caller"] == "Manta" and u["type"] == "BND")]
    fp_exon = "no"
    tg = {n: {k: d.target_row(n, k) for k in ("SNV", "INDEL")} for n, _ in TARGETED}
    tg_snv_k = sum(v["SNV"][0] for v in tg.values())
    tg_snv_n = sum(v["SNV"][1] for v in tg.values())
    ts = {n: d.targeted[n]["qc"]["readQC"]["targetSummary"] for n, _ in TARGETED}
    tg_rows = "".join(
        f"<tr><td>{esc(label)}</td><td>{ts[n]['meanDepth']:.0f}×</td><td>{ts[n]['uniformityPct']:.1f} %</td><td>{ts[n]['fold80BasePenalty']:.2f}</td>"
        f"<td>{tg[n]['SNV'][0]} / {tg[n]['SNV'][1]}</td><td>{tg[n]['SNV'][2]}</td><td>{tg[n]['INDEL'][0]} / {tg[n]['INDEL'][1]}</td><td>{tg[n]['INDEL'][2]}</td>"
        f"<td>{d.targeted[n]['qc']['mosaicCandidates']}</td></tr>" for n, label in TARGETED)
    mos = [m for m in d.mosaic if m["spikedEvents"]]
    mos_neg = next(m for m in d.mosaic if not m["spikedEvents"])
    mos_fp = sum(m["falsePositiveMosaicCandidates"] + m["falsePositiveGermlineCalls"] for m in d.mosaic)
    mos_rows = "".join(
        f"<tr><td>{m['spikedAlleleFraction'] * 100:.0f} %</td><td>{m['reportedAsMosaicCandidate']}</td><td>{m['calledByGermlineCaller']}</td>"
        f"<td>{m['reportedEither']} / {m['spikedEvents']}</td><td>{m['falsePositiveMosaicCandidates'] + m['falsePositiveGermlineCalls']}</td></tr>" for m in mos)
    mos_n = mos[0]["spikedEvents"]
    mos_at = {round(m["spikedAlleleFraction"] * 100): m for m in mos}
    h_es = d.hprc_small("exons_pm20", "SNV")
    h_ei = d.hprc_small("exons_pm20", "INDEL")
    h_gs = d.hprc_small("gene_loci", "SNV")
    h_gi = d.hprc_small("gene_loci", "INDEL")
    h_n = len(d.hprc["samples"])
    tc_h, tc_hr = d.tandem(d.tandem_hprc, "carrier"), d.tandem(d.tandem_hprc, "reference")
    tc_o, tc_or = d.tandem(d.tandem_ont, "carrier"), d.tandem(d.tandem_ont, "reference")
    h_alu, h_vntr, h_tr = (d.hprc["bySite"][s] for s in ("chr13:32348295", "chr17:43100595", "chr17:43128264"))
    hprc_rows = "".join(
        f"<tr><td>{reg}</td><td>{kind}</td><td>{v[1]:,}</td><td>{v[0]:,}</td><td>{v[2]:,}</td><td>{v[3]:,}</td><td>{100 * v[0] / v[1]:.1f} %</td></tr>"
        for reg, kind, v in (("Exons ±20 bp", "SNV", h_es), ("Exons ±20 bp", "Indel", h_ei), ("Gene loci", "SNV", h_gs), ("Gene loci", "Indel", h_gi)))
    ont_pairs = len(d.ont)
    ont_agree = sum((r["longReads"] == "present") == (r["pipeline"] != "-") for r in d.ont)
    ont_samples = len({r["sample"] for r in d.ont})
    ont_panel_wrong = [r for r in d.ont if (r["longReads"] == "present") != (r["panel"] == "present")]
    carriers_rows = "".join(
        f"<tr><td>{esc(r['sample'])}</td><td>{esc(r['population'])}</td><td><i>{esc(r['gene'])}</i></td><td><code>{esc(r['hgvs'])}</code></td>"
        f"<td>{esc(r['consequence'])}</td><td>{'yes' if r['detected'] else 'no'}, {esc(r['genotype'] or '')}</td><td>{esc(r['depth'])}</td><td>{esc(r['AD'])}</td></tr>"
        for r in c["carriers"])
    sv_rows = ""
    for s in d.sv_sites:
        t = d.sv_site(s["id"])
        kind = "Alu insertion" if s["type"] == "ALU" else "Deletion"
        sv_rows += (f"<tr><td><i>{esc(s['gene'])}</i></td><td>{kind}, {int(s['length']):,} bp</td><td><code>{esc(s['chrom'])}:{int(s['pos']):,}</code></td>"
                    f"<td>{esc(s['callers_in_panel'])}</td><td>{t['carriersReported']} / {t['carriers']}</td>"
                    f"<td>{(str(t['carriersWithReadEndSignal']) + ' / ' + str(t['carriers'])) if s['type'] == 'ALU' else 'n/a'}</td>"
                    f"<td>{t['controlsReported'] + (t['controlsWithReadEndSignal'] if s['type'] == 'ALU' else 0)} / {t['controls']}</td></tr>")
    style = f"""
{fonts_css(fonts)}
@page {{ size: A4; margin: 17mm 16mm 18mm 16mm; }}
* {{ box-sizing: border-box; }}
body {{ font-family: {SANS}; font-size: 9.2pt; line-height: 1.42; color: {INK}; margin: 0; }}
.mast {{ display: flex; justify-content: space-between; align-items: center; border-bottom: 1.2px solid {INK}; padding-bottom: 6pt; margin-bottom: 12pt; }}
.mast .kind {{ font-family: {MONO}; font-size: 7.8pt; color: {SUB}; text-align: right; line-height: 1.35; }}
h1 {{ font-family: {MONO}; font-weight: 500; font-size: 16.5pt; line-height: 1.22; margin: 0 0 8pt; letter-spacing: -0.2pt; }}
h2 {{ font-family: {MONO}; font-weight: 500; font-size: 10.6pt; margin: 13pt 0 4pt; break-after: avoid; }}
h3 {{ font-family: {MONO}; font-weight: 500; font-size: 9.2pt; margin: 9pt 0 2pt; break-after: avoid; }}
p {{ margin: 0 0 5.5pt; text-align: justify; hyphens: auto; }}
.meta {{ font-size: 9.4pt; margin-bottom: 3pt; }} .meta strong {{ font-weight: 500; }}
.small {{ font-family: {MONO}; font-size: 7.6pt; color: {SUB}; line-height: 1.4; }}
.top {{ display: grid; grid-template-columns: 1fr 58mm; gap: 6mm; margin: 10pt 0 6pt; }}
.summary {{ border-top: 1px solid {LINE}; padding-top: 6pt; }}
.summary p {{ font-size: 9pt; }}
.highlights {{ background: #F3F1EA; border-radius: 5px; padding: 8pt 9pt; font-size: 8.3pt; line-height: 1.38; }}
.highlights h2 {{ margin: 0 0 4pt; font-size: 9pt; }} .highlights ul {{ margin: 0; padding-left: 11pt; }} .highlights li {{ margin-bottom: 3.5pt; }}
.ruo {{ border: 1px solid {WARN}; border-left-width: 3px; border-radius: 4px; padding: 5pt 8pt; font-size: 8.2pt; margin: 6pt 0 8pt; }}
.ruo strong {{ font-family: {MONO}; font-weight: 500; color: {WARN}; }}
.cols {{ column-count: 2; column-gap: 7mm; }}
figure {{ margin: 8pt 0 9pt; break-inside: avoid; column-span: all; }}
figure .frame {{ border: 1px solid {LINE}; border-radius: 5px; padding: 6pt; background: #FFFFFF; }}
figcaption, caption {{ font-size: 8pt; line-height: 1.35; color: {INK}; text-align: left; margin-top: 4pt; }}
figcaption strong, caption strong {{ font-family: {MONO}; font-weight: 500; }}
table {{ border-collapse: collapse; width: 100%; font-size: 7.9pt; margin: 4pt 0 9pt; break-inside: avoid; }}
.wide {{ column-span: all; }}
caption {{ caption-side: top; margin: 0 0 3pt; }}
th, td {{ border-bottom: 0.6px solid {LINE}; padding: 2.4pt 4pt; text-align: left; vertical-align: top; }}
th {{ font-family: {MONO}; font-weight: 500; font-size: 7.4pt; color: {SUB}; border-bottom: 1px solid {INK}; }}
code {{ font-family: {MONO}; font-size: 7.9pt; }}
pre {{ font-family: {MONO}; font-size: 7.4pt; background: #F3F1EA; padding: 5pt 6pt; border-radius: 4px; white-space: pre-wrap; }}
ul {{ margin: 0 0 5.5pt; padding-left: 12pt; }} li {{ margin-bottom: 2pt; }}
a.cite {{ color: inherit; text-decoration: none; }}
ol.refs {{ font-size: 7.6pt; padding-left: 14pt; line-height: 1.33; }} ol.refs li {{ margin-bottom: 1.5pt; }}
"""
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>brca12-wgs: a reproducible Nextflow pipeline for BRCA1/BRCA2 germline analysis from genome, exome and panel sequencing on GRCh38</title>
<style>{style}</style></head><body>

<div class="mast"><div>{LOGO}</div><div class="kind">Software report · preprint<br>Version {esc(version)} · {esc(date)}</div></div>

<h1>brca12-wgs: a reproducible Nextflow pipeline for <i>BRCA1</i>/<i>BRCA2</i> germline analysis from genome, exome and panel sequencing on GRCh38</h1>
<div class="meta"><strong>Peter Jungwirth, MSc</strong> · PJ Labs, Austria</div>
<div class="small">Software: https://github.com/pjlabsdev/brca12-wgs · DOI: <span id="doi">{esc(doi)}</span> · Licence: MIT</div>

<div class="top">
<div class="summary">
<p><strong>Summary.</strong> Germline pathogenic variants in <i>BRCA1</i> and <i>BRCA2</i> cause hereditary breast and ovarian cancer syndrome. The relevant spectrum spans single-nucleotide variants, small insertions and deletions, exon-level copy-number changes and larger rearrangements, yet many open BRCA workflows target amplicon panels on GRCh37 and report small variants only. <i>brca12-wgs</i> is a Nextflow pipeline that takes FASTQ files or a BAM/CRAM from whole-genome, exome, hybrid-capture panel or amplicon sequencing and produces a structured BRCA1/BRCA2 report on GRCh38. It combines BWA-MEM, GATK HaplotypeCaller or DeepVariant, Manta, CNVkit with per-exon depth ratios, Ensembl VEP with MANE Select transcripts, ClinVar and mosdepth, and lists automated hints on single ACMG/AMP criteria following the ClinGen ENIGMA specification without assigning a classification. We measured performance on four kinds of whole-genome data: Genome in a Bottle truth sets ({snv_k} of {snv_n} SNVs and {ind_k} of {ind_n} indels), ten 1000 Genomes carriers of pathogenic variants ({c['detected']} of {c['total']} detected), {small_total} small and {n_large} large variants spiked into real reads, and real Alu insertions and deletions in 1000 Genomes samples, checked in long reads. On {h_n} samples with assembly-based truth from the Human Pangenome Reference Consortium it found {h_es[0]} of {h_es[1]} exonic SNVs without a false positive, and on two exomes and an amplicon panel {tg_snv_k} of {tg_snv_n} SNVs inside the targets. Spiked variants at 20 % allele fraction were all reported, with mosaic candidates from Mutect2 below the germline caller's range. Real Alu insertions were called in only {alu_k} of {alu_n} carriers, but an insertion signature from clipped read ends, added for this reason, was present in {alu_sig} of {alu_n}. Manta reported all {n_large} spiked large events, whereas single-sample CNVkit segmentation reported {cnv_after} of {n_large} after a correction made during this work ({cnv_before} before). The software is for research use; results require confirmation by a clinically validated method.</p>
<p><strong>Availability.</strong> Source code, test data, validation scripts and results are available under the MIT licence at https://github.com/pjlabsdev/brca12-wgs.</p>
</div>
<div class="highlights"><h2>Highlights</h2><ul>
<li>One command from WGS reads to a BRCA1/BRCA2 report with small variants, structural variants, copy number and exon coverage on GRCh38.</li>
<li>{c['detected']} of {c['total']} real carriers of pathogenic variants from the 1000 Genomes Project were detected with the correct genotype.</li>
<li>{small_found} of {small_total} spiked SNVs and indels up to 40 bp and {manta} of {n_large} spiked deletions, duplications and an Alu insertion were reported.</li>
<li>Real Alu insertions were called in {alu_k} of {alu_n} carriers; a read-end signature added for this reason was present in {alu_sig} of {alu_n} and in {alu_ctrl_sig} of {alu_ctrl} non-carrier evaluations.</li>
<li>Against assembly-based truth in {h_n} HPRC samples: {h_es[0]} of {h_es[1]} SNVs in the exons ±20 bp, no false positive.</li>
<li>Exome and amplicon data: {tg_snv_k} of {tg_snv_n} SNVs inside the targets; the misses lie in a UTR stretch without reads.</li>
<li>Mosaic candidates (Mutect2) and the germline caller together reported every spiked variant at 20 % allele fraction, with no false positive.</li>
<li>The validation exposed and fixed three defects, including missed copy-number segments.</li>
<li>Research use only. No classification is assigned.</li>
</ul></div>
</div>

<div class="ruo"><strong>Research use only.</strong> brca12-wgs is not a medical device and has not been clinically validated. Its output must not be used for diagnosis, risk assessment or treatment decisions without confirmation in an accredited laboratory (for example Sanger sequencing for small variants and MLPA for copy-number changes). ACMG/AMP entries are automated hints on individual criteria, not variant classifications. A result without findings does not exclude a hereditary BRCA1/BRCA2 predisposition.</div>

<figure><div class="frame">{fig_workflow()}</div>
<figcaption><strong>Figure 1. Workflow.</strong> Reads are checked against GRCh38 and aligned. Four evidence layers are computed independently and kept separate in the report. Hints on single ACMG/AMP criteria are attached to annotated small variants; no combination rule is applied.</figcaption></figure>

<div class="cols">
<h2>1 Introduction</h2>
<p>Pathogenic germline variants in the tumour suppressor genes <i>BRCA1</i> (chr17, minus strand) and <i>BRCA2</i> (chr13, plus strand) substantially increase lifetime risks of breast, ovarian, prostate and pancreatic cancer [@kuchenbaecker]. Most pathogenic alleles are truncating SNVs and small indels, but large genomic rearrangements such as single- or multi-exon deletions and duplications account for a relevant fraction of findings, particularly in <i>BRCA1</i>, which is rich in Alu repeats and has a partially duplicated pseudogene (<i>BRCA1P1</i>) near its 5′ end [@sluiter]. Insertions of mobile elements are rare but include a founder variant, the Alu insertion c.156_157insAlu in <i>BRCA2</i> exon 3 [@teugels, @machado].</p>
<p>Whole-genome sequencing (WGS) captures all of these variant classes in a single assay, but extracting a reliable gene-level result requires several specialised tools, consistent coordinates and explicit statements about where the data are not sufficient. Openly available workflows such as BRCA-analyzer [@brcaanalyzer] were developed for amplicon panels on GRCh37/hg19, depend on licence-restricted annotation and do not report copy-number or structural variants. General-purpose germline pipelines such as Sarek [@sarek] cover WGS but are not tailored to a concise, gene-focused report with exon-level coverage and conservative interpretation hints. Because many laboratories sequence exomes or gene panels rather than genomes, the same report should also be available for targeted data.</p>
<p><i>brca12-wgs</i> addresses this gap with three design goals: (i) GRCh38 only, with automatic rejection of other assemblies; (ii) separate, traceable evidence layers for small variants, copy-number candidates, structural variants and coverage; and (iii) reproducibility through pinned software environments, recorded tool versions and a verbatim command log.</p>

<h2>2 Implementation</h2>
<p>The pipeline is written in Nextflow DSL2 [@nextflow] and runs locally with conda environments, a Docker or Singularity image, or tools on the system path (Figure 1, Table 1).</p>

<table><caption><strong>Table 1.</strong> Processes, tools and outputs.</caption>
<tr><th>Step</th><th>Tool (pinned version)</th><th>Purpose</th></tr>
<tr><td>Input check</td><td>Python</td><td>GRCh38 detection from chr13/chr17 lengths; region and exon BED harmonisation</td></tr>
<tr><td>Alignment</td><td>BWA-MEM 0.7.18 [@bwa], samtools 1.21 [@samtools]</td><td>Alignment, fixmate, sort, duplicate marking (not for amplicons), primer clipping with <code>ampliconclip</code></td></tr>
<tr><td>Read QC</td><td>samtools stats, mosdepth</td><td>Mapping, pairing, duplicates, error rate; per target depth, uniformity, fold-80, on-target reads</td></tr>
<tr><td>SNVs/indels</td><td>GATK 4.6.1.0 HaplotypeCaller [@gatk] or DeepVariant [@deepvariant]</td><td>Calling in both genes ± 200 bp</td></tr>
<tr><td>Mosaic candidates</td><td>GATK Mutect2 [@mutect2] (optional)</td><td>Calls at 2–30 % allele fraction not made by the germline caller</td></tr>
<tr><td>Structural variants</td><td>Manta 1.6.0 [@manta]</td><td>Breakpoints within ± 1 Mb of both genes</td></tr>
<tr><td>Copy number</td><td>CNVkit 0.9.11 [@cnvkit]</td><td>WGS, hybrid or amplicon mode, flat or pooled reference; per-exon depth ratios</td></tr>
<tr><td>Annotation</td><td>Ensembl VEP 113 [@vep], bcftools</td><td>MANE Select [@mane], HGVS, gnomAD [@gnomad], ClinVar [@clinvar]</td></tr>
<tr><td>Coverage</td><td>mosdepth 0.3.10 [@mosdepth]</td><td>Per exon at MAPQ ≥ 20 and ≥ 0; 10/20/30×</td></tr>
<tr><td>Report</td><td>Python (standard library)</td><td>HTML, Markdown, PDF, Excel (English, German or Russian), JSON</td></tr>
</table>

<h3>2.1 Reference and regions</h3>
<p>The input check reads the FASTA index and accepts the reference only if chr13 and chr17 have the GRCh38 lengths; GRCh37 is rejected with an explicit message. Gene loci are fixed to chr17:43,044,295–43,170,327 (<i>BRCA1</i>) and chr13:32,315,077–32,400,268 (<i>BRCA2</i>). Exons are taken from the MANE Select transcripts [@mane] NM_007294.4 (23 exons) and NM_000059.4 (27 exons). For <i>BRCA1</i> the report also shows the legacy exon numbering, which omits exon 4.</p>

<h3>2.2 Variant calling</h3>
<p>Reads are aligned with BWA-MEM [@bwa]; sorting, mate fixing and duplicate marking use samtools [@samtools]. Small variants are called with GATK HaplotypeCaller [@gatk] restricted to the padded gene regions, or optionally with DeepVariant [@deepvariant]. Manta [@manta] is run with call regions extending 1 Mb around both genes, so that rearrangements with one breakpoint inside a gene are retained, including breakend (BND) pairs. Two breakends a few bases apart whose mates lie elsewhere are flagged as a possible insertion of foreign sequence such as a mobile element. In addition, soft-clipped read ends inside the genes are searched for the signature of an insertion: reads from the left clipped at one position, reads from the right clipped at a position up to 30 bp before it (the target-site duplication), at least three clipped reads on each side and a clipped fraction of at least 15 %. Such a pair is listed as an insertion signature when the clipped sequence contains k-mers of the Alu consensus; it is a hint from the reads, not a variant call. A second signature marks an extra copy of a tandem repeat unit: reads clipped at the end of a unit whose clipped part matches the start of that unit in the reference, and reads clipped at its start whose clipped part matches its end, both pointing to the same unit length of 50 to 1,000 bp, with at least two reads at each border. For copy number, CNVkit [@cnvkit] is run in WGS mode with 1 kb bins over chr13 and chr17 and a flat reference; bins with almost no coverage are excluded from segmentation (Section 3.3). Segments overlapping a gene with an estimated copy number different from two are reported as candidates with the affected exons. Independently, the mean depth of each exon is divided by the median depth of 100 kb windows on the same chromosome; ratios below 0.65 or above 1.35 are listed as dosage hints. Because only a single sample is available, both are labelled as candidates that require MLPA [@mlpa] or another validated method.</p>

<h3>2.3 Annotation and ACMG/AMP hints</h3>
<p>VEP [@vep] annotates consequences, HGVS notation and exon or intron number on the MANE Select transcript. ClinVar [@clinvar] significance, review status and condition come from a user-supplied ClinVar VCF, so that the ClinVar version is under the user's control. The pipeline then evaluates single ACMG/AMP criteria [@acmg] following version 1.2 of the ClinGen ENIGMA BRCA1/BRCA2 specification [@enigma]: BA1, BS1 and BS1_Supporting from gnomAD [@gnomad] non-cancer filtering allele frequencies in non-founder populations, PM2_Supporting for substitutions absent from both gnomAD data sets, PP3 and BP4 from BayesDel [@bayesdel] and SpliceAI [@spliceai] inside the defined functional domains, BP1_Strong outside them and BP7 for silent and deep intronic variants. The gnomAD, SpliceAI and BayesDel files are supplied by the user; without them the corresponding hints are reported as open or as approximations. Truncating and canonical splice variants are reported as PVS1 candidates, restricted in the last exon [@pvs1]; the exon-specific PVS1 decision tree is not yet applied. No combination rule is used and no classification is assigned; the JSON field <code>clinicalClassification</code> is always null.</p>

<h3>2.4 Coverage and uncertain regions</h3>
<p>Depth per exon is computed with mosdepth [@mosdepth]. An exon is flagged as uncertain if its mean depth is below 20×, if fewer than 95 % of its bases reach 20×, or if more than 10 % of reads have mapping quality below 20, which identifies regions affected by pseudogene homology or repeats. All intervals below 10× are listed. Small variants receive quality flags for non-PASS filters, low depth, low genotype quality, heterozygous allele fractions outside 25–75 %, low mapping quality and location in an uncertain exon. All thresholds are parameters.</p>

<h3>2.5 Reporting and reproducibility</h3>
<p>The report has separate sections for <i>BRCA1</i> and <i>BRCA2</i> small variants, copy-number candidates, structural variants, coverage per exon and uncertain intervals, followed by tool versions, the invocation and limitations. It is produced as HTML and Markdown in English or German and as JSON (schema <code>genome.brca12.v1</code>). Tool versions are captured at run time, every task copies its command script into <code>brca12_commands.sh</code>, and conda environments pin all tool versions. Unit tests and a stub run of the complete workflow run in GitHub Actions on every change.</p>

<h3>2.6 Targeted sequencing and low allele fractions</h3>
<p>With <code>--assay exome</code>, <code>panel</code> or <code>amplicon</code> and a target BED, the target intervals inside both genes are kept and exons without a target are reported as not assessable. Amplicon reads are not marked as duplicates; with a primer BED, primers are removed with <code>samtools ampliconclip</code> and the share of reads with a clipped primer is reported. CNVkit runs in hybrid or amplicon mode, optionally with a pooled normal reference, and the exon dosage baseline becomes the median exon depth of the sample. Single-end reads are accepted (Manta is skipped), and <code>--platform iontorrent</code> adds a warning to every indel. A read QC step reports mapping, pairing, duplicate and error rates and, for targets, mean depth, uniformity, fold-80 base penalty, on-target reads and weak targets. <code>--mosaic</code> adds Mutect2 candidates at low allele fraction (Section 3.6).</p>

<h2>3 Validation</h2>
<p>Eight data sets were used (Figures 2 and 3). All scripts, sample lists and result files are in the <code>validation/</code> directory of the repository.</p>

<h3>3.1 Genome in a Bottle truth sets</h3>
<p>Small-variant calls were compared with the GIAB NIST v4.2.1 GRCh38 truth sets for HG001 and HG002 [@giab1, @giab2]. Illumina HiSeq 2×148 bp read pairs that the GIAB alignments placed within either gene ±10 kb were extracted, downsampled to about 36× and 30×, converted back to FASTQ and processed from FASTQ. Calls were matched with RTG Tools <code>vcfeval</code> [@vcfeval], which compares haplotypes and requires the genotype to agree, inside the GIAB high-confidence regions (91.9 % and 95.1 % of the 211,225 bp of both loci). The pipeline recovered {snv_k} of {snv_n} SNVs and {ind_k} of {ind_n} indels, and all {ex_n} variants in the exons ±20 bp without a false positive. The three discordant calls were intronic; two were genotype errors in homopolymers at 12× and 14× that carried the pipeline's own quality flags.</p>

<h3>3.2 Real carriers of pathogenic variants</h3>
<p>The 30× joint call set of the 3,202 samples of the 1000 Genomes Project [@kgp] was searched for loss-of-function variants in the MANE coding sequence and its canonical splice sites. Ten individuals carry eight such variants, seven of them annotated as pathogenic in the ClinVar release used by that call set (Table 2). Read pairs over both genes ±10 kb were taken from the public CRAM files and processed from FASTQ. All {c['detected']} variants were called heterozygous with genotype quality 99, including a 5 bp deletion that the 1000 Genomes call set itself marks as low confidence. As a secondary check, calls in the exons ±20 bp agreed with the 1000 Genomes genotypes of the same individuals for {cs_k} of {cs_n} SNVs and {ci_k} of {ci_n} indels, with {cs_fp} and {ci_fp} additional calls; that call set is a comparison, not a truth set.</p>

<table class="wide"><caption><strong>Table 2.</strong> Carriers of pathogenic <i>BRCA1</i>/<i>BRCA2</i> variants in the 1000 Genomes Project and the result of the pipeline. AD, reads supporting reference and variant allele.</caption>
<tr><th>Sample</th><th>Pop.</th><th>Gene</th><th>Variant (MANE)</th><th>Consequence</th><th>Detected, genotype</th><th>Depth</th><th>AD</th></tr>
{carriers_rows}
</table>

</div>

<figure><div class="frame">{fig_gene_maps(d)}</div>
<figcaption><strong>Figure 2. Variants used for validation, drawn on both gene loci.</strong> Black boxes are MANE Select exons. Above the gene: pathogenic variants of real carriers. Below: spiked small variants (one tick each), spiked large events and the positions of real structural variants.</figcaption></figure>

<figure><div class="frame">{fig_recall(d)}</div>
<figcaption><strong>Figure 3. Detection by data set and variant class.</strong> Points are the fraction of known variants that the pipeline reported; bars are 95 % Wilson intervals [@wilson]. Green, all found; amber and red, partly found. Small numbers give wide intervals.</figcaption></figure>

<div class="cols">
<h3>3.3 Spike-in experiment</h3>
<p>No openly available short-read WGS data set of a carrier of a pathogenic exon deletion, duplication or Alu insertion was found; the GeT-RM reference materials [@getrm] are distributed as DNA. We therefore spiked synthetic heterozygous variants into real reads. The background was GIAB HG002 on chr13 and chr17 at about 30×, which gives CNVkit a real chromosome-wide baseline. Inside both genes ±20 kb every real read pair was replaced by a pair with the same position, fragment length, read lengths and base qualities, read from one of two haplotypes: the reference, or the reference with the spiked events. Sequencing errors were drawn from the real base qualities. Coverage profile and fragment sizes therefore stayed those of the real run. The replaced pairs were aligned with BWA-MEM and merged with the untouched background, and the pipeline was run with real Manta and CNVkit.</p>
<p>Two samples carried {small_total} small variants: one per coding exon segment with rotating types, SNVs at canonical splice sites, and five variants modelled on well-known pathogenic alleles (<i>BRCA1</i> c.68_69del, c.181T&gt;G, c.5266dup; <i>BRCA2</i> c.771_775del, c.5946del), whose coordinates were derived from a coding-sequence model that is checked by translation. All {small_found} were called with the correct genotype, with {fp_exon} false positive in the exons (Figure 3). Outside the exons, every sample including the negative control showed the same two calls at one intronic <i>BRCA2</i> site, which come from real HG002 reads with an unmapped mate that the replacement does not cover, and the samples with large deletions showed short artefactual calls directly at some breakpoints.</p>
<p>Six samples carried {n_large} large events from 120 bp to a whole-gene deletion, two tandem duplications and an Alu insertion in <i>BRCA2</i> exon 3 (Figure 4). Manta reported {manta} of {n_large} with exact breakpoints and heterozygous genotype; the Alu insertion appeared as a pair of breakends. The read-end signature (Section 2.2) listed the spiked Alu insertion and nothing in the other eight spike-in samples. With the settings of version 1.0.2, CNVkit reported {cnv_before} of {n_large}: bins without coverage dominated the segmentation statistics, and even whole-gene deletions with a clear half-depth signal were merged into their neighbours (Figure 5). Excluding such bins (<code>--drop-low-coverage</code>) raised this to {cnv_after} of {n_large}; the smallest event reported by CNVkit was {cnv_min}. In addition, exon dosage hints had been suppressed for exactly the affected exons, because a heterozygous deletion also makes an exon fall below the depth threshold. After the correction, {hints} of the {n_large - 1} deletions and duplications produced a dosage hint for an affected exon. A third defect stopped the pipeline when the input BAM was named like the sample.</p>
</div>

<figure><div class="frame">{fig_large(d)}</div>
<figcaption><strong>Figure 4. Spiked large events and the evidence layer that reported them.</strong> CNVkit 1.0.2 and 1.1.0 differ only in the exclusion of bins without coverage. Single-sample depth methods are not expected to resolve events of one or a few 1 kb bins; n/a marks layers that cannot apply.</figcaption></figure>

<figure><div class="frame">{fig_depth(d)}</div>
<figcaption><strong>Figure 5. Read-depth signal of spiked heterozygous deletions.</strong> Each point is a 1 kb bin (CNVkit log2 ratio). The dashed line marks the expected value for one copy (−1). Shaded: spiked interval. The signal is present in the bins; whether it is reported depends on the segmentation.</figcaption></figure>

<div class="cols">
<h3>3.4 Real structural variants</h3>
<p>The phased SNV, indel and SV panel of the 1000 Genomes Project [@kgp] contains two Alu insertions in <i>BRCA2</i> and two deletions in <i>BRCA1</i>, all intronic or otherwise not disease-causing. They allow the same kind of event to be tested in real reads (Table 3). Carriers and six samples that are homozygous reference at all four sites were processed as in Section 3.2, with one addition: mates that the original alignment had placed elsewhere in the genome were retrieved as well, because read pairs with a mate inside the inserted element are otherwise lost: without them only {{SV_NOMATES}}. At the two Alu sites the panel genotypes were looked up for all 33 samples of Sections 3.2 and 3.4, which added three heterozygous carriers of the common insertion; the ten samples of Section 3.2 were rerun with retrieved mates for this purpose. The panel calls are statistically phased short-read calls (the Alu insertions come from MELT [@melt] within the GATK-SV pipeline [@gatksv]) and serve as a comparison, not as truth.</p>

</div>

<table><caption><strong>Table 3.</strong> Real structural variants in 1000 Genomes samples. Called: by Manta or the small-variant caller at the site. Signature: insertion signature from clipped read ends. Controls: samples with panel genotype 0|0 and any call or signature at the site.</caption>
<tr><th>Gene</th><th>Event</th><th>Position</th><th>Panel caller</th><th>Called</th><th>Signature</th><th>Controls</th></tr>
{sv_rows}
</table>

<figure><div class="frame">{fig_alu(d)}</div>
<figcaption><strong>Figure 6. Read ends at the two real Alu insertion sites.</strong> Each symbol is one sample at one site: the fraction of reads over the site that are clipped there by at least 10 bases, grouped by the genotype in the 1000 Genomes panel. Filled: the pipeline listed an insertion signature with Alu sequence. The raw fraction counts every read within 15 bp of the site and is therefore lower than the clipped fraction computed by the detector.</figcaption></figure>

<div class="cols">
<p>{{SV_TEXT}}</p>

<h3>3.5 Exome and amplicon data</h3>
<p>Three public GIAB data sets from targeted sequencing [@giabdata] were run from FASTQ with <code>--assay exome</code> or <code>--assay amplicon</code> (Table 4): an Agilent SureSelect v5 exome of HG002 (HiSeq 2×125 bp), an Illumina TruSeq exome of HG001 and an Ion AmpliSeq exome of HG002 (Ion Proton, single-end, run with <code>--platform iontorrent</code>). Reads over both genes ±10 kb were taken from the producers' hg19 alignments; vendor target intervals were moved to GRCh38 with one offset per locus, and for the Agilent kit, whose BED is not public, the MANE exons ±50 bp served as targets. Truth was GIAB v4.2.1 inside the targets and the high-confidence regions. Inside the targets {tg_snv_k} of {tg_snv_n} SNVs were found. All missed SNVs and indels of the two exomes lie in the 3′ UTR of <i>BRCA2</i> exon 27, where the source alignments have 0 to 2 reads; the report lists this stretch as an uncertain interval. The three SNV false positives of the Ion data had allele fractions of 10 to 20 % and carried the allele-balance flag; its 98 false indels are homopolymer errors of the platform, and with <code>--platform iontorrent</code> every indel carries a warning. Target QC reported the expected spread: uniformity from 83.9 to 91.6 % and the 1.27 % error rate of the Ion run as a flag.</p>
</div>

<table class="wide"><caption><strong>Table 4.</strong> Exome and amplicon data against GIAB v4.2.1, inside the targets. Uniformity: share of target bases with at least 0.2 × the mean depth. Fold-80: mean depth divided by the 20th percentile. Mosaic: candidates reported with <code>--mosaic</code>.</caption>
<tr><th>Data set</th><th>Mean depth</th><th>Uniformity</th><th>Fold-80</th><th>SNVs</th><th>SNV FP</th><th>Indels</th><th>Indel FP</th><th>Mosaic</th></tr>
{tg_rows}
</table>

<div class="cols">
<h3>3.6 Low allele fractions</h3>
<p>With <code>--mosaic</code> the pipeline runs GATK Mutect2 [@mutect2] in tumour-only mode with FilterMutectCalls over the same regions and reports PASS calls with an allele fraction of 2 to 30 %, at least three supporting reads and at least 20× depth that the germline caller did not report, as mosaic candidates. The {mos_n} small variants of Section 3.3 were spiked into the real HG002 reads at about 30× with 5, 10, 15, 20 and 30 % of the read pairs drawn from the variant haplotype (Table 5, Figure 7). At 20 % and above every variant was reported, by the germline caller or as a mosaic candidate; at 10 % {mos_at[10]['reportedEither']} of {mos_n} and at 5 % {mos_at[5]['reportedEither']} of {mos_n}. At 30× a 5 % variant is carried by about 1.5 reads on average, so low fractions need the depth of panels or exomes. The negative control and all spiked samples gave {mos_fp} false positives. The candidate filters were set after a first run on the exome and amplicon data of Section 3.5, where they removed four TruSeq candidates at 8× with two supporting reads each and 23 Ion indels; that is therefore not an independent test.</p>
</div>

<table><caption><strong>Table 5.</strong> Spiked variants at low allele fractions in real HG002 reads (about 30×). Mosaic: reported as mosaic candidate; germline: called by HaplotypeCaller (with the allele-balance flag below 25 %).</caption>
<tr><th>Fraction</th><th>Mosaic</th><th>Germline</th><th>Either</th><th>False positives</th></tr>
{mos_rows}
</table>

<figure><div class="frame">{fig_mosaic(d)}</div>
<figcaption><strong>Figure 7. Spiked variants by allele fraction.</strong> Each bar is one sample with {mos_n} spiked small variants: reported as mosaic candidate (blue), called by the germline caller (green) or not reported (grey).</figcaption></figure>

<div class="cols">
<h3>3.7 Assembly-based truth</h3>
<p>The {h_n} individuals of the Human Pangenome Reference Consortium [@hprc] that are also part of the 1000 Genomes 30× data [@kgp] provide a truth that does not depend on short reads: the haplotype-resolved long-read assemblies, here as the GRCh38-based VCF of the Minigraph-Cactus pangenome [@mgcactus]. It covers both genes completely, introns and repeats included. Read pairs were prepared as in Section 3.4 and compared per sample with <code>vcfeval</code>; positions where an assembly has no allele were excluded. In the exons ±20 bp {h_es[0]} of {h_es[1]} SNVs were found without a false positive (Table 6). {h_ei[3] - 2} of the {h_ei[3]} missed and 9 of the {h_ei[2]} false exon indels lie in one poly-T stretch of the <i>BRCA1</i> 3′ UTR (chr17:43,044,804–43,044,817), where assemblies and short reads disagree on the length of the homopolymer; outside it the exon indels were 41 of 43 with one false positive. Over the whole gene loci the figures include every intronic repeat and were not curated.</p>
<p>Variants of 50 bp or more formed three sites. The common intronic Alu insertion in <i>BRCA2</i> showed the read-end signature in {h_alu['carriersWithReadEndSignal']} of {h_alu['carriers']} carriers and in {h_alu['referenceWithReadEndSignal']} of {h_alu['reference']} non-carriers, while Manta called it in {h_alu['carriersCalled']}. A VNTR in <i>BRCA1</i> intron 6 with alleles from −80 to +113 bp was called in {h_vntr['carriersCalled']} of {h_vntr['carriers']} carriers, with {h_vntr['referenceCalled']} calls in {h_vntr['reference']} non-carriers. A 126 bp insertion that adds one copy of an adjacent tandem repeat was called in {h_tr['carriersCalled']} of {h_tr['carriers']} carriers: the extra copy aligns to the existing one with soft-clipped ends only. The tandem-copy signature (Section 2.2) was written after this result and found the insertion in {tc_h[0]} of {tc_h[1]} carriers and {tc_hr[0]} of {tc_hr[1]} non-carriers; its only other signature in the 39 samples was a 73 bp copy at the VNTR in a carrier of 68 and 76 bp alleles. As an independent check it ran on the samples of Sections 3.2 and 3.4 that have long reads (Section 3.8): {tc_o[0]} of {tc_o[1]} long-read carriers and {tc_or[0]} of {tc_or[1]} non-carriers had the signature. Carriers without it had fewer than two supporting clipped reads at one or both borders.</p>
</div>

<table><caption><strong>Table 6.</strong> {h_n} HPRC samples against assembly-based truth, small variants below 50 bp.</caption>
<tr><th>Region</th><th>Type</th><th>Truth</th><th>TP</th><th>FP</th><th>FN</th><th>Recall</th></tr>
{hprc_rows}
</table>

<div class="cols">
<h3>3.8 Long reads</h3>
<p>The structural variants of Section 3.4 come from short-read callsets. For the {ont_samples} of their samples that are in the 1KG-ONT cohort [@ont1kg], every site was genotyped in the Oxford Nanopore reads: a read spanning the site supports the variant when its alignment carries an insertion or deletion of at least half the expected length there. The pipeline result agreed with the long reads in {ont_agree} of {ont_pairs} sample-site pairs. {{ONT_PANEL}}</p>
{{ONT_ALL}}

<h3>3.9 Limits of this validation</h3>
<ul>
<li>Reads for Sections 3.1, 3.2, 3.4, 3.5 and 3.7 were taken from existing alignments near the genes and realigned against chr13 and chr17 only. This removes competition from paralogous sequence elsewhere, which may overstate precision. A whole-genome run against the complete reference is the next step.</li>
<li>Spiked variants are synthetic. They test whether a variant that is present in the reads is called and reported, not how real pathogenic rearrangements look. The spiked Alu element was copied from the reduced reference, which makes its detection easier than a real mobile element insertion.</li>
<li>No Ensembl VEP cache was available on the validation machine, so annotation and the ACMG/AMP hints were not part of these runs. In Sections 3.1, 3.2 and 3.4 CNVkit was replaced by a neutral stand-in.</li>
<li>Numbers are small. Ten carriers and twelve large events give wide confidence intervals (Figure 3). The mosaic filters were tuned on the data of Section 3.5.</li>
<li>The read-end signature was written after the missed Alu insertions had been inspected, so its result on these carriers is not an independent test.</li>
<li>No real sample with a pathogenic exon deletion, duplication or Alu insertion was tested.</li>
</ul>

<h2>4 Comparison with existing workflows</h2>
<table><caption><strong>Table 7.</strong> Feature comparison with BRCA-analyzer [@brcaanalyzer], as described in its repository.</caption>
<tr><th></th><th>brca12-wgs 1.1</th><th>BRCA-analyzer</th></tr>
<tr><td>Data type</td><td>WGS, exome, hybrid-capture and amplicon panels; paired or single-end; Illumina or Ion Torrent</td><td>Amplicon panels (MiSeq)</td></tr>
<tr><td>Input</td><td>FASTQ, BAM or CRAM</td><td>FASTQ</td></tr>
<tr><td>Assembly</td><td>GRCh38, enforced; GRCh37 position for every small variant</td><td>hg19</td></tr>
<tr><td>Primer removal</td><td>samtools ampliconclip, share of clipped reads</td><td>cutPrimers, share of trimmed reads</td></tr>
<tr><td>Small variants</td><td>GATK 4 HaplotypeCaller or DeepVariant</td><td>Pisces, GATK 3 UnifiedGenotyper or FreeBayes</td></tr>
<tr><td>Low allele fraction</td><td>Mutect2 mosaic candidates, tested by spike-ins</td><td>Pisces calls low fractions</td></tr>
<tr><td>Copy number / SV</td><td>CNVkit and exon dosage / Manta, Alu and tandem-copy read-end signatures</td><td>Not included</td></tr>
<tr><td>Annotation</td><td>VEP, MANE Select, gnomAD, user ClinVar; BIC names</td><td>ANNOVAR (licence), snpEff, bundled ClinVar and BIC</td></tr>
<tr><td>Interpretation support</td><td>ACMG/AMP criterion hints (ENIGMA v1.2), no classification</td><td>Class from ClinVar/BIC</td></tr>
<tr><td>QC</td><td>Read QC, per exon and per target coverage, uniformity, fold-80, intervals below 10×</td><td>Read counts, amplicon coverage and uniformity</td></tr>
<tr><td>Batches</td><td>Sample sheet, cohort table, contamination hint</td><td>Joint table, contamination hint</td></tr>
<tr><td>Output</td><td>HTML, Markdown, PDF, Excel, JSON; English, German, Russian</td><td>Excel, Word; English, Russian</td></tr>
<tr><td>Reproducibility</td><td>Nextflow, pinned conda, Docker/Singularity, command log, tests in CI</td><td>Python script, Docker image</td></tr>
<tr><td>Open validation</td><td>GIAB WGS, exomes and amplicons, real carriers, spike-ins, HPRC assemblies, long reads</td><td>Example data</td></tr>
<tr><td>Clinical validation</td><td>None (research use only)</td><td>More than 900 patient samples with Sanger confirmation, according to its documentation</td></tr>
<tr><td>Licence</td><td>MIT</td><td>GPL-3.0</td></tr>
</table>

<h2>5 Limitations</h2>
<ul>
<li>At 30× WGS, variants below about 15 % allele fraction are often missed even with <code>--mosaic</code> (Section 3.6); low fractions need deeper targeted data. Sensitivity is also limited for complex rearrangements and in regions with pseudogene homology such as <i>BRCA1</i> exons 1–2.</li>
<li>Extra copies of a tandem repeat unit are reported only as read-end hints, found in about two of three carriers of a common 126 bp insertion (Section 3.7); the rule was developed on the same samples and tested independently on only six carriers. Indel lengths in long homopolymers can differ from assembly-based truth.</li>
<li>Ion Torrent indels are dominated by homopolymer errors and need confirmation; exons without a target or without reads in targeted data are reported as not assessable, not as negative.</li>
<li>Copy-number results come from a single sample with a flat reference and are candidates only. In the spike-in experiment CNVkit alone missed every event below {cnv_min}; detection of exon-level events currently rests on Manta, which needs breakpoints inside mappable sequence.</li>
<li>Mobile element insertions have no dedicated caller. Manta missed real Alu insertions in several carriers; the read-end signature recognises Alu sequence only, was developed on the same two insertion sites it is reported for here, and has not been tested on an independent set or on L1 and SVA elements (Section 3.4).</li>
<li>The ACMG/AMP hints cover a subset of criteria. The PVS1 decision tree, functional (PS3/BS3), segregation (PP1), case-control (PS4) and same-codon (PS1/PM5) evidence are not evaluated.</li>
<li>The pipeline has not been validated for clinical use.</li>
</ul>

<h2>6 Availability and requirements</h2>
<p><strong>Repository:</strong> https://github.com/pjlabsdev/brca12-wgs · <strong>Archived version:</strong> DOI {esc(doi)} · <strong>Operating systems:</strong> Linux, macOS · <strong>Requirements:</strong> Nextflow ≥ 23.10, Java 17, conda/mamba, Docker, Singularity or the listed tools, GRCh38 FASTA with index and dictionary, Ensembl VEP GRCh38 cache, ClinVar VCF · <strong>Licence:</strong> MIT.</p>
<pre>nextflow run main.nf -profile conda \\
  --sample S --reads_r1 S_R1.fastq.gz --reads_r2 S_R2.fastq.gz \\
  --reference GRCh38.fa --vep_cache ~/.vep --clinvar_vcf clinvar.vcf.gz</pre>

<h2>Acknowledgements</h2>
<p>The software and this report were developed with the assistance of an AI coding assistant (Claude, Anthropic). The author reviewed the design and is responsible for the content. This work used data of the Genome in a Bottle Consortium and of the 1000 Genomes Project (New York Genome Center high-coverage sequencing), the exome and amplicon data of the GIAB producers (Oslo University Hospital, University of Nebraska, Thermo Fisher Scientific), the assemblies of the Human Pangenome Reference Consortium and the long reads of the 1KG-ONT cohort. The 1KG-ONT data were generated at the Institute of Molecular Pathology (Vienna, Austria) with funds provided by Boehringer-Ingelheim.</p>
<h2>Conflict of interest</h2>
<p>None declared.</p>

<h2>References</h2>
{{REFS}}
</div>
</body></html>
"""


def missed_text(d):
    """Carriers according to the panel without a read-end signature: say what the reads show."""
    missed = [r for r in d.sv["rows"] if r["type"] == "ALU" and r["panelGenotype"] != "0|0" and not r["readEndSignal"]]
    if not missed:
        return ""
    counts = ", ".join(f"{r['readsAtSite']['clippedAtSite']} of {r['readsAtSite']['reads']} reads in {r['sample']}" for r in missed)
    zero = [r for r in d.sv["rows"] if r["type"] == "ALU" and r["panelGenotype"] == "0|0"]
    top = max(r["readsAtSite"]["clippedAtSite"] for r in zero)
    who = "The one carrier without a signature is" if len(missed) == 1 else f"The {len(missed)} carriers without a signature are"
    return (f"{who} listed as heterozygous in the panel, but only few reads are clipped at "
            f"the site ({counts}), below the threshold of three per side and far from the half expected for a heterozygous insertion; "
            f"samples with panel genotype 0|0 had at most {top}. Whether the panel genotype or the reads are right "
            "cannot be decided here.")


def number_citations(page):
    """Replace [@key, @key] by numbers in order of first appearance and build the reference list."""
    order = []

    def cite(match):
        out = []
        for key in re.findall(r"@([a-z0-9]+)", match.group(0)):
            if key not in REFS:
                raise KeyError(f"unknown reference {key}")
            if key not in order:
                order.append(key)
            n = order.index(key) + 1
            out.append(f'<a class="cite" href="#ref-{n}">{n}</a>')
        return "[" + ", ".join(out) + "]"

    page = re.sub(r"\[@[a-z0-9]+(?:, @[a-z0-9]+)*\]", cite, page)
    unused = [k for k in REFS if k not in order]
    if unused:
        raise ValueError(f"references never cited: {unused}")
    items = "".join(f'<li id="ref-{n}">{REFS[k]}</li>\n' for n, k in enumerate(order, 1))
    return page.replace("{REFS}", f'<ol class="refs">\n{items}</ol>')


def sv_text(d):
    """Sentences of Section 3.4 that depend on the outcome."""
    big, small = d.sv_site("HGSV_216998"), d.sv_site("HGSV_216996")
    alu_c, alu_r = d.sv_site("HGSV_181048"), d.sv_site("HGSV_181052")
    controls_hit = sum(s["controlsReported"] + s["controlsWithReadEndSignal"] for s in d.sv["bySite"].values())
    rows = [r for r in d.sv["rows"] if r["site"] == "HGSV_216998" and r["manta"]]
    shift = sorted({m.get("startDifference") for r in rows for m in r["manta"] if m.get("startDifference") is not None})
    def span(hom):
        v = [x["clippedFraction"] for r in d.sv["rows"] if r["type"] == "ALU" and (r["panelGenotype"] == "1|1") == hom
             and r["panelGenotype"] != "0|0" for x in r["readEndSignal"]]
        return f"{min(v):.2f} to {max(v):.2f}" if v else "n/a"

    het, hom = span(False), span(True)
    parts = [f"The 1,081 bp deletion was reported by Manta in {big['carriersReported']} of {big['carriers']} carriers"
             + (", with the genotype of the panel including the homozygous carrier" if big["carriersReported"] == big["carriers"] else "")
             + (f"; its start lies {abs(shift[0])} bp from the panel coordinate, whose breakpoints are not base-exact" if shift else "") + ".",
             f"The 271 bp deletion was reported in {small['carriersReported']} of {small['carriers']} carriers."
             + (" In the carrier without a call (HG01796) the depth inside the deleted interval equals the depth next to it "
                "(32x and 33x; 13x and 31x in the carrier HG02601), so the reads of this sample do not show the deletion."
                if small["carriersReported"] == 2 and small["carriers"] == 3 else ""),
             f"The Alu insertions were called in {alu_c['carriersReported']} of {alu_c['carriers']} carriers of the common and "
             f"{alu_r['carriersReported']} of {alu_r['carriers']} carriers of the rare insertion, by Manta as insertions (about 320 bp "
             "where the length was resolved) or by the small-variant caller as a partial insertion.",
             f"The read-end signature was present in {alu_c['carriersWithReadEndSignal']} of {alu_c['carriers']} and "
             f"{alu_r['carriersWithReadEndSignal']} of {alu_r['carriers']} carriers, with clipped fractions of {het} in heterozygous "
             f"and {hom} in homozygous carriers according to the panel.",
             missed_text(d),
             ("No sample with panel genotype 0|0 had a call or signature at any of the four sites." if controls_hit == 0 else
              f"{controls_hit} calls or signatures were made in samples with panel genotype 0|0.")]
    return " ".join(x for x in parts if x)


def ont_panel_text(d):
    """Pairs where the long reads and the 1000 Genomes panel disagree."""
    wrong = [r for r in d.ont if (r["longReads"] == "present") != (r["panel"] == "present")]
    if not wrong:
        return "The long reads also agreed with the panel genotypes."
    parts = []
    for r in wrong:
        n = r["supportingReads"].split("/")[1]
        parts.append(f"{r['sample']} at {r['site']} ({'absent' if r['longReads'] == 'absent' else 'present'} in all {n} spanning long reads, "
                     f"{'carrier' if r['panel'] == 'present' else 'non-carrier'} in the panel)")
    return ("The long reads disagreed with the panel only for " + "; ".join(parts) +
            ". This is the carrier of the 271 bp deletion without a call in Section 3.4, whose short-read depth showed no deletion either: "
            "the panel genotype is most likely wrong.")


def ont_all_text(d):
    """Discovery of all events of 50 bp or more in the long reads of the whole cohort."""
    a = d.ont_all
    if not a:
        return ""
    sites = a["sites"]
    by = {s["id"]: s for s in sites}
    def site(chrom, pos, kind):
        return min((s for s in sites if s["chrom"] == chrom and s["type"] == kind), key=lambda s: abs(s["pos"] - pos))
    def share(s):
        n = len(s["carriers"]) + s["absent"]
        return f"{len(s['carriers']):,} of {n:,}"
    alu, tr, dl = site("chr13", 32348296, "INS"), site("chr17", 43128264, "INS"), site("chr17", 43071626, "DEL")
    single = sum(1 for s in sites if len(s["carriers"]) == 1)
    return (f"<p>In addition, all insertions and deletions of 50 bp or more inside both gene loci were collected from the long reads of all "
            f"{a['samples']:,} samples of the cohort and every sample was genotyped at every site ({len(sites)} sites, none in an exon ±20 bp; "
            f"{single} seen in a single sample). The intronic <i>BRCA2</i> Alu insertion was present in {share(alu)} evaluable samples, "
            f"the 1,081 bp <i>BRCA1</i> deletion in {share(dl)} and the 126 bp tandem-repeat insertion of Section 3.7 in {share(tr)}. "
            "The last one is a common allele, which short-read callers miss and the tandem-copy signature finds only in part (Section 3.7).</p>")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fonts")
    p.add_argument("--version", default="1.1.0 (in preparation)")
    p.add_argument("--date", default="13 November 2026")
    p.add_argument("--doi", default="10.5281/zenodo.24173638 (version 1.0.1)")
    args = p.parse_args(argv)
    d = Data()
    page = build(d, args.fonts, args.version, args.date, args.doi)
    nomates = load(os.path.join(RES, "sv_carriers_summary.without_mates.json"))["bySite"]
    alu_nomates = sum(nomates[s]["carriersReported"] for s in ("HGSV_181048", "HGSV_181052"))
    alu_total = sum(nomates[s]["carriers"] for s in ("HGSV_181048", "HGSV_181052"))
    sig_nomates = sum(nomates[s]["carriersWithReadEndSignal"] for s in ("HGSV_181048", "HGSV_181052"))
    page = page.replace("{SV_NOMATES}", f"{alu_nomates} of {alu_total} Alu insertions were called and {sig_nomates} of {alu_total} "
                                        "showed the read-end signature")
    page = page.replace("{SV_TEXT}", sv_text(d))
    page = page.replace("{ONT_PANEL}", ont_panel_text(d))
    page = page.replace("{ONT_ALL}", ont_all_text(d))
    page = number_citations(page)
    sys.stdout.write(page)
    return 0


if __name__ == "__main__":
    sys.exit(main())

# Validation against assembly-based truth (HPRC)

Research use only. Not a clinical validation.

## What is tested

39 individuals of the Human Pangenome Reference Consortium, year 1, are in the 1000 Genomes NYGC
30x short-read data. For them the truth comes from haplotype-resolved long-read assemblies (HPRC
v1.1 minigraph-cactus pangenome, GRCh38 VCF decomposed with vcfwave). The truth is independent of
short-read callers and covers the whole of BRCA1 and BRCA2, including introns, repeats and
variants of 50 bp or more.

- `hprc_samples_1000g.tsv`: the 39 samples and their 1000 Genomes CRAM files.
- `run_hprc.sh WORKDIR`: cuts the gene regions out of the pangenome VCF, runs the pipeline from
  read pairs (as in `pathogenic/run_sv_carriers.sh`) and evaluates.
- `hprc_eval.py`: small variants (below 50 bp) per sample with RTG vcfeval; positions with a
  missing allele in the truth and the surroundings of larger variants are excluded. Variants of
  50 bp or more are grouped into sites and looked up in the Manta calls, the small-variant calls
  and the read-end insertion signatures.
- `results/hprc_summary.{tsv,json}`: totals, per-sample counts and the large-variant sites.

Stand-ins: VEP passthrough (no cache) and neutral CNVkit segments; neither changes the calls.

## Results

Small variants, 39 samples together:

| Region | Type | Truth | TP | FP | FN |
|---|---|---:|---:|---:|---:|
| Exons ± 20 bp | SNV | 590 | 590 | 0 | 0 |
| Exons ± 20 bp | Indel | 86 | 68 | 10 | 18 |
| Gene loci | SNV | 8,095 | 8,026 | 100 | 69 |
| Gene loci | Indel | 3,312 | 2,494 | 392 | 818 |

16 of the 18 missed exon indels and 9 of the 10 false exon indels lie in one poly-T stretch in the
BRCA1 3'UTR (chr17:43044804-43044817). There the assemblies and the short reads disagree on the
length of the homopolymer; the coding sequence is not involved. Outside this stretch the exon
indels are 41/43 with 1 false positive (the poly-T stretch itself: 27/43 found, 9 false). The gene-locus figures
include all introns and repeats and have not been curated site by site.

Variants of 50 bp or more:

| Site | Event | Carriers found | Calls in non-carriers |
|---|---|---|---|
| chr13:32348295 (BRCA2 intron) | Alu insertion, 318 bp | 6/6 with read-end signal, 0/6 Manta call | 0/33 |
| chr17:43100595 (BRCA1 intron) | VNTR, -80 to +113 bp | 8/14 called | 5/25 |
| chr17:43128264 (BRCA1 intron) | 126 bp tandem-repeat insertion | 0/19 | 0/20 |

The 126 bp insertion is a further copy of an adjacent tandem repeat. Its reads align to the
existing copy with soft-clipped ends only, so neither Manta nor the read-end signature (which
needs clipped sequence with Alu k-mers) reports it. This is a known blind spot of short reads
and is listed in the limitations of the report.

## Sources

Liao WW et al. A draft human pangenome reference. Nature 2023;617:312-324.
Byrska-Bishop M et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes
Project cohort including 602 trios. Cell 2022;185:3426-3440.

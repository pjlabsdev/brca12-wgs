# GIAB benchmark (HG001, HG002)

Small-variant calls of brca12-wgs compared with the Genome in a Bottle (GIAB) NIST v4.2.1 truth sets
for HG001 (NA12878) and HG002 (NA24385) on GRCh38, inside BRCA1 and BRCA2.

Run on 2026-11-13 with brca12-wgs at the state of this branch (including the review fixes of 1.0.2),
GATK HaplotypeCaller 4.6.1.0, BWA 0.7.18, samtools/bcftools 1.21, mosdepth 0.3.10, Manta 1.6.0,
Nextflow 24.10.4 and RTG Tools 3.13 (`vcfeval`).

## Results

Matching with RTG `vcfeval` (haplotype-aware, the genotype must agree). Evaluated only inside the
GIAB high-confidence regions: 91.9 % (HG001) and 95.1 % (HG002) of the 211,225 bp of the two gene
loci. 95 % confidence intervals (Wilson) in brackets.

| Sample | Region | Type | Truth | TP | FP | FN | Recall | Precision |
|---|---|---|---:|---:|---:|---:|---|---|
| HG001 | gene loci | SNV | 269 | 269 | 1 | 0 | 1.000 [0.986, 1.000] | 0.996 [0.979, 0.999] |
| HG001 | gene loci | indel | 49 | 48 | 1 | 1 | 0.980 [0.893, 0.996] | 0.980 [0.893, 0.996] |
| HG001 | exons ±20 bp | SNV | 20 | 20 | 0 | 0 | 1.000 [0.839, 1.000] | 1.000 [0.839, 1.000] |
| HG001 | exons ±20 bp | indel | 2 | 2 | 0 | 0 | 1.000 | 1.000 |
| HG002 | gene loci | SNV | 95 | 95 | 0 | 0 | 1.000 [0.961, 1.000] | 1.000 [0.961, 1.000] |
| HG002 | gene loci | indel | 20 | 19 | 1 | 1 | 0.950 [0.764, 0.991] | 0.950 [0.764, 0.991] |
| HG002 | exons ±20 bp | SNV | 11 | 11 | 0 | 0 | 1.000 [0.741, 1.000] | 1.000 [0.741, 1.000] |
| HG002 | exons ±20 bp | indel | 1 | 1 | 0 | 0 | 1.000 | 1.000 |

Across both samples: 364 of 364 truth SNVs and 67 of 69 truth indels found; in the coding exons
±20 bp all 34 truth variants were found without a false positive.

**The three discordant calls**

- HG001 chr13:32344166 (BRCA2 intron, poly-A): truth is a multiallelic `GAA>GAAA,G` (1/2), the
  pipeline called only the deletion (0/1). It counts as one FP and one FN. Depth 12x, GQ 1.
- HG002 chr13:32349814 (BRCA2 intron, homopolymer): truth homozygous 1-bp deletion, called
  heterozygous. One FP and one FN (wrong genotype). Depth 14x, GQ 7.
- HG001 chr17:43164525 G>T (BRCA1 upstream intron, not coding): heterozygous call (34x, GQ 99) that
  is absent from the truth set.

Both indel discordances carry the pipeline's own quality flags in the report (DP < 15, GQ < 20,
unusual allele balance), so a reader is warned. The SNV false positive carries no flag.

**Coverage per exon:** HG001 (≈36x) had no uncertain exon (exon mean depth 23 to 53x). HG002
(≈30x) had three exons flagged as uncertain: BRCA1 exon 1, BRCA1 exon 21 (legacy 22, mean 17.5x)
and BRCA2 exon 1. No exon had more than 10 % reads with MAPQ < 20.

**Structural variants and CNVs:** Manta reported no SV in either sample, as expected for these
reference samples. CNVs were not benchmarked (see below).

Files: `results/HG001.summary.tsv`, `results/HG002.summary.tsv` (+ `.json` with confidence
intervals), `results/exon_coverage.json`.

## How it was run

`validation/run_giab.sh WORKDIR` does everything:

1. **Reference:** chr13 and chr17 of the GIAB GRCh38 no-alt analysis set, with a BWA index,
   `.dict` and RTG SDF.
2. **Reads:** GIAB/NHGRI Illumina HiSeq 2×148 bp data from the public GIAB S3 mirror.
   - Read pairs mapped to BRCA1/BRCA2 ±10 kb are extracted and downsampled: HG001 300x to ≈36x,
     HG002 60x to ≈30x.
   - They are written back to FASTQ, so the pipeline runs from FASTQ, BWA-MEM included.
3. **Pipeline:** the full Nextflow pipeline with real BWA, samtools, GATK, bcftools, mosdepth and
   Manta.
4. **Evaluation:** `validation/giab_benchmark.py` runs `rtg vcfeval` on
   `calls/<sample>.small.vcf.gz`. Evaluation regions are GIAB high-confidence ∩ gene loci, and
   GIAB high-confidence ∩ MANE exons ±20 bp.

## Limitations of this benchmark

- **Targeted read set.** Only reads that GIAB's novoalign mapping placed near BRCA1/BRCA2 are used,
  and they are realigned to chr13 + chr17 only.
  - Reads from paralogous sequence elsewhere in the genome cannot compete for alignment, which can
    flatter precision.
  - Reads that the original mapping placed elsewhere are missing, which can lower recall.
  - A whole-genome run against the full GRCh38 reference is the next step.
- **Annotation and CNVs not benchmarked.** No Ensembl VEP cache was reachable from the benchmark
  machine.
  - VEP was replaced by `vep_passthrough.sh`, which leaves variants unannotated.
  - CNVkit was replaced by the test stand-in, which reports neutral segments, because a read set
    restricted to the genes has no meaningful genome-wide depth baseline.
  - The `cnvkit` and `ensembl-vep` versions in these reports are therefore not real.
- **Few variants.** Two samples and a small region give few events, especially exonic indels
  (3 in total). The confidence intervals are wide.
- **Non-confident regions not evaluated.** About 5 to 8 % of the gene loci lie outside the GIAB
  high-confidence regions.
- **No pathogenic variants.** Neither sample carries a known pathogenic BRCA1/BRCA2 variant. Real
  pathogenic SNVs/indels, large deletions/duplications (MLPA-confirmed) and Alu insertions need
  characterised samples, for example from Coriell or the GeT-RM BRCA reference materials.

# Validation on exome and amplicon data

Research use only. `run_targeted.sh` re-aligns three public GIAB data sets from FASTQ to GRCh38 and
runs brca12-wgs with `--assay exome` or `--assay amplicon` and `--mosaic`; `validation/giab_benchmark.py`
compares the small-variant calls with GIAB v4.2.1 (vcfeval) in three region sets: the gene loci, the
MANE exons ±20 bp and the targets, each intersected with the GIAB confident regions.

| Data set | Producer, kit | Reads | Targets |
|---|---|---|---|
| HG002_exome_agilent | Oslo University Hospital, Agilent SureSelect v5 | HiSeq 2x125, paired | MANE exons ±50 bp (kit BED not public) |
| HG001_exome_truseq | University of Nebraska, Illumina TruSeq exome | paired | vendor BED, lifted to GRCh38 |
| HG002_ampliseq_ion | Thermo Fisher, Ion AmpliSeq Exome, one of four runs | Ion Proton, single-end | AmpliSeq effective regions, lifted |

hg19 intervals are moved to GRCh38 with one offset per locus (BRCA1 +1,847,983, BRCA2 −574,137;
Ensembl assembly mapping, one block each). Annotation is a passthrough (no VEP cache) and not part of
the benchmark. Results: `results/*.benchmark.tsv` (vcfeval counts) and `results/*.qc.json` (read and
target QC, uncertain exons, low-coverage intervals, number of mosaic candidates).

## Results inside the targets

| Data set | SNVs | SNV FP | Indels | Indel FP |
|---|---|---|---|---|
| HG002 Agilent exome | 10/12 | 0 | 0/1 | 0 |
| HG001 TruSeq exome | 15/17 | 0 | 1/2 | 1 |
| HG002 Ion AmpliSeq | 8/8 | 3 | 0/0 | 98 |

- All four missed SNVs and both missed indels lie in the 3' UTR of BRCA2 exon 27
  (chr13:32399139-32399885), where the original exome BAMs have 0 to 2 reads although the TruSeq BED
  lists the region. The report lists this region as an uncertain interval (depth below 10x).
- The TruSeq indel FP is chr13:32399599 AT>A called homozygous from 2 reads (truth heterozygous).
- Ion Torrent: the 3 SNV FPs have alternative-allele fractions of 10 to 20 % and carry the
  allele-balance flag; the 98 indel FPs are homopolymer errors of the platform. With
  `--platform iontorrent` every indel carries a warning and mosaic candidates are SNVs only.
- No mosaic candidate in any of the three data sets with the default filters (≥ 3 supporting reads,
  ≥ 20x). These filters were chosen after a first run on the same data (4 candidates for TruSeq at
  8x with 2 supporting reads each, 23 Ion indels), so this is not an independent test of them.
- Target QC: Agilent mean target depth 216x, uniformity 83.9 %, fold-80 3.43; TruSeq 91x, 88.7 %,
  2.22; AmpliSeq 200x, 91.6 %, 2.53, error rate 1.27 % (flagged).

Outside the targets (introns) exome data have few reads, so the gene-locus numbers in the TSV files
are low by design and not a measure of the caller.

## Reproduce

```bash
validation/run_giab.sh WORK            # reference chr13+chr17, SDF and GIAB truth (once)
validation/targeted/run_targeted.sh WORK/targeted WORK/data
```

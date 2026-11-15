# Preregistration: GIAB HG003 to HG007 before release 1.1.0

Frozen on 2026-11-15, before any read of these five samples was run through the pipeline.
Code under test: `main` at commit ce825de (Delly 1.3.3 as default SV caller).

## Why

The published GIAB benchmark covers HG001 and HG002 only. Before release 1.1.0 the pipeline is
tested on five more GIAB reference genomes, from two further ancestries (Ashkenazi trio parents,
Han Chinese trio), that were never used for tuning.

## Samples

| Sample | Reads (GIAB S3) | Downsampling | Truth |
|---|---|---|---|
| HG003 | NHGRI novoalign `HG003.GRCh38.60x.1.bam` | `-s 42.5` (~30x) | v4.2.1, `_noinconsistent.bed` |
| HG004 | NHGRI novoalign `HG004.GRCh38.60x.1.bam` | `-s 42.5` (~30x) | v4.2.1, `_noinconsistent.bed` |
| HG005 | NHGRI novoalign `...minus_alts.300x.bam` | `-s 42.1` (~30x) | v4.2.1, `benchmark.bed` |
| HG006 | NHGRI novoalign `...minus_alts.100x.bam` | `-s 42.3` (~30x) | v3.3.2 `_noinconsistent` (v4.2.1 is not on the GIAB S3 mirror) |
| HG007 | NHGRI novoalign `...minus_alts.100x.bam` | `-s 42.3` (~30x) | v3.3.2 `_noinconsistent` (as HG006) |

HG008 has no germline benchmark and HG009/HG010 are not public; they are not part of this test.

## Arm A: targeted benchmark (cloud), identical to the published HG001/HG002 method

`validation/giab_more/run.sh`, the same steps as `validation/run_giab.sh`: read pairs over
BRCA1/BRCA2 ±10 kb from the BAM, downsampled, back to FASTQ, pipeline from FASTQ (BWA-MEM,
GATK HaplotypeCaller, Delly), chr13+chr17 reference; VEP and CNVkit are stand-ins.
Comparison with `validation/giab_benchmark.py` (hap.py-compatible RTG vcfeval, GIAB confident
regions inside BRCA1/BRCA2).

## Arm B: whole genome on the Mac

GIAB/Google NovaSeq PCR-free 30x FASTQ of HG003 to HG007, whole-genome GRCh38 alignment, full
pipeline with real VEP 116.2 and CNVkit. Same comparison script and the same thresholds. Truth
v4.2.1 for all five where available from NCBI. Arm B runs when the downloads (~266 GB) are done.

## Thresholds (per sample, inside GIAB confident regions of BRCA1/BRCA2), fixed now

- SNVs: recall ≥ 0.99 and precision ≥ 0.99.
- Indels: recall ≥ 0.90 and precision ≥ 0.90.
- Exons ±20 bp: every truth variant found, no false positive.
- Structural variants: every Delly or exon copy-number call inside the genes is reported. There is
  no SV truth set here, so a call is not a failure by itself; each call is checked in the GIAB
  long reads (PacBio HiFi) where available.

The test passes if every sample meets the small-variant thresholds in arm A. Every miss and false
positive is listed with its cause. No parameter is changed after the first result; if the test
fails, the cause is fixed in a separate change and the five samples are reported both before and
after, clearly labelled.

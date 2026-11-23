# Arm B in the cloud: HG003 to HG007 with new NovaSeq reads (BRCA1/BRCA2 only)

Arm B of the GIAB benchmark (`validation/giab_more/PREREGISTRATION.md`,
`validation/giab_more/ADDENDUM_arm_b.md`) planned whole-genome runs of HG003 to HG007 on pj's Mac from
the Google NovaSeq 6000 PCR-free 30x FASTQ. The Mac run stopped repeatedly (external disk dropping
out, closed windows; HG003 reached 15 of 22 alignment chunks). On 2026-11-23 pj chose to run arm B
in the cloud with only the reads over BRCA1/BRCA2, because whole genomes do not fit the cloud
environment (4 cores, about 4.5 GB free disk, background jobs of at most 2 hours).

## Method (`run.sh`)

- Reads: the same NovaSeq 30x reads, taken from Google's public GRCh38 BAMs of those reads
  (`brain-genomics-public/research/sequencing/grch38/bam/novaseq/wgs_pcr_free/30x`) through the
  index: read pairs over BRCA1 and BRCA2 ±10 kb, turned back into FASTQ, so the pipeline runs from
  FASTQ (BWA-MEM alignment included), as in arm A. No downsampling (the data are already 30x).
- Reference: chr13 + chr17 of the GIAB GRCh38 no-alt analysis set (as arm A).
- Truth: GIAB v4.2.1 for all five samples (as in the addendum; arm A used v3.3.2 for HG006 and HG007).
- Comparison: `validation/giab_benchmark.py` (rtg vcfeval) in the GIAB confident regions of the gene
  loci and of the exons ±20 bp; thresholds of the preregistration (SNV recall and precision ≥ 0.99,
  indel ≥ 0.90, every exonic variant found and no false positive in the exons).
- Code: the frozen arm B version `main` a705662, and as a second run the current `main` 4122e71
  (release 1.1.1 changes plus the result files in the form of hla-wgs).

## Deviations from the plan

1. Only the reads over BRCA1/BRCA2 ±10 kb, aligned to chr13 + chr17, instead of whole-genome
   alignment of all reads to the full GRCh38 with decoys. Reads from elsewhere in the genome that
   would map here (or reads from here that would map elsewhere) are not represented.
2. VEP and CNVkit are stand-ins (passthrough, neutral segments), as in arm A: annotation and
   CNVkit copy number are not tested. Delly runs for real.
3. Truth for HG006/HG007 is v4.2.1 (as the addendum says), not v3.3.2 as in arm A.

So this arm tests the small-variant calling in BRCA1/BRCA2 on new reads of the same five people
("same people, new reads"), not the whole-genome path.

## Result

Both code versions give exactly the same benchmark result (`results/summary_a705662.tsv`,
`results/summary_main.tsv`).

| Sample | SNVs TP/truth, FP | Indels TP/truth, FP | Exons ±20 bp |
|---|---|---|---|
| HG003 | 263/263, 1 FP | 47/47, 0 FP | 22/22, 0 FP |
| HG004 | 86/86, 0 FP | 21/21, 0 FP | 12/12, 0 FP |
| HG005 | 255/255, 0 FP | 32/32, 0 FP | 23/23, 0 FP |
| HG006 | 253/253, 0 FP | 46/46, 0 FP | 24/24, 0 FP |
| HG007 | 117/117, 0 FP | 21/21, 0 FP | 13/13, 0 FP |
| All | 974/974, 1 FP | 167/167, 0 FP | 94/94, 0 FP |

**The preregistered test is passed by all five samples.** The single false positive (HG003,
chr17:43084901 C>G, intron, 9 of 49 reads) lies at the same position as HG003's false positive in
arm A, which came from different reads (HiSeq); it may be a variant missing from the truth set, but
this is not shown here. No structural variant or copy number change was called. In arm A, HG004
(precision 0.977) and HG005 (recall 0.988) had missed the SNV threshold; with these reads and the
poly-G fix neither happens.

With the current `main` every sample gets `result` = `none` in both genes (`results/HG00?.tsv`); the
known common alleles (the 126 bp copy upstream of BRCA1 exon 1 in four samples, the BRCA2 intronic
Alu insertion in HG007) are listed but no longer raise attention, while a705662 still counted them.
Coverage is `limited` in at least one gene of every sample: these NovaSeq reads give
24 to 33x mean exon depth, and exons below 20x or with less than 95 % of bases at 20x are marked
uncertain (HG006 BRCA2: 23.6x, 18 of 27 exons). A "no finding" there does not exclude a variant in
those exons, and the page says so.

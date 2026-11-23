# Validation of brca12-wgs 2.0 (in development)

Small variants from the program's own caller (no GATK), compared with RTG vcfeval in the gene loci
and in the coding exons +/- 20 bp, against the same data and truth sets as the 1.1 pipeline (GATK
HaplotypeCaller). Research use only.

## GIAB HG001 (1000 Genomes NA12878, 30x)

| region | type | truth | 2.0 TP | 2.0 FP | 2.0 FN | 1.1 TP | 1.1 FP | 1.1 FN |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| gene loci | SNV | 269 | 269 | 1 | 0 | 269 | 1 | 0 |
| gene loci | indel | 49 | 49 | 0 | 0 | 48 | 1 | 1 |
| exons +/- 20 bp | SNV | 20 | 20 | 0 | 0 | 20 | 0 | 0 |
| exons +/- 20 bp | indel | 2 | 2 | 0 | 0 | 2 | 0 | 0 |

The one SNV false positive (chr17:43164525 G>T, intronic, 34x, both strands) is the same in 1.1: a
well-supported call that the GIAB truth does not contain.

## HPRC, 39 samples (assembly-based truth, 1000 Genomes 30x reads)

| region | type | truth | 2.0 TP | 2.0 FP | 2.0 FN | 1.1 TP | 1.1 FP | 1.1 FN |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| gene loci | SNV | 8095 | 8048 | 110 | 47 | 8026 | 100 | 69 |
| gene loci | indel | 3312 | 2520 | 388 | 792 | 2494 | 392 | 818 |
| exons +/- 20 bp | SNV | 590 | 590 | 0 | 0 | 590 | 0 | 0 |
| exons +/- 20 bp | indel | 86 | 68 | 11 | 18 | 68 | 10 | 18 |

Most indel errors of both versions are homopolymer lengths, where short reads and the assembly-based
truth often disagree.

## Known pathogenic carriers (1000 Genomes)

The whole analysis (`brca12-wgs INPUT OUTPREFIX --reference`) on six carriers: each variant was
found, named correctly and gave `result = pathogenic` for its gene and `none` for the other gene.

| sample | gene | variant | 2.0 statement |
|---|---|---|---|
| HG02164 | BRCA1 | c.34C>T | pathogenic: NM_007294.4:c.34C>T p.(Gln12Ter), heterozygous |
| HG00365 | BRCA1 | c.4357+1G>A | pathogenic: NM_007294.4:c.4357+1G>A, heterozygous |
| HG02274 | BRCA2 | c.5146_5149del | pathogenic: NM_000059.4:c.5146_5149del p.(Tyr1716LysfsTer8), heterozygous |
| NA19130 | BRCA2 | c.5241_5242insTA | pathogenic: NM_000059.4:c.5241_5242insTA p.(Ser1748Ter), heterozygous |
| HG01377 | BRCA2 | c.7063G>T | pathogenic: NM_000059.4:c.7063G>T p.(Glu2355Ter), heterozygous |
| HG00611 | BRCA2 | c.7697dup | pathogenic: NM_000059.4:c.7697dup p.(Asp2566GlufsTer5), heterozygous |

## HGVS and consequence

Against Ensembl VEP 116 (`--hgvs --shift_hgvs 1`) on all 37,125 ClinVar alleles in BRCA1/BRCA2:
c. equal 99.97 %, p. equal 99.96 %, most severe consequence equal 99.70 %
(details in docs/output-2.md).

## Exon copy number

Same results as 1.1 on the same input: on a 30x whole genome all 58 interval depths equal
`samtools depth -a -Q 20` and all 50 exon rows (log2 deviation, copy number, posterior) equal the
Python implementation, with and without allele balance; with half of the read pairs removed over
BRCA1 exons 11-13 and BRCA2 exons 14-16, both implementations call the same segments with the same
quality, and none in the unmodified genome. The detection rates of the method are those of the 1.1
validation (validation/exon_cnv). Depth for the panel takes under 1 s instead of about 4 min.

## Structural variants and insertion signatures

The Rust module gives the same JSON as 1.1's `brca12_insertion_signals.py --gene-fasta
--read-pair-svs` on all 53 samples tested (39 HPRC, NA12878, six 1000 Genomes carriers, seven
spike-ins): 46 + 7 files identical, among them 7 known Alu insertions and 23 known 126 bp copies.

Spike-ins as in validation/pathogenic (HG002 chr13/chr17 background, simulated heterozygous events
in the gene windows, BWA-MEM), events found by the read-pair/split-read module:

| sample | event | found |
|---|---|---|
| L1 | BRCA1 deletion 510 bp, exon 21 | yes, exact breakpoints |
| L1 | BRCA2 c.156_157insAlu | yes, as the known founder allele (gene result pathogenic) |
| L2 | BRCA1 deletion 3,835 bp, exon 12 | yes |
| L2 | BRCA2 deletion 8.5 kb, exons 14-16 | yes, exact |
| L3 | BRCA1 tandem duplication 6 kb, exon 12 | yes |
| L3 | BRCA2 deletion 1.5 kb, exon 21 | yes, exact |
| L4 | BRCA1 deletion 11 kb, exons 1-2 | yes |
| L4 | BRCA2 deletion 120 bp in exon 11 | yes |
| L5 | BRCA1 whole-gene deletion (135 kb) | no breakpoint in the gene region; the exon copy number shows all exons at one copy |
| L5 | BRCA2 tandem duplication 5 kb, exons 19-20 | yes |
| L6 | BRCA1 deletion 300 bp in exon 10 | yes |
| L6 | BRCA2 whole-gene deletion (95 kb) | no breakpoint in the gene region; exon copy number one copy for exons 1-7 and 13-18 |
| negative control | none | no signature |

1.1 for comparison: Manta 12/12, CNVkit 3/12. The exon copy number cannot be judged on these
spike-ins: the gene windows are simulated reads in HiSeq/novoalign data, which the panel of 1000
Genomes NovaSeq genomes does not fit (false duplication calls also in the negative control).

## Not covered yet

Gene-level CNV segmentation as CNVkit did (CNVkit found 3 of the 12 spike-ins in 1.1; whole-gene
events are seen by the exon copy number); an assay check for the exon copy-number panel; large insertions that the read
aligner soft-clips (for example a 32 bp intronic insertion that GATK assembles) are not found yet.

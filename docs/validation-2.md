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

## Not covered yet

Exon copy number, CNV and structural variants (as in 1.1, phase 3); large insertions that the read
aligner soft-clips (for example a 32 bp intronic insertion that GATK assembles) are not found yet.

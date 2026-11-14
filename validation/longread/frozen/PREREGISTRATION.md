# Frozen test of known-allele genotyping (written before the data were fetched)

Research use only. Written on 2026-11-14 before any short-read data of the samples below were
downloaded or looked at. Nothing in the method may change after this point; any later change is
reported as a post-hoc analysis next to the frozen result.

## Method under test

`bin/brca12_insertion_signals.py --gene-fasta` at commit 6a03091e8a90239979e99e953570af52ce20e6c6
(branch `known-allele-realignment`), unchanged defaults:

- Allele-aware realignment: reads starting 400 bp before to 300 bp after the site (unmapped mates
  placed there included, at least 50 bp) are aligned to the reference haplotype and to the allele
  haplotype (700 bp flanks, HPRC allele sequences). A read supports a haplotype when it fits it with
  at least 2 fewer edits; reads with more than 8 % edits to both are ignored. Call at 2 or more
  allele-supporting reads.
- 126 bp tandem copy (chr17:43128264): reported on realignment alone, or on 2 of: clipped reads at
  both unit borders, marker base chr17:43128365 C, linked haplotype (3 of 5 tag SNPs).
- BRCA2 intronic Alu insertion (chr13:32348295): reported on realignment; zygosity homozygous when
  reference-supporting reads are at most 10 % of the informative reads.

## Samples

200 samples of the 1KG-ONT cohort (Schloissnig et al. 2025) with 1000 Genomes 30x short reads that
were not part of any earlier analysis of this repository (the 72 samples of `validation/hprc` and
`validation/pathogenic` and the 120 samples of `results/tandem_copy_independent.tsv`), drawn with
`random.seed(20261114); random.sample(pool, 200)` from the 834 eligible samples (sorted list of the
1,019 ONT CRAM index names that are in the two 1000 Genomes sequence indices). List:
`samples.tsv`, SHA-256 d1ad88e0c0b583e33640cdfa30c54d2ad72a2964fa2228228f6c2498d37d1854.

Short reads: original 1000 Genomes alignments (`samtools view` of chr17:43126000-43150000 and
chr13:32346000-32351000).

## Truth

Long-read genotypes of `ont_discover.py` (`results/all_disc.json`, sites ONT_BRCA1_INS_43128267 and
ONT_BRCA2_INS_32348296): carrier when listed as carrier, otherwise non-carrier. The discovery run
does not list which samples were not evaluable at a site (17 to 26 of 1,019 per site), so a few
non-carriers may be unevaluable samples; this can only lower the measured specificity or
sensitivity, not raise it. Zygosity: homozygous when at least 90 % of the spanning long reads carry
the insertion.

## Endpoints

Primary, per allele: sensitivity and specificity of the reported known allele with exact 95 %
confidence intervals (Clopper-Pearson). Secondary: the same for the read-end signature alone (the
previous method) on the same samples, and zygosity agreement with the long reads. Evaluation script:
`../frozen_eval.py` (in this commit). Samples whose short reads cannot be fetched are listed, not
replaced.

# Long-read check of structural variants (1KG-ONT)

Research use only. Not a clinical validation.

The structural-variant sites of `pathogenic/sv_sites_1000g.tsv` come from short-read callsets of
the 1000 Genomes panel. Here they are checked in Oxford Nanopore long reads of the same
individuals (1KG-ONT Vienna cohort, Schloissnig et al. 2025, GRCh38 CRAMs read remotely, only the
two gene regions).

- `ont_genotype.py`: a read counts when it spans the site with a margin on both sides; it supports
  the variant when its alignment carries an insertion (Alu) or deletion of at least half the
  expected length near the site. Present when at least 2 spanning reads and at least 10 % of
  them support it.
- `ont_discover.py`: finds every insertion and deletion of 50 bp or more in BRCA1/BRCA2 in all
  cohort samples, clusters them into sites and genotypes every sample at every site.
- `compare_shortread.py`: joins the long-read genotypes with the panel genotypes and the
  pipeline results of `pathogenic/run_sv_carriers.sh`.

## Results for the carrier samples (`results/`)

11 of the carrier and control samples of `pathogenic/sv_carriers_1000g.tsv` have long reads
(`results/ours_ont.txt`).

- The pipeline (Manta call or read-end insertion signal) agrees with the long reads in 26 of 26
  sample-site pairs (`results/ours_vs_shortread.tsv`).
- Long reads confirm the intronic BRCA2 Alu insertion in HG00117, HG00637 and HG01377 and the
  1,081 bp BRCA1 deletion (980 bp in the reads) in HG00637 and HG01990.
- The 271 bp BRCA1 deletion that the panel lists for HG01796 is absent in all 26 spanning long
  reads; the pipeline did not report it either. The panel genotype is most likely wrong.

## All 1,019 samples (`results/all_disc.*`)

`ont_discover.py` collected every insertion and deletion of 50 bp or more inside both gene loci from
the long reads of all 1,019 cohort samples (no sample failed) and genotyped every sample at every
site: 18 sites, none in an exon ±20 bp, 9 seen in a single sample.

- BRCA2 intronic Alu insertion (chr13:32348296, 316 bp): 136 of 1,002 evaluable samples.
- BRCA1 1,081 bp intronic deletion (980 bp in the reads): 13 of 993.
- BRCA1 VNTR at chr17:43100560-43100705: several alleles from -77 to +568 bp.
- 126 bp tandem-repeat insertion upstream of BRCA1 exon 1 (chr17:43128267): 552 of 1,004. The
  short-read callers miss it; the tandem-copy signature finds it in 5 of 6 long-read
  carriers among the 11 tested samples and in none of the 5 non-carriers
  (`results/tandem_copy.tsv`, `results/tandem_copy_truth.tsv`).

## Independent test of the tandem-copy signature (`results/tandem_copy_independent*.tsv`)

The rules of the signature were set on the 39 HPRC and 11 ONT samples above. As an independent
test, 120 further 1KG-ONT samples with 1000 Genomes 30x data were taken (none used before): 67
carry the 126 bp insertion in the long reads (`ont_genotype.py`), 53 do not. The short reads are the
original 1000 Genomes alignments of chr17:43126000-43150000; the evaluation is
`../hprc/tandem_copy_eval.py --sam`.

| Signal | Carriers (67) | Non-carriers (53) |
|---|---:|---:|
| Clipped reads at both unit borders | 60 | 1 |
| Marker base chr17:43128365 C | 65 | 2 |
| Linked haplotype (at least 3 of 5 tag SNPs) | 67 | 1 |
| Reported: at least 2 of 3 | 67 | 1 |

The one non-carrier with a report (HG02687) shows all three signals (4+4 clipped reads, 4 of 43
reads with C, 5 of 5 tag SNPs); it has only 8 spanning long reads, 1 of them with the insertion, so
it may well be a carrier missed by the long-read genotype. The haplotype count gives the zygosity:
all 18 carriers with two haplotype copies show the insertion in every spanning long read, and 48
of the 49 with one copy show it in 22 to 83 % of them.

Sites outside the gene loci (reads reach beyond them) are dropped, and a shorter insertion cluster
with the same carriers next to a longer one (an insertion aligned in two pieces) is merged into it.

## Preregistered test of known-allele genotyping (`frozen/`, `results/frozen_*`)

Known alleles from the assemblies (the 126 bp copy and the intronic BRCA2 Alu insertion) are also
genotyped by allele-aware realignment (`bin/brca12_insertion_signals.py`): each read near the site
is aligned to the reference haplotype and to the haplotype with the allele, and two reads that fit
the allele with at least 2 fewer edits make a call. Because all rules so far were developed on the
samples above, the final method was frozen and the test written down and pushed before any data were
fetched (`frozen/PREREGISTRATION.md`, commit 9acce4d): 200 further 1KG-ONT samples drawn with a fixed
seed from the 834 not used before, original 1000 Genomes alignments, long-read genotypes of
`results/all_disc.*` as truth. Evaluation: `frozen_eval.py`, per sample `results/frozen_test.tsv`.

| Allele | Method | Carriers found | Non-carriers with a call |
|---|---|---:|---:|
| 126 bp tandem copy | frozen method | 116/116 (95 % CI 96.9-100 %) | 4/84 |
| 126 bp tandem copy | read-end signature alone | 97/116 | 2/84 |
| BRCA2 Alu insertion | frozen method | 25/25 (95 % CI 86.3-100 %) | 0/175 |
| BRCA2 Alu insertion | read-end signature alone | 19/25 | 0/175 |

Genotype agreement with the long reads (homozygous when at least 90 % of spanning long reads carry
the insertion): Alu 25/25, tandem copy 113/116. Post hoc, the 4 non-carriers with a call were
genotyped again from the long reads (`results/frozen_audit_ont.tsv`): they have only 3 to 6
spanning long reads, and 3 of them have at least one long read with the insertion. They may be gaps
in the truth, but the frozen result counts them as false calls.

By superpopulation (1000 Genomes labels), carriers found for the tandem copy: AFR 27/27, AMR
20/20, EAS 21/21, EUR 17/17, SAS 31/31 (non-carrier calls 0/29, 2/12, 1/17, 1/16, 0/10); for the Alu
insertion every carrier in every group (2 to 8 per group) and no call in any non-carrier.

## Sources

Schloissnig S et al. Structural variation in 1,019 diverse humans based on long-read sequencing.
Nature 2025. Data generated at the Institute of Molecular Pathology (Vienna, Austria) with funds
provided by Boehringer-Ingelheim.

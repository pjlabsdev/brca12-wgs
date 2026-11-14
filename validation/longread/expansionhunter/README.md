# ExpansionHunter comparison for the 126 bp tandem copy (BRCA1)

Research use only. Not a clinical validation.

## What is compared

The 126 bp insertion upstream of BRCA1 exon 1 adds a unit with the sequence of the second
reference copy (GRCh38 chr17:43128247-43128372, with C instead of T at chr17:43128365). As an
established repeat genotyper, ExpansionHunter (Dolzhenko et al. 2019, see `references.md`) was run
on the same short reads and scored against the same long-read and assembly truth as the
pipeline's known-allele realignment (`bin/brca12_insertion_signals.py`).

- ExpansionHunter 5.0.0.
- Locus (`eh_catalog.json`): `LocusStructure` `(unit)*` with the 126 bp unit of chr17:43128247-43128372,
  `ReferenceRegion` chr17:43128246-43128372 (one copy in GRCh38). A carrier has two copies on one
  haplotype (genotype 1/2), a homozygous carrier 2/2.
- Input: the original 1000 Genomes 30x alignments (NYGC, GRCh38) of chr17:43127600-43128900, cut
  from the CRAM files (`eh_run.sh SAMPLES.tsv REFERENCE.fa OUTDIR [TRUTH.tsv]`).
- Truth: HPRC assemblies for the 39 HPRC samples (`../../hprc/results/tandem_copy_truth.tsv`) and
  1KG-ONT long reads for the independent set (`../results/tandem_copy_independent_truth.tsv`).

## Results

| Data set | Carriers called | Non-carriers called |
|---|---:|---:|
| HPRC, 39 samples | 14 of 19 | 0 of 20 |
| Independent 1KG-ONT set, 113 of 120 samples | 54 of 64 | 0 of 49 |

The independent set covers the 113 of 120 samples with a run accession in the 1000 Genomes 2504
sequence index; HG03373, NA12376, HG00418, HG00702, NA10852, HG04132 and HG03794 were not run.
ExpansionHunter returned no 2/2 genotype, although 18 carriers of the independent set carry two
copies of the linked haplotype (homozygous in the long reads).

Per-sample calls: `tandem126_hprc39_eh_vs_realign.tsv`, `tandem126_ont_independent113_eh_vs_realign.tsv`
(columns: sample, truth, genotype, confidence interval, spanning and in-repeat read counts as
reported by ExpansionHunter).

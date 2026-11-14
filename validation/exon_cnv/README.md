# Exon copy number against a reference cohort (prototype)

Research use only. Not a clinical validation. All events below are simulated (read pairs of one
pseudo-haplotype removed or duplicated inside the event), except where stated.

## Method (`cnvmodel.py`)

- Per sample and exon: mean depth over exon ± 200 bp (MAPQ ≥ 20, no duplicates), divided by the
  median depth of eight single-copy control windows on chr13/chr17 (`fetch.py`), log2.
- Reference panel: 80 unrelated 1000 Genomes 30x samples (`panel.tsv`): per exon median and
  MAD-based spread.
- Hidden Markov model over the exons of each gene in genomic order, states CN 0-4, Gaussian
  emissions with the exon's panel spread; segments with CN ≠ 2 and quality ≥ 10
  (quality = −10 log10 of the posterior probability of CN 2) are reported.
- Comparator: the current pipeline rule (exon depth / chromosome median outside 0.65-1.35).

## Design

Samples were drawn from the 2,504-sample 1000 Genomes index with seed 20261114: panel (80),
tuning (20), test (60). Parameters were chosen on the tuning set only and frozen in
`frozen_rules.json` (commit 500b3ac, pushed before any test data were fetched). Test command:
`evaluate.py ../../assets/brca12_exons.bed DATA panel.tsv test.tsv 2 --pad 200 --p-change 0.05 --min-q 10`.
Per test sample and gene: deletions of 1, 2, 3 and 5 exons and of the whole gene, duplications of
1, 2 and 3 exons, and the BRCA1 exon 1-2 deletion and duplication through the pseudogene BRCA1P1
(chr17:43123829-43167978).

## Results on the frozen test set (60 samples)

| Simulated event | Cohort model | Current rule |
|---|---:|---:|
| Deletion, 1 exon | 106/120 | 110/120 |
| Deletion, 2 exons | 115/120 | 117/120 |
| Deletion, 3 exons | 112/120 | 120/120 |
| Deletion, 5 exons | 118/120 | 120/120 |
| Deletion, whole gene | 118/120 | 120/120 |
| Duplication, 1 exon | 37/120 | 86/120 |
| Duplication, 2 exons | 81/120 | 114/120 |
| Duplication, 3 exons | 91/120 | 118/120 |
| BRCA1 exon 1-2 deletion via BRCA1P1 | 60/60 | 60/60 |
| BRCA1 exon 1-2 duplication via BRCA1P1 | 51/60 | 60/60 |
| **Calls in the 60 unmodified genomes** | **0** | **87 hints** |
| Calls outside the event in event samples | 0 | not counted |

The cohort model trades sensitivity for duplications against specificity: the current rule
reports 1.45 hints per unmodified genome, the cohort model none. Single-exon duplications
(depth ×1.5 over a few hundred bases) are near the limit of depth at 30x.

Tuning set (20 samples, `results/tuning_summary.txt`): 305 of 360 simulated events, no false call;
current rule 342 of 360 with 28 hints in the 20 unmodified genomes.

## Real event found during tuning

NA18949 (1000 Genomes, JPT, tuning set) carries a heterozygous deletion chr17:43069973-43075601
(5,629 bp, BRCA1 MANE exons 14 and 15): depth halved over the interval, 10 read pairs with
~6 kb insert spanning it, and soft-clipped reads at 43075602 whose clipped part ends at 43069972.
It is not in the 1000 Genomes phased SV panel; no long reads are available for this sample. It was
found by the cohort model as a deletion of exons 14-15 (quality 99) and is excluded from the
false-call counts. Research use only.

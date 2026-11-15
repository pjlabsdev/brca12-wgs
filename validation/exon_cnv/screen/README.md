# Real-data screen: exon copy number in 340 unmodified 1000 Genomes genomes

Research use only. Not a clinical validation.

- Samples: `screen1.tsv`, positions 220-559 of the 2,504-sample 1000 Genomes index shuffled with seed
  20261114 (disjoint from panel, tuning and both test sets), 26 populations.
- Model: exon copy number with the panel asset and default rules of main a98a4b2 (call Q >= 10,
  hint 5 <= Q < 10), depth only (no allele balance; the screen has no variant calls).
- Depth: `screen_fetch.py` (exon +/- 200 bp, reads overlapping the window, MAPQ >= 20, no duplicates;
  8 control windows), from the original 1000 Genomes CRAM alignments.
- Read-level check: `screen_eval.py` by the rule in `CONFIRMATION_RULE.md`, written before the screen
  (SHA-256 of the rule file as stored on 2026-11-14: 447a2565d7ba177e1703f9f7c87919110251ff58fd7495c07b37ff60adbc045c).
  The duplication arm of that rule counts outward pairs anywhere in a 40 kb window and is therefore
  permissive; `posthoc_local_outward_pairs` (outward pairs at the two ends of the segment) was added
  afterwards as a descriptive measure only.

## Result

340 genomes, 3 segments (`screen1_results.tsv`), none confirmed:

| Sample | Tier | Gene | Event | Exons | Quality | Depth ratio (reads) | Status |
|---|---|---|---|---|---:|---:|---|
| HG01092 | call | BRCA2 | deletion | exon 26 | 13 | 0.82 | not confirmed |
| HG01917 | hint | BRCA1 | duplication | exon 8 (legacy 9) | 5 | 1.24 | not confirmed |
| HG02139 | hint | BRCA1 | deletion | exon 3 | 5 | 0.81 | not confirmed |

False calls in real genomes: 1 call in 340 (0.3 %) and 2 hints in 340. No real exon copy-number
change was found in these 340 genomes. The one real event known so far, the heterozygous BRCA1
exon 14-15 deletion in NA18949 (tuning set, chr17:43069973-43075601), is also called by Manta 1.6.0
in a full pipeline run (PASS, 13 read pairs, 7 split reads; breakpoints chr17:43069972-43075601).

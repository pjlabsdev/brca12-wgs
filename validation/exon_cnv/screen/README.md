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

## Coordinator cross-check

Independent check by the project coordinator on the 1000 Genomes short reads and, where available,
1KG-ONT Vienna long reads (MAPQ >= 20). Positive control NA18949 exon 14-15 deletion: 26 long-insert
pairs, 5 edge clips, depth 0.59.

| Sample | Event | Short reads | ONT long reads | Verdict |
|---|---|---|---|---|
| HG01917 | BRCA1 exon 8 duplication (hint) | depth 1.54, 0 outward pairs, 0 clips | 16 spanning reads, 0 with an event | refuted |
| HG02139 | BRCA1 exon 3 deletion (hint) | depth 0.53, 0 pairs, 0 clips | 21 spanning reads, 0 with an event | refuted |
| HG01092 | BRCA2 exon 26 deletion (call) | depth 0.70, 0 long-insert pairs, 0 edge clips | no data | unconfirmed, probably false |
| NA18949 | BRCA1 exon 19 duplication (full pipeline run, Q20) | depth 1.52 over about 400 bp (other genomes 0.87-1.16), 0 outward pairs, 1 clip | no data | not a tandem duplication, unconfirmed |

Two of the three screen hits are refuted by long reads; none is supported by breakpoint evidence.
Short-window depth deviations without breakpoint evidence account for all false hits seen so far.

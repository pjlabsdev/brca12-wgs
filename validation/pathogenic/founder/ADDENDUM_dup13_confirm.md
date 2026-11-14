# Addendum: confirmation test of the BRCA1 exon 13 duplication (written before the test was run)

Research use only. Written on 2026-11-14 after the frozen test of `ADDENDUM_dup13.md` and before any
of the confirmation reads below were simulated or the real genomes were analysed with the changed
method.

## Why

The frozen test found 9 of 30 simulated carriers (`results/dup13_test_summary.txt`). Cause, found on
the frozen simulated carriers and on SPIKE_NEG, SPIKE_L2 and SPIKE_L3: (1) BWA places many reads of
the new junction at the start of the duplicated segment, outside the read window around its end;
(2) the 80 bases around both ends are Alu sequence that differs at only 6 positions, so few reads are
informative. The frozen result stands as reported; this is a changed method with a new test.

## Changed method (this commit)

- Realignment for a known duplication also takes reads starting 400 bp before to 300 bp after the
  segment start, and both haplotypes span the whole segment plus 700 bp flanks, so every reference
  junction is in both and only the new junction differs.
- New evidence: read pairs facing outward across the new junction (reverse read from 50 bp before
  to 700 bp after the segment start, forward mate from 850 bp before the segment end to the end,
  mapping quality at least 20).
- A duplication is called with 2 or more junction reads, or with 1 junction read and 3 or more
  outward pairs. Insertions and deletions keep the frozen rule.

## Test sets

1. 20 new simulated carriers (SIMULATED reads; `events_dup13_confirm.tsv`, spike.py seeds 11 to 30,
   `run_dup13_confirm.sh`), at all reads, 1/2 and 1/3.
2. The 9 spike-in scenarios as simulated non-carriers (SPIKE_NEG, SPIKE_L2 and SPIKE_L3 were looked
   at while finding the cause and are reported separately as well).
3. The 200 real samples of `../../longread/frozen/samples.tsv`, re-analysed with the changed method
   (`realneg_dup13.sh` into a new output directory).

Endpoints and evaluation as in `ADDENDUM_dup13.md` (`eval_dup13.py OUTDIR dup13_confirm`).

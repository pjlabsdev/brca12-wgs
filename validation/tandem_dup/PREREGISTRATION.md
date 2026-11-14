# Frozen test: tandem duplications from outward-facing read pairs (written before the test was run)

Research use only. Written on 2026-11-14 before any of the test reads below were simulated, fetched
or analysed with the new rule. The rule was tuned only on SPIKE_L3 and SPIKE_L5 of the earlier
spike-in test (reported separately below, not as test carriers), on the known-allele duplications
PSIDUP and DUP13C, and on 20 HPRC genomes of `../hprc/hprc_samples_1000g.tsv` that are not in the
test set. The method and the evaluation may not change after this commit; any later change is
reported as a post-hoc analysis next to the frozen result.

## Why

The read-depth exon CNV call (`../exon_cnv/`) finds most deletions but only about 6 in 10
duplications, and the allele realignment only finds duplications that are in the catalog. A tandem
duplication of any size leaves read pairs that face outward: the reverse read lies at the start of
the duplicated segment and its forward mate at the end. Short reads that cross the new junction are
clipped and carry a supplementary alignment to the other end. Neither signal depends on read depth.

## Method under test

`outward_pairs` and `tandem_duplications` in `bin/brca12_insertion_signals.py` at this commit:

- reads in the BRCA1/BRCA2 regions with mapping quality 20 or more, primary, not duplicates;
- an outward pair: the read is reverse, its mate forward on the same chromosome and further right;
- pairs are clustered when both ends lie within 600 bp of the cluster ends;
- clipped reads whose SA tag joins the two cluster ends are added to the cluster (split reads);
- call when outward pairs plus split reads (distinct read names) are 4 or more and the span is
  1,000 to 200,000 bp; a call within 600 bp of a catalog duplication is reported once, as the
  catalog allele. Breakpoints are approximate (read-pair span).

## Test sets

1. Simulated carriers (SIMULATED reads): 40 heterozygous tandem duplications from `design.py`
   (seed 202611142, `events.tsv`), each containing at least one coding exon, 22 in BRCA1 and 18 in
   BRCA2, sizes log-uniform 1 to 20 kb, random breakpoints. Each one is spiked into real HG002 read
   positions, fragment sizes and base qualities with `../pathogenic/spike.py simulate`, aligned with
   BWA-MEM, and analysed with all reads (about 30x), 1/2 and 1/3 of the reads (`run.sh`).
2. The same reads analysed with Manta 1.6 (call regions: the gene regions plus 20 kb). A Manta call
   counts when it is PASS and one breakpoint lies within 1 kb of a true breakpoint.
3. Simulated non-carriers: SPIKE_NEG, S1, S2, L1, L2, L4, L6 of the earlier spike-in test at the
   three depths. SPIKE_L3 and SPIKE_L5 (tandem duplications seen in tuning) are reported separately.
4. Real genomes: the 200 samples of `../longread/frozen/samples.tsv` (1000 Genomes 30x short reads of
   both gene regions, `realneg.sh`). No tandem duplication in BRCA1/2 is known in these samples.
   Every call is checked against the sample's ONT long reads where available and reported either
   way, as a true finding or a false call.

## Endpoints

- Sensitivity per depth with 95 % Clopper-Pearson intervals: a call matches when both ends are
  within 1 kb of the true breakpoints. Also by size (1-3 kb, 3-8 kb, 8-20 kb), for us and Manta.
- False calls: calls that match no simulated event, in simulated carriers and non-carriers, and all
  calls in real genomes that the long reads do not confirm.
- Success criterion, fixed now: at least 90 % sensitivity with all reads and no more than one false
  call per 50 real genomes. If it is not met, the rule stays out of the report and the result is
  published as a negative result.

# Tandem duplications from outward-facing read pairs

Research use only. Simulated carriers are SIMULATED reads (real HG002 read positions, fragment sizes
and base qualities), not patient samples. The rule and the evaluation were frozen and pushed in
`PREREGISTRATION.md` (commit 604c1e2) before the test reads were simulated or fetched.

## Result (`results/tandem_dup_summary.txt`, `results/tandem_dup_test.tsv`)

| Set | Read pairs (this scan) | Manta 1.6 |
|---|---|---|
| 40 simulated exon-spanning duplications, 1-20 kb, all reads (about 30x) | 40/40 | 40/40 |
| same, 1/2 of the reads | 40/40 | 38/40 |
| same, 1/3 of the reads | 39/40 | 31/40 |
| all depths, 1-3 kb / 3-8 kb / 8-20 kb | 47/48, 24/24, 48/48 | 42/48, 23/24, 44/48 |
| false calls in 147 simulated analyses (carriers and 7 non-carrier scenarios) | 0 | not counted |
| real 1000 Genomes genomes with a call | 1 of 200 (not confirmed by ONT) | not run |

Sensitivity over all depths: 119 of 120 (95 % CI 95.4-100 %), Manta 109 of 120 (84.2-95.3 %). The
one miss (TDUP_34, 2,683 bp, 1/3 of the reads) had 3 supporting reads, one below the threshold.
SPIKE_L3 and SPIKE_L5, the two duplications of the earlier spike-in test, were seen during tuning
and are reported separately: 3 of 3 and 2 of 3 depths.

The only call in real genomes, HG03802 chr17:43111922-43114097 (16 outward pairs, 3 split reads),
was checked with `ont_check.py` against the sample's 1000 Genomes ONT reads: 18 long reads span the
segment plus 500 bp on each side, none holds an insertion or split alignment of about 2.2 kb
(`results/ont_check.tsv`). For a heterozygous duplication about half of them would. It counts as a
false call: 1 in 200 real genomes, within the preregistered limit of 1 per 50. The success
criterion (at least 90 % with all reads, at most one false call per 50 real genomes) is met.

## Limits

- Simulated duplications are clean tandem copies with random breakpoints. Real duplications often
  join two Alu elements, where reads across the junction can map ambiguously and fewer pairs pass the
  mapping-quality filter; the BRCA1 exon 13 founder duplication is such a case (see
  `../pathogenic/founder/`). No real carrier genome was available.
- Breakpoints are approximate (read-pair span); zygosity cannot be read from the reads.
- Duplications below 1 kb, inverted or non-tandem (inserted elsewhere) duplications are not covered.
- Results are hints, not calls; confirmation with MLPA or another laboratory method is needed.

## Reproduce

```
THREADS=4 ./run.sh SPIKEDIR REFERENCE OUTDIR      # SPIKEDIR from ../pathogenic/run_spikeins.sh
./realneg.sh REFERENCE OUTDIR
./eval.py OUTDIR
./ont_check.py HG03802 chr17 43111922 43114097 REFERENCE
```

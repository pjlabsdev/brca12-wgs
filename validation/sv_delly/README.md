# Delly 1.3.3 instead of Manta: frozen head-to-head

Research use only. Simulated carriers are SIMULATED reads (real HG002 read positions, fragment sizes
and base qualities), not patient samples. Rules and evaluation were frozen in `PREREGISTRATION.md`
(commit b2a359c) and `ADDENDUM_filter.md` (10d7d9e, pipeline window filter, written before any output was
opened) before the run. Full numbers: `results/sv_summary.txt`, `results/sv_test.tsv`, `results/sv_real.tsv`.

## Simulated (new test set plus the tandem duplications of `../tandem_dup`)

| Events (3 depths each) | Delly | Manta 1.6 | read-pair research arm |
|---|---|---|---|
| 40 deletions, 50 bp-20 kb | 108/120 | 105/120 | 105/120 |
| ... 50-299 bp | 31/39 | 31/39 | 29/39 |
| ... 300-999 bp | 11/15 | 10/15 | 10/15 |
| ... 1-20 kb | 66/66 | 64/66 | 66/66 |
| ... with 1/3 of the reads | 30/40 | 27/40 | 28/40 |
| 20 inversions, 1-20 kb | 60/60 | 60/60 | 60/60 |
| 40 tandem duplications, 1-20 kb | 118/120 | 109/120 | 119/120 |
| ... with 1/3 of the reads | 38/40 | 31/40 | 39/40 |
| false calls, 309 simulated analyses | 0 (21 without the window filter) | 0 | 0 |

The nine earlier spike-in scenarios (partly seen in tuning): deletions 25/27 for Delly and Manta, tandem
duplications 6/6 for both. Delly does not report the Alu insertion of SPIKE_L1 (0/3, Manta 2/3); in the
pipeline insertions stay with the read-end signature and the allele catalog, which found it.

## Real genomes (1000 Genomes 30x, truth from 1KG-ONT long reads)

- 200 frozen genomes: 9 long-read deletions of 50 bp or more. Delly found the two copies of the 980 bp BRCA1
  intronic deletion (HG03439, NA19819). It missed four 51 bp deletions in the VNTR at chr17:43100607, a
  55 bp, a 403 bp and a 120 bp deletion; for the last one (HG03802) it reported a 2.1 kb deletion
  chr17:43112256-43114343 instead, which the long reads do not support (unconfirmed call; the same sample
  also gave the only false call of the tandem-duplication scan at this locus). One unconfirmed call in
  200 genomes, within the preregistered limit of 1 per 50. Without the window filter, 86 samples had a
  call, almost all a common 31 bp deletion at chr13:32386567 that HaplotypeCaller reports as a small indel.
- Real carrier list (14): all 9 long-read-confirmed carriers of the 980 bp deletion were found. Delly also
  called it in HG01323, HG01393 and HG01395 and called the 271 bp deletion in HG02602; these four are
  1000 Genomes SV-panel carriers without long reads, so the calls agree with the panel but count as
  unconfirmed in the frozen evaluation. HG01796, a panel carrier of the 271 bp deletion that the long reads
  did not confirm, has no call.

## Verdict

The preregistered criterion is met: Delly is at least as sensitive as Manta for deletions (108 vs 105) and
inversions (60 vs 60), better for tandem duplications at low depth, with no more false calls (0 vs 0) and
1 unconfirmed call in 200 real genomes. Delly becomes the default SV caller; Manta stays selectable
(`--sv_caller manta`). The read-pair research arm is as good on simulated data but stays off: it adds no
event type that Delly misses.

## Limits

- Simulated events are clean, with random breakpoints; real rearrangements between Alu elements can be
  harder to call for both tools. No real carrier of a pathogenic exon rearrangement was available.
- Deletions below about 300 bp and inside tandem repeats are called less reliably; small indels stay
  with HaplotypeCaller.
- The cloud test ran the static Linux build of Delly v1.3.3; the pipeline on macOS uses the bioconda
  osx-arm64 build of the same version.
- Harness fix after the first run (`ADDENDUM_filter.md`): the real-genome BAM headers had to be restricted
  to chr13/chr17. Method and evaluation unchanged.

## Reproduce

```
THREADS=4 DELLY=delly ./run.sh SPIKEDIR TDUPDIR REFERENCE OUTDIR   # SPIKEDIR from ../pathogenic, TDUPDIR from ../tandem_dup
DELLY=delly ./realneg.sh REFERENCE OUTDIR
./eval.py OUTDIR
```

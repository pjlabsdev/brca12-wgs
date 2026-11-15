# Frozen test: Delly 1.3.3 as the structural-variant caller instead of Manta (written before the test was run)

Research use only. Written on 2026-11-15 before any of the test reads below were simulated, fetched or
analysed with Delly. Before this commit Delly was run only on tuning material: SPIKE_L1, L2, L4, L6 and
SPIKE_NEG of the earlier spike-in test and two simulated tuning inversions (not in `events.tsv`). The
read-pair research arm was tuned on the same material plus 20 HPRC genomes of `../hprc/hprc_samples_1000g.tsv`
and the real deletion carriers HG02322, HG00637 and HG02601, which are not in the test sets. The method and
the evaluation may not change after this commit; any later change is reported as a post-hoc analysis next
to the frozen result.

## Why

Manta 1.6 (2019, no longer developed) has no working build for Apple Silicon; on macOS it runs only in
Docker. Delly 1.3.3 is maintained and runs natively. The pipeline switches to Delly if it is at least as
good as Manta on the same reads.

## Methods under test

- **Delly** (pipeline rule): `delly call -g REF -x EXCLUDE -o calls.bcf BAM` with default parameters, where
  EXCLUDE lists every reference contig except chr13 and chr17 (empty for the chr13/chr17 test reference),
  then `bcftools view -f PASS -i 'GT!="0/0" && GT!="./." && (SVTYPE=="DEL" || SVTYPE=="DUP" ||
  SVTYPE=="INV" || SVTYPE=="BND")'`. A record is reported when it overlaps a gene locus
  (`assets/brca12_regions.bed`). IMPRECISE records are kept. Version: static Linux binary of the v1.3.3
  release (sha256 5aa7f16b7a7234cde7a58884e9655967e24f5de857e226ab7d52beded3212145).
- **Manta 1.6** as in the pipeline so far (PASS records, call regions = gene loci ±20 kb), same locus filter.
- **Read-pair research arm** (`bin/brca12_insertion_signals.py --read-pair-svs`, off in the pipeline):
  deletions from forward/reverse pairs with insert above median + 6 robust spreads and same-strand split
  reads; inversions from same-strand pairs and opposite-strand split reads; at least 4 distinct reads.

## Test sets

1. New simulated carriers (SIMULATED reads, `design.py`, seed 202611151, `events.tsv`): 40 heterozygous
   deletions, size log-uniform 50 bp to 20 kb, each removing at least part of a coding exon; 20 heterozygous
   inversions, 1 to 20 kb, each containing a coding exon. Spiked into real HG002 read positions, fragment
   sizes and base qualities with `../pathogenic/spike.py` (inversion support added in this commit),
   BWA-MEM, analysed with all reads (about 30x), 1/2 and 1/3.
2. The 40 simulated tandem duplications of `../tandem_dup` at the same three depths.
3. The nine scenarios of the earlier spike-in test at three depths. Some were seen in tuning, so they are
   reported separately; SPIKE_NEG, S1 and S2 (no large event) count towards false calls.
4. Real genomes, 1000 Genomes 30x short reads of both gene regions: the 200 samples of
   `../longread/frozen/samples.tsv`, and the 14 samples of `real_carriers.tsv` (8 ONT-confirmed carriers of
   the BRCA1 980 bp intronic deletion, 4 further 1000 Genomes SV-panel carriers of it, and 2 panel carriers
   of the 271 bp BRCA1 deletion, one of which, HG01796, the long reads did not confirm). Truth: the
   long-read events of 50 bp or more of `../longread/results/all_disc.*` (1KG-ONT). Manta is not run on the
   real genomes; its real-data behaviour is in `../pathogenic/` and `../longread/`.

## Matching and endpoints

- A call matches a true event of the same type (an inversion or duplication may also be matched by a BND
  pair) when both breakpoints are within max(20 bp, min(1 kb, half the event size)).
- Sensitivity per type and depth, with 95 % Clopper-Pearson intervals; deletions also by size
  (50-299 bp, 300-999 bp, 1-20 kb).
- False calls: reported records that match no simulated event (sets 1, 2 and SPIKE_NEG/S1/S2), and in real
  genomes records that match no long-read event of the same sample (reported as unconfirmed).
- Success criterion, fixed now: Delly is at least as sensitive as Manta for deletions and for inversions
  over all depths (point estimate, not lower by more than 2 events), and has no more false calls than Manta
  in the simulated sets plus at most one unconfirmed call per 50 real genomes. If not, Manta stays the
  default where it runs, and the result is published either way.

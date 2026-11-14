# Pathogenic founder rearrangements in the known-allele catalog

Research use only. Every test here uses SIMULATED carriers (spike-ins into real HG002 read
positions, fragment sizes and base qualities, `../spike.py`) because no public short-read genome of
a carrier is available; the non-carriers include 200 real 1000 Genomes genomes. Each test was
written down and pushed before it was run (`PREREGISTRATION.md`, `ADDENDUM_*.md`).

Only alleles with published base-pair breakpoints that pass a plausibility check against GRCh38
are in the catalog (`bin/brca12_insertion_signals.py`, `KNOWN_ALLELES`, marked in `PATHOGENIC`):

| Allele | Published as | GRCh38 | Check |
|---|---|---|---|
| BRCA2 c.156_157insAlu (Portugal) | AluYa5, same orientation as BRCA2, 14 bp TSD `AAGAATCTGAACAT` (Teugels et al. 2005) | after chr13:32319165 | TSD matches the reference ending at c.156 |
| BRCA1 exon 13 duplication | c.4186-1787_4357+4122dup (Edwinsdotter Ardnor et al. 2019) | chr17:43078282-43084362, 6,081 bp | ends in two antisense Alus, 6 differences in 80 bp |
| BRCA1 exon 13 deletion (Netherlands) | IVS12-1643del3835 (Zeegers et al. 2004) = c.4186-1643_4357+2020del | chr17:43080384-43084218 | ends in two antisense Alus |

Not added: the Dutch 510 bp exon 22 deletion (`IVS21-36del510`); one end lies in no Alu element
and the ends share no homology, which does not fit the published Alu-mediated mechanism. The exact
element of c.156_157insAlu is not published; the catalog uses the Dfam AluYa5 consensus
(`AluYa5_Dfam.fa`).

## Results

Exact 95 % intervals in the summary files. "Reads 1/3" is about 10x.

| Allele | Simulated carriers found | Manta, same reads | Simulated non-carriers | Real non-carriers |
|---|---|---|---|---|
| c.156_157insAlu | 89/90 (30 elements 0 to 10 % from the consensus, 3 depths) | 5/90 | 0/24 | 0/200 |
| Exon 13 deletion 3,835 bp | 29/30 | 27/30 | 0/27 | 0/200 |
| Exon 13 duplication, frozen | 9/30 | 13/30 | 0/27 | 0/200 |
| Exon 13 duplication, confirmation (changed method, new seeds) | 58/60 | 37/60 | 0/27 | 0/200 |

- c.156_157insAlu: the read-end Alu signature alone found 28/90 (28/30 at full depth, none at 1/2
  and 1/3). The zygosity was right in all 89 calls. SPIKE_L1, whose element is a genomic AluY copy,
  was found at all three depths.
- Exon 13 deletion: SPIKE_L2 carries a deletion of the same size 116 bp away and was not called as
  this allele. The zygosity was right in all 29 calls.
- Exon 13 duplication: the frozen method missed most carriers because BWA places junction reads at
  the segment start, outside the read window, and the Alu homology at both ends leaves few
  informative bases. The changed method (`ADDENDUM_dup13_confirm.md`) reads both ends, lets both
  haplotypes span the segment and counts outward-facing read pairs; the frozen result stands as
  reported. A tandem duplication keeps every reference junction, so the zygosity is reported as
  unknown. The 2 Manta calls in non-carriers are its correct deletion call in SPIKE_L2.

## Limits

- No real carrier was available; all sensitivities come from simulated reads, aligned to a
  reference of chr13 and chr17 only (as in `../README.md`, experiment 2). The real non-carriers come
  from the original whole-genome alignments.
- The simulated elements of c.156_157insAlu vary in sequence and poly-A length, not in TSD length or
  position; the robustness test on the common intronic Alu showed that the insertion position must be
  right within about 5 bp.
- Breakpoints inside a long Alu homology are only defined up to that homology; the catalog uses the
  published positions.
- Manta ran on the same window BAMs with a 40 kb call region; it may do better on a whole genome.

## Files

- `design.py`, `events*.tsv`, `elements/`: simulated carriers.
- `run.sh`, `run_dup13.sh`, `run_dup13_confirm.sh`, `run_del3835.sh`: simulation, three depths,
  catalog realignment, Manta. `realneg*.sh`: the 200 real genomes.
- `eval.py`, `eval_dup13.py`, `eval_del3835.py`: evaluation; `results/`: per-test tables and summaries.

Reproduce: `run.sh SPIKEDIR REFERENCE OUTDIR` with SPIKEDIR from `../run_spikeins.sh`, then
`realneg.sh REFERENCE OUTDIR` and `eval.py OUTDIR` (likewise for the other alleles).

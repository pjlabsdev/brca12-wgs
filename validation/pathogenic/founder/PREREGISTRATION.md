# Frozen test: BRCA2 c.156_157insAlu (written before the test was run)

Research use only. Written on 2026-11-14 before any of the reads below were simulated, fetched or
analysed with the new catalog entry. The method and the evaluation may not change after this
commit; any later change is reported as a post-hoc analysis next to the frozen result.

## The allele

BRCA2 c.156_157insAlu (NM_000059.4), the Portuguese founder mutation (Teugels et al. 2005, Hum Mutat,
Mutation in Brief 837; Machado et al. 2007; Peixoto et al. 2009). Published structure (Teugels et al. 2005, Fig. 2):
an AluYa5 element in the same orientation as BRCA2, inserted between c.156 and c.157 in exon 3 and
flanked by a 14 bp target site duplication `AAGAATCTGAACAT`. On GRCh38 the duplicated bases are
chr13:32319152-32319165 (checked against the reference; c.156 = chr13:32319165, c.157 =
chr13:32319166, from `../cds.py`). The base-by-base sequence of the element is not published, so
the catalog carries the Dfam AluYa5 consensus with its poly-A tail (`AluYa5_Dfam.fa`) followed by the
second copy of the target site duplication. Catalog id `alu_brca2_c156`, marked pathogenic.

## Method under test

`bin/brca12_insertion_signals.py --gene-fasta` at this commit, unchanged rules of
`../../longread/frozen/PREREGISTRATION.md` (reads starting 400 bp before to 300 bp after the site,
700 bp flanks, a read supports a haplotype with at least 2 fewer edits, reads above 8 % edits
ignored, call at 2 or more allele-supporting reads; heterozygous unless reference reads are at most
10 % of the informative reads).

## Test sets

1. Simulated carriers (SIMULATED reads): 30 heterozygous carriers from `design.py` (seed 20261115):
   element divergence from the consensus 0, 2, 5, 8 and 10 % (80 % substitutions, 20 % 1 bp
   indels) times poly-A tail 15 or 30 times 3 replicates. Each element is spiked into real HG002
   read positions, fragment sizes and base qualities with `../spike.py simulate --alu-sequence`
   (14 bp target site duplication as published), aligned with BWA-MEM.
2. SPIKE_L1 of the earlier spike-in test (SIMULATED; element = a genomic AluY copy, chr17:43152573-
   43152885), as an extra carrier. Manta can map this element to its reference copy, so its Manta
   result is reported separately.
3. Simulated non-carriers: the other 8 spike-in scenarios (SPIKE_NEG, S1, S2, L2 to L6).
4. Real non-carriers: the 200 samples of `../../longread/frozen/samples.tsv` (1000 Genomes 30x short
   reads of chr13:32315077-32325000). No carrier of the founder allele is known in 1000 Genomes;
   any call is checked against the sample's ONT long reads and reported either way.

Sets 1 to 3 are analysed with all reads (about 30x) and down-sampled to 1/2 and 1/3
(`samtools view -s`, seeds 9 and 10).

## Endpoints

- Primary: sensitivity of the catalog call (pathogenic known allele within 30 bp of
  chr13:32319165) on set 1, all depths pooled and per depth, with exact 95 % intervals; false
  calls on sets 3 and 4.
- Secondary: the same for the read-end Alu signature (the earlier method, same script) and for
  Manta (any PASS record within 50 bp, sets 1 to 3); sensitivity by element divergence; zygosity
  (truth: heterozygous).

Scripts: `design.py`, `run.sh`, `realneg.sh`, `eval.py` (all in this commit).

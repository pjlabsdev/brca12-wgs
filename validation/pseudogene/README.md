# BRCA1 exons 1 and 2 and the pseudogene BRCA1P1

Research use only. The events below are SIMULATED (spike-ins into real HG002 reads); their
breakpoints are modelled, not taken from a patient.

## Mapping

BRCA1P1 lies 30 kb upstream of BRCA1 on the same strand (chr17:43.154-43.175 Mb in GRCh38). Blocks
homologous to the BRCA1 5' region: BRCA1 chr17:43118000-43122500 with BRCA1P1 43154800-43158800, and
the exon 2 region chr17:43123500-43125000 with BRCA1P1 43167700-43169400. The sequences differ
enough that 500 bp pieces of the BRCA1 5' region map uniquely (MAPQ 60). In the real HG002 reads
of the spike-in background, exons 1 to 3 have 0 % reads below MAPQ 20. In GRCh38 the pseudogene
does not cause wrong placement of short reads in exons 1 and 2. Exon 1 has low depth for another
reason (GC content; depth ratio 0.7 in the negative control).

## Pseudogene-mediated rearrangements

The recurrent BRCA1 deletion of exons 1 and 2 arises by recombination between BRCA1 intron 2 and
the homologous part of BRCA1P1 (Puget et al. 2002). Model (`events.tsv`): crossover in a 30 bp
stretch that is identical in both copies (BRCA1 chr17:43123814, BRCA1P1 43167965). The deletion
removes chr17:43123829-43167978 (44,150 bp: BRCA1 exons 1 and 2, NBR2 and the 3' part of BRCA1P1).
The reciprocal duplication of the same segment is modelled as well. Both are run through the whole
pipeline (`run.sh`).

| Event (SIMULATED) | Manta | CNVkit (flat reference) | Exon dosage hint |
|---|---|---|---|
| Deletion of exons 1 and 2 | DEL 0/1, imprecise (no split reads; ±419 bp) | deletion, exons 1 and 2 | exon 1 0.36, exon 2 0.57 |
| Duplication of exons 1 and 2 | no call | no call (bins in the segment log2 0.66, merged into a 2 Mb segment) | exon 2 1.65 only (exon 1 1.01) |

The deletion is found by all three layers. Its junction lies inside identical sequence, so there are
no split reads and the breakpoints are imprecise. The duplication is missed by Manta and CNVkit.
Only an exon dosage hint remains, which is not trustworthy alone: the negative control shows
similar hints in BRCA1 exons 5 and 6. This is a test case for copy-number calling against a
reference cohort.

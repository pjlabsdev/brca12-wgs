# Addendum: frozen test of the BRCA1 exon 13 duplication (written before the test was run)

Research use only. Written on 2026-11-14 before any of the reads below were simulated, fetched or
analysed with this catalog entry. Same rules as `PREREGISTRATION.md`.

## The allele

BRCA1 exon 13 duplication (legacy numbering; exon 12 of NM_007294.4), a founder mutation (The BRCA1
Exon 13 Duplication Screening Group 2000, Am J Hum Genet 67:207), published with base-pair
breakpoints as c.4186-1787_4357+4122dup (Edwinsdotter Ardnor et al. 2019, Fam Cancer 18:37; the
transcript version is not stated, the exon structure of NM_007294 is the same in all versions). On
GRCh38, c.4186 = chr17:43082575 and c.4357 = chr17:43082404 are the exon boundaries (`../cds.py`),
so the duplicated segment is chr17:43078282-43084362 (6,081 bp). Plausibility check before the
test: both ends lie in antisense Alu elements, and the 80 bases around the two ends are identical
except for 6 positions, as expected for the Alu-Alu recombination described for this allele
(Puget et al. 1999, Cancer Res). Within that homology the exact crossover cannot be fixed from the literature;
the catalog uses the published position.

Catalog id `dup13_brca1` (inserted sequence = copy of the reference segment after its end, marked
pathogenic). A tandem duplication keeps every reference junction, so no short read is specific to
the reference allele; the report states the zygosity as unknown.

## Test sets

1. Simulated carriers (SIMULATED reads): 10 heterozygous carriers (`events_dup13.tsv`, spike.py seeds
   1 to 10) in the real HG002 read positions, fragment sizes and base qualities, at all reads, 1/2
   and 1/3 (`run_dup13.sh`).
2. Simulated non-carriers: the 9 spike-in scenarios of `../large_events.tsv`, at the same depths.
   SPIKE_L3 carries a different 6 kb exon 12 duplication with synthetic breakpoints
   (chr17:43079800-43085799) and must not be called as this allele.
3. Real non-carriers: the 200 samples of `../../longread/frozen/samples.tsv`, 1000 Genomes 30x
   short reads of chr17:43076000-43087000 (`realneg_dup13.sh`). Any call is checked against the
   sample's ONT long reads and reported either way.

## Endpoints

Primary: sensitivity on set 1 (pooled and per depth, exact 95 % intervals); false calls on sets 2
and 3. A call is a pathogenic known duplication whose start and end are within 50 bp of the
published ends. Secondary: Manta (any PASS record within 50 bp of either end, sets 1 and 2).
Evaluation: `eval_dup13.py`.

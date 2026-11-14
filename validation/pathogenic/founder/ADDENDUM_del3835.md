# Addendum: frozen test of the BRCA1 exon 13 deletion of 3,835 bp (written before the test was run)

Research use only. Written on 2026-11-14 before any of the reads below were simulated, fetched or
analysed with this catalog entry. Same rules as `PREREGISTRATION.md`.

## The allele

Dutch founder deletion of BRCA1 exon 13 (legacy numbering; exon 12 of NM_007294.4), 3,835 bp,
Alu-mediated (Petrij-Bosch et al. 1997, Nat Genet 17:341), listed as `IVS12-1643del3835` (Zeegers et
al. 2004, Eur J Hum Genet 12:591, Table 2), i.e. c.4186-1643_4357+2020del. On GRCh38 c.4186 =
chr17:43082575 (`../cds.py`), so the deleted segment is chr17:43080384-43084218. Plausibility check
before the test: both ends lie in antisense Alu elements (13 and 23 edits over 121 bp to the AluYa5
consensus); the 80 bases around the two ends differ at 14 positions.

Not added: the Dutch 510 bp deletion of exon 22 (`IVS21-36del510`, c.5333-36_5406+400del). Its
upstream end (chr17:43049230) lies in no Alu element and the two ends share no homology, which does
not fit the published Alu-mediated mechanism; without the junction sequence the coordinates could
not be checked, so the allele stays out of the catalog.

Catalog id `del3835_brca1`, marked pathogenic.

## Test sets

1. Simulated carriers (SIMULATED reads): 10 heterozygous carriers (`events_del3835.tsv`, spike.py
   seeds 1 to 10) in the real HG002 read positions, fragment sizes and base qualities, at all reads,
   1/2 and 1/3 (`run_del3835.sh`).
2. Simulated non-carriers: the 9 spike-in scenarios of `../large_events.tsv` at the same depths.
   SPIKE_L2 carries a 3,835 bp deletion of the same exon with synthetic breakpoints 116 bp away
   (chr17:43080500-43084334) and must not be called as this allele.
3. Real non-carriers: the 200 samples of `../../longread/frozen/samples.tsv`, 1000 Genomes 30x
   short reads of chr17:43076000-43087000 (`realneg_del3835.sh`). Any call is checked against the
   sample's ONT long reads and reported either way.

## Endpoints

Primary: sensitivity on set 1 (pooled and per depth, exact 95 % intervals); false calls on sets 2
and 3; zygosity (truth: heterozygous). A call is a pathogenic known deletion whose start and end
are within 50 bp of the published ends. Secondary: Manta (any PASS record within 50 bp of either
end, sets 1 and 2). Evaluation: `eval_del3835.py`.

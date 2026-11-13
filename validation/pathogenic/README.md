# Validation with pathogenic variants, large deletions and Alu insertions

Research use only. Nothing here is a clinical validation.

The GIAB benchmark in `validation/README.md` contains no pathogenic variant, no exon deletion and
no mobile element insertion. This folder adds three experiments that do. All were run on
2026-11-13 with the tools listed in `validation/README.md`, plus CNVkit 0.9.11 for the spike-ins.

| Experiment | Data | What it shows |
|---|---|---|
| 1. Real carriers | 10 samples of the 1000 Genomes Project with a pathogenic SNV or indel | a real pathogenic variant in real reads is called |
| 2. Spike-ins | 185 small and 12 large variants written into real HG002 reads | detection by size class, per caller, with exact truth |
| 3. Real SVs | two Alu insertions and two deletions in 1000 Genomes samples (not disease-causing) | the same kind of event in real reads |

GeT-RM reference materials for BRCA1/BRCA2 (Barker et al. 2009) are distributed as DNA. No public
short-read WGS data set of a carrier of a pathogenic exon deletion, duplication or Alu insertion
was found, which is why experiment 2 uses simulated variants.

## 1. Real carriers of pathogenic variants (1000 Genomes)

Carriers were found in the joint calls of the 3,202 high-coverage genomes (Byrska-Bishop et al.
2022): variants in BRCA1/BRCA2 that are nonsense, frameshift or canonical splice variants and
pathogenic in ClinVar, plus one nonsense variant without a ClinVar entry (`carriers_1000g.tsv`).
Reads of the NYGC 30x CRAM inside the genes ±10 kb were converted back to FASTQ and run through the
pipeline (`run_carriers.sh`).

| Sample | Gene | Variant | Called | Depth | Allele depths |
|---|---|---|---|---:|---|
| HG02164 | BRCA1 | c.34C>T, p.Gln12Ter | 0/1, GQ 99 | 47 | 23,24 |
| HG03929 | BRCA1 | c.3403C>T, p.Gln1135Ter | 0/1, GQ 99 | 26 | 14,12 |
| HG00365 | BRCA1 | c.4357+1G>A | 0/1, GQ 99 | 33 | 16,17 |
| NA19397 | BRCA2 | c.2279_2283del | 0/1, GQ 99 | 25 | 10,15 |
| HG02274 | BRCA2 | c.5146_5149del | 0/1, GQ 99 | 25 | 20,5 |
| NA19130 | BRCA2 | c.5241_5242insTA | 0/1, GQ 99 | 36 | 21,15 |
| NA19132 | BRCA2 | c.5241_5242insTA | 0/1, GQ 99 | 28 | 9,19 |
| HG01377 | BRCA2 | c.7063G>T, p.Glu2355Ter | 0/1, GQ 99 | 32 | 20,12 |
| HG00611 | BRCA2 | c.7697dup | 0/1, GQ 99 | 46 | 27,19 |
| HG00612 | BRCA2 | c.7697dup | 0/1, GQ 99 | 25 | 20,5 |

10 of 10 detected as heterozygous. NA19397 is a case where the NYGC call itself did not pass VQSR.

Concordance of all calls in these ten samples with the NYGC calls (comparison, not truth): in the
exons ±20 bp 141 of 141 SNVs (2 extra) and 19 of 21 indels (2 extra); over the whole gene loci
1,772 of 1,788 SNVs and 677 of 783 indels. The indel differences outside the exons were not
examined one by one. Files: `results/carriers_summary.{tsv,json}`.

## 2. Spike-ins: simulated variants in real reads

**These are not patient samples.** The background is real HG002 data (chr13 and chr17, about 30x).
Inside each gene ±20 kb every real read pair is replaced by a simulated pair with the same position,
fragment length, read lengths and base qualities, read from one of two haplotypes: the reference or
a copy that carries the events. Sequencing errors are drawn from the real base qualities. The
simulated reads are aligned with BWA-MEM, merged with the untouched background and run through the
pipeline with Manta and CNVkit (`spike.py`, `run_spikeins.sh`, `spike_eval.py`). HGVS positions were
checked by translating the coding sequence of both transcripts (`cds.py`).

**Small variants** (`small_events.tsv`, samples SPIKE_S1 and SPIKE_S2, all heterozygous). SPIKE_S1
contains five well-known founder variants: BRCA1 c.68_69del, c.181T>G, c.5266dup and BRCA2
c.771_775del, c.5946del.

| Class | Spiked | Found |
|---|---:|---:|
| SNV | 90 | 90 |
| Deletion 1 bp | 11 | 11 |
| Deletion 2 to 5 bp | 18 | 18 |
| Deletion 6 to 20 bp | 16 | 16 |
| Deletion 21 to 50 bp | 8 | 8 |
| Insertion or duplication 1 bp | 10 | 10 |
| Insertion or duplication 2 to 5 bp | 16 | 16 |
| Insertion or duplication 6 to 20 bp | 8 | 8 |
| Insertion or duplication 21 to 50 bp | 8 | 8 |
| Total | 185 | 185 |

No false positive call in the exons. Every sample, including the negative control, shows two calls
at chr13:32367575 and 32367579 that come from residual real HG002 reads (pairs with an unmapped
mate are not replaced).

**Large events** (`large_events.tsv`, samples SPIKE_L1 to SPIKE_L6, one event per gene, all
heterozygous, breakpoints synthetic).

| Sample | Gene | Event | Manta | CNVkit | Exon dosage hint |
|---|---|---|---|---|---|
| L1 | BRCA1 | deletion 510 bp, exon 21 | yes | no | yes |
| L1 | BRCA2 | Alu insertion in exon 3 (327 bp) | yes, as breakend pair | no | not applicable |
| L2 | BRCA1 | deletion 3,835 bp, exon 12 | yes | no | yes |
| L2 | BRCA2 | deletion 8.5 kb, exons 14 to 16 | yes | no | yes |
| L3 | BRCA1 | tandem duplication 6 kb, exon 12 | yes | no | yes |
| L3 | BRCA2 | deletion 1.5 kb, exon 21 | yes | no | yes |
| L4 | BRCA1 | deletion 11 kb, exons 1 and 2 | yes | yes | yes |
| L4 | BRCA2 | deletion 120 bp in exon 11 | yes | no | no |
| L5 | BRCA1 | deletion of the whole gene (135 kb) | yes | yes | yes |
| L5 | BRCA2 | tandem duplication 5 kb, exons 19 and 20 | yes | no | yes |
| L6 | BRCA1 | deletion 300 bp in exon 10 | yes | no | no |
| L6 | BRCA2 | deletion of the whole gene (95 kb) | yes | yes | yes |

- Manta: 12 of 12 with exact breakpoints and heterozygous genotype.
- CNVkit: 0 of 12 before the correction made in this work (`--drop-low-coverage`), 3 of 12 after
  it. Single-sample CNVkit with 1 kb bins and a flat reference does not find events below about
  10 kb. `results/spike_summary.before_fix.json` keeps the result before the correction.
- Exon dosage hint: 9 of 11 deletions and duplications. The negative control shows five hints
  without any event (BRCA1 exons 5, 6 and 21, BRCA2 exons 8 and 19), so a hint alone means little.
- The 120 bp deletion was also found by the small-variant caller.
- The read-end signature (see below) reports the spiked Alu insertion and nothing in the other
  eight spike-in samples.
- Negative control: no SV, no CNV candidate.

The spiked Alu element was copied from an Alu on chr17 of the reduced reference. Its reads can be
placed there, which makes detection easier than for a real insertion.

Files: `results/spike_summary.json`, `results/depth_bins.tsv`.

## 3. Real structural variants (1000 Genomes)

@@SV_SECTION@@

## Limits

- Reads were taken from existing alignments near the genes and realigned against chr13 and chr17
  only. There is no competition from similar sequence elsewhere in the genome, which may overstate
  precision. A whole-genome run against the complete reference is still missing.
- Spiked variants are synthetic. They show whether a variant that is present in the reads is called
  and reported, not how real pathogenic rearrangements look.
- No real sample with a pathogenic exon deletion, duplication or Alu insertion was tested.
- VEP annotation and the ACMG/AMP hints were not part of these runs (no VEP cache on the validation
  machine). In experiments 1 and 3 CNVkit was replaced by a neutral stand-in.
- The numbers are small: ten carriers, twelve large events, eleven Alu carriers.
- The read-end signature was written after the missed Alu insertions had been inspected. Its
  result on the same carriers is not an independent test.

## Reproduce

```bash
validation/run_giab.sh WORKDIR                         # reference, tools, GIAB runs
validation/pathogenic/run_carriers.sh WORKDIR          # experiment 1
CNVKIT=/path/to/cnvkit.py validation/pathogenic/run_spikeins.sh WORKDIR   # experiment 2
MATES=1 validation/pathogenic/run_sv_carriers.sh WORKDIR                  # experiment 3
```

## References

- Byrska-Bishop M et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project
  cohort including 602 trios. Cell 2022;185:3426-3440.
- Barker SD et al. Development and characterization of reference materials for MTHFR, SERPINA1,
  RET, BRCA1, and BRCA2 genetic testing. J Mol Diagn 2009;11:553-561.
- Chen X et al. Manta: rapid detection of structural variants and indels for germline and cancer
  sequencing applications. Bioinformatics 2016;32:1220-1222.
- Talevich E et al. CNVkit: genome-wide copy number detection and visualization from targeted DNA
  sequencing. PLoS Comput Biol 2016;12:e1004873.

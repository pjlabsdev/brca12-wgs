# Long-read check of structural variants (1KG-ONT)

Research use only. Not a clinical validation.

The structural-variant sites of `pathogenic/sv_sites_1000g.tsv` come from short-read callsets of
the 1000 Genomes panel. Here they are checked in Oxford Nanopore long reads of the same
individuals (1KG-ONT Vienna cohort, Schloissnig et al. 2025, GRCh38 CRAMs read remotely, only the
two gene regions).

- `ont_genotype.py`: a read counts when it spans the site with a margin on both sides; it supports
  the variant when its alignment carries an insertion (Alu) or deletion of at least half the
  expected length near the site. Present when at least 2 spanning reads and at least 10 % of
  them support it.
- `ont_discover.py`: finds every insertion and deletion of 50 bp or more in BRCA1/BRCA2 in all
  cohort samples, clusters them into sites and genotypes every sample at every site.
- `compare_shortread.py`: joins the long-read genotypes with the panel genotypes and the
  pipeline results of `pathogenic/run_sv_carriers.sh`.

## Results for the carrier samples (`results/`)

11 of the carrier and control samples of `pathogenic/sv_carriers_1000g.tsv` have long reads
(`results/ours_ont.txt`).

- The pipeline (Manta call or read-end insertion signal) agrees with the long reads in 26 of 26
  sample-site pairs (`results/ours_vs_shortread.tsv`).
- Long reads confirm the intronic BRCA2 Alu insertion in HG00117, HG00637 and HG01377 and the
  1,081 bp BRCA1 deletion (980 bp in the reads) in HG00637 and HG01990.
- The 271 bp BRCA1 deletion that the panel lists for HG01796 is absent in all 26 spanning long
  reads; the pipeline did not report it either. The panel genotype is most likely wrong.

## All 1,019 samples (`results/all_disc.*`)

`ont_discover.py` collected every insertion and deletion of 50 bp or more inside both gene loci from
the long reads of all 1,019 cohort samples (no sample failed) and genotyped every sample at every
site: 18 sites, none in an exon ±20 bp, 9 seen in a single sample.

- BRCA2 intronic Alu insertion (chr13:32348296, 316 bp): 136 of 1,002 evaluable samples.
- BRCA1 1,081 bp intronic deletion (980 bp in the reads): 13 of 993.
- BRCA1 VNTR at chr17:43100560-43100705: several alleles from -77 to +568 bp.
- 126 bp tandem-repeat insertion upstream of BRCA1 exon 1 (chr17:43128267): 552 of 1,004. The
  short-read callers miss it; the tandem-copy signature finds it in 5 of 6 long-read
  carriers among the 11 tested samples and in none of the 5 non-carriers
  (`results/tandem_copy.tsv`, `results/tandem_copy_truth.tsv`).

Sites outside the gene loci (reads reach beyond them) are dropped, and a shorter insertion cluster
with the same carriers next to a longer one (an insertion aligned in two pieces) is merged into it.

## Sources

Schloissnig S et al. Structural variation in 1,019 diverse humans based on long-read sequencing.
Nature 2025. Data generated at the Institute of Molecular Pathology (Vienna, Austria) with funds
provided by Boehringer-Ingelheim.

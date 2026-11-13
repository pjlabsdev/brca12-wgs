# Changelog

## Unreleased

- CNVkit WGS mode now gets one BED for `--targets` and `--access` (the accessible part of chr13/chr17);
  CNVkit rejects two different files there, so the CNV step failed on real runs.
- `envs/brca12.yml` no longer lists Manta, which needs Python 2.7 and made the main environment
  unsolvable; Manta keeps its own environment (`envs/manta.yml`).
- Manta breakend pairs (linked by `MATEID`) are reported once, from the end inside BRCA1/BRCA2,
  instead of twice with an off-gene coordinate.
- Genotypes with a missing allele (`0/.`, `./1`) are shown as "not determined".

## 1.0.1 (2026-11-12)

- Archive release for the Zenodo DOI (Zenodo licence id in lowercase SPDX form). No pipeline changes.

## 1.0.0 (2026-11-12)

- First public release: FASTQ or BAM/CRAM to BRCA1/BRCA2 report on GRCh38.
- SNVs/indels (GATK HaplotypeCaller or DeepVariant), SVs (Manta), CNVs (CNVkit + exon dosage).
- Ensembl VEP with MANE Select, gnomAD and ClinVar; ACMG/AMP hints with ENIGMA thresholds.
- Exon coverage with uncertain exons and low-depth intervals; HTML/Markdown report in English or German.
- Pinned conda environments, command log and tool versions for reproducibility.

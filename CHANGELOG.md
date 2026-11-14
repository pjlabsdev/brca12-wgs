# Changelog

## Unreleased (1.1.0)

- ACMG/AMP hints now follow the ClinGen ENIGMA BRCA1/BRCA2 specifications v1.2:
  - BA1, BS1 and the new BS1_Supporting use the filtering allele frequency of gnomAD non-cancer
    non-founder populations (`--gnomad_exomes_vcf`, `--gnomad_genomes_vcf`). Without these files
    the hints are labelled as approximations.
  - PM2_Supporting is no longer given to insertions, deletions or delins, and needs both gnomAD files.
  - PP3, BP4, BP1_Strong and BP7 are evaluated with BayesDel no-AF inside the functional domains
    and SpliceAI (`--bayesdel_tsv`, `--spliceai_vcf`) instead of being reported as open.
  - The blanket "BP4/BP7 candidate" for every non-coding variant is gone.
- New JSON fields per variant: `gnomadFAF`, `spliceAI`, `bayesDel`; top level `evidenceSources`.

## 1.0.2 (2026-11-13)

- GIAB HG001/HG002 benchmark of the BRCA1/BRCA2 small-variant calls (`validation/`).
- DOI badge and citation in the README, `doi:` in CITATION.cff.
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

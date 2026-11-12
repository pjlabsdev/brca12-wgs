# Changelog

## 1.0.0 (2026-11-12)

- First public release: FASTQ or BAM/CRAM to BRCA1/BRCA2 report on GRCh38.
- SNVs/indels (GATK HaplotypeCaller or DeepVariant), SVs (Manta), CNVs (CNVkit + exon dosage).
- Ensembl VEP with MANE Select, gnomAD and ClinVar; ACMG/AMP hints with ENIGMA thresholds.
- Exon coverage with uncertain exons and low-depth intervals; HTML/Markdown report in English or German.
- Pinned conda environments, command log and tool versions for reproducibility.

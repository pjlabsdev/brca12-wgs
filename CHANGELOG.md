# Changelog

## Unreleased (1.1.0)

- Exome, hybrid-capture panel and amplicon data (`--assay`, `--targets_bed`): no duplicate marking
  and no per-start downsampling for amplicons, primer clipping with `samtools ampliconclip`
  (`--primers_bed`), CNVkit hybrid/amplicon modes, optional pooled CNV reference
  (`--cnv_reference`, also for WGS), Manta exome mode, exons outside the targets marked uncertain.
- Read and target QC (`qc/read_qc.json`, report section, Excel sheets): reads, mapping, pairing,
  duplicates, insert size, error rate; per target depth, uniformity, fold-80, on-target rate,
  weakly covered targets; share of reads with a clipped primer.
- Mosaic candidates (`--mosaic`): GATK Mutect2 tumour-only with FilterMutectCalls, allele fraction
  limits, minimum supporting reads and depth; VEP-annotated; validated with spike-ins at 5 to 30 %.
- Single-end FASTQ (only `--reads_r1`; Manta skipped) and `--platform iontorrent` (indels flagged,
  mosaic candidates SNVs only).
- Container: `Dockerfile` with all tools at the pinned versions, profiles `docker` and
  `singularity`, image built and smoke-tested in CI. `perl-list-moreutils` added to the conda
  environment (missing VEP dependency).
- Validation on exome and amplicon data (`validation/targeted/`): GIAB HG002 Agilent exome, HG001
  TruSeq exome and HG002 Ion AmpliSeq exome against GIAB v4.2.1; mosaic spike-ins
  (`validation/pathogenic/run_mosaic.sh`, `mosaic_eval.py`).
- README: feature comparison with BRCA-analyzer.
- Validation against assembly-based truth for 39 HPRC samples (`validation/hprc/`) and a long-read
  check of the structural variants in the 1KG-ONT cohort (`validation/longread/`, all 1,019 samples).
- Tandem-copy signature (`brca12_insertion_signals.py --gene-fasta`): extra copies of a repeat unit
  (50 to 1,000 bp) from clipped reads at both unit borders, reported next to the Alu signature.
  Finds the common 126 bp insertion upstream of BRCA1 exon 1 in 16 of 19 HPRC and 4 of 6 ONT carriers,
  none in 25 non-carriers.
- Paper: sections on targeted data, low allele fractions, assembly truth and long reads; Figure 7.

- GRCh37 position and BIC-style legacy names for every small variant (`grch37Pos`, `bicNames`,
  shown next to HGVS), for comparison with older reports; checked against 4,464 BRCA Exchange BIC
  names (89.3 % identical, 96.9 % same position and type).
- Several samples per run: `bin/brca12_batch.py` reads a sample sheet (FASTQ pair or BAM/CRAM per
  sample), runs each sample and writes a cohort summary (`bin/brca12_cohort.py`: samples, all calls
  with carrier counts, contamination hints when a low-fraction call is clearly present in another
  sample of the batch, `cohort.xlsx`).
- Excel export of every result (`brca12_report.xlsx`, process `XLSX_EXPORT`, no extra packages).
- HTML report with figures per gene: mean depth per exon with threshold lines and an exon dosage
  plot (MLPA-like view of possible exon losses or gains), inline SVG.
- Russian report language (`--report_lang ru`; machine-assisted, not reviewed by a native speaker).
- Optional PDF report (`--report_pdf`, `--chrome`): the HTML report printed with headless
  Chrome/Chromium, A4 with header (sample, RUO) and page numbers; skipped with a note when no
  browser is available.
- ACMG/AMP hints now follow the ClinGen ENIGMA BRCA1/BRCA2 specifications v1.2:
  - BA1, BS1 and the new BS1_Supporting use the filtering allele frequency of gnomAD non-cancer
    non-founder populations (`--gnomad_exomes_vcf`, `--gnomad_genomes_vcf`). Without these files
    the hints are labelled as approximations.
  - PM2_Supporting is no longer given to insertions, deletions or delins, and needs both gnomAD files.
  - PP3, BP4, BP1_Strong and BP7 are evaluated with BayesDel no-AF inside the functional domains
    and SpliceAI (`--bayesdel_tsv`, `--spliceai_vcf`) instead of being reported as open.
  - The blanket "BP4/BP7 candidate" for every non-coding variant is gone.
- New JSON fields per variant: `gnomadFAF`, `spliceAI`, `bayesDel`; top level `evidenceSources`.
- Validation with pathogenic variants (`validation/pathogenic/`): ten real carriers from the
  1000 Genomes Project, 185 small and 12 large variants spiked into real reads, and real Alu
  insertions and deletions in 1000 Genomes samples. Results and limits are in
  `validation/pathogenic/README.md`.
- CNVkit segmentation now runs with `--drop-low-coverage`. Without it, bins without coverage
  dominated the segmentation and no spiked deletion, not even a whole-gene deletion, was reported.
- Exon dosage hints are no longer suppressed for exons that are also flagged as uncertain; a deleted
  exon is both. The hint says that loss and a coverage problem cannot be told apart.
- New: insertion signatures from soft-clipped read ends inside the genes
  (`bin/brca12_insertion_signals.py`, JSON `insertionSignals[]`, listed in the SV table). Only
  signatures whose clipped sequence contains Alu consensus k-mers are listed. Hints, not calls.
- Manta insertions whose inserted sequence carries an Alu signature, and pairs of breakends at the
  same site, are flagged as possible mobile element insertions.
- `--bam` input whose file name equals `<sample>.bam` no longer fails in the CNVkit step.
- PVS1 and PM5_PTC from the ENIGMA decision trees (Table 4) per exon, splice site and initiation
  codon instead of a generic candidate; PS3/BS3 from calibrated functional assays (Table 9); PP4/BP5
  from published combined likelihood ratios. Tables bundled in `assets/enigma` with source and
  checksums (`bin/brca12_make_enigma.py`).
- Exon deletions and duplications (CNV candidates, exon dosage hints) show the Table 4 code.
- Conda environments could not be solved: CNVkit 0.9.11 needs Python 3.10 (was 3.11), and the Manta
  environment no longer pins htslib (conflicts with Python 2.7).

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

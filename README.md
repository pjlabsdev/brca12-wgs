# brca12-wgs

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24173638.svg)](https://doi.org/10.5281/zenodo.24173638)

**BRCA1/BRCA2 germline analysis from whole-genome sequencing, GRCh38 only.**
One Nextflow pipeline from FASTQ (or an existing BAM/CRAM) to a readable report: SNVs and small
indels, copy-number changes, structural variants, ClinVar, ACMG/AMP hints and exon-level coverage
with uncertain regions flagged.

> **Research use only (RUO).** Results must not be interpreted medically without clinically
> validated confirmation (accredited laboratory, e.g. Sanger sequencing or MLPA). ACMG/AMP entries
> are automated hints on individual criteria, never a classification. An unremarkable result does
> not rule out a hereditary BRCA1/BRCA2 predisposition.

## Why another BRCA pipeline?

Older BRCA workflows such as [BRCA-analyzer](https://github.com/aakechin/BRCA-analyzer) target
amplicon panels on hg19 and call SNVs/indels only. brca12-wgs is built for WGS on GRCh38:

| | brca12-wgs | amplicon-era workflows |
|---|---|---|
| Data | WGS (FASTQ or BAM/CRAM) | amplicon panels |
| Reference | GRCh38/hg38 (checked from contig lengths) | hg19 |
| SNVs/indels | GATK 4 HaplotypeCaller or DeepVariant | GATK 3 |
| CNVs | CNVkit + per-exon dosage | none |
| Structural variants | Manta | none |
| Annotation | Ensembl VEP, MANE Select, gnomAD, your ClinVar VCF | ANNOVAR (licence required) |
| ACMG/AMP | ENIGMA-threshold hints, clearly not a classification | none |
| Coverage | per exon at 10/20/30x, MAPQ-0 fraction, low-depth intervals | per amplicon |
| Reproducibility | pinned conda envs, every command logged, tool versions in report | manual |
| Output | HTML, Markdown and machine-readable JSON | Excel, Word |

## Quick start

```bash
git clone https://github.com/<owner>/brca12-wgs && cd brca12-wgs

# From FASTQ
nextflow run main.nf -profile conda \
  --sample sample --reads_r1 sample_R1.fastq.gz --reads_r2 sample_R2.fastq.gz \
  --reference GRCh38.fa --vep_cache ~/.vep --clinvar_vcf clinvar.vcf.gz \
  --outdir brca12_results --threads 8

# From an existing BAM/CRAM (alignment is skipped)
nextflow run main.nf -profile conda --sample sample --bam sample.bam \
  --reference GRCh38.fa --vep_cache ~/.vep --clinvar_vcf clinvar.vcf.gz
```

Requirements: Nextflow ≥ 23.10, Java 17, conda/mamba (or the tools on `PATH`). The reference needs
`.fai`, `.dict` and the BWA index next to the FASTA. VEP needs an offline GRCh38 cache
(`vep_install -a cf -s homo_sapiens -y GRCh38`). ClinVar: `clinvar.vcf.gz` + `.tbi` from NCBI (GRCh38).

## What it does

| Step | Tool | Output |
|---|---|---|
| Input check (GRCh38 via chr13/chr17 lengths, BED inside the loci) | `bin/brca12_prepare.py` | `qc/reference_check.json` |
| Alignment | BWA-MEM, samtools fixmate/sort/markdup/index | `alignment/` |
| SNVs/indels | GATK HaplotypeCaller (default) or DeepVariant (`--snv_caller deepvariant`) | `calls/<sample>.small.vcf.gz` |
| Structural variants | Manta, breakpoints within ±1 Mb of the genes | `calls/<sample>.manta.diploidSV.vcf.gz` |
| CNVs | CNVkit WGS, flat reference over chr13/chr17, plus per-exon dosage | `calls/cnvkit/` |
| Per-gene extraction | bcftools norm/view with `brca12_regions.bed` | `variants/BRCA1.*`, `variants/BRCA2.*` |
| Annotation | Ensembl VEP (MANE, gnomAD, HGVS) + ClinVar (`--custom`) | `variants/*.vep.vcf.gz` |
| Exon coverage | mosdepth (MAPQ ≥ 20 and ≥ 0), samtools depth | `coverage/` |
| Report | `bin/brca12_report.py` | `brca12_report.html/.md`, `brca12_results.json`, `brca12_commands.sh` |

The report has separate sections for BRCA1 SNVs/indels, BRCA2 SNVs/indels, BRCA1/BRCA2 CNVs,
BRCA1/BRCA2 structural variants, coverage per exon and uncertain regions. Use `--report_lang de`
for a German report. The JSON (schema `genome.brca12.v1`) is meant for downstream tools; see
[docs/output.md](docs/output.md).

Gene loci (GRCh38, 1-based): BRCA1 chr17:43044295-43170327 (minus strand), BRCA2
chr13:32315077-32400268 (plus strand). Exons: MANE Select ENST00000357654 / NM_007294.4 (23 exons;
BRCA1 also shows the legacy numbering without exon 4) and ENST00000380152 / NM_000059.4 (27 exons),
generated from Ensembl with `bin/brca12_make_exons.py`.

## Flagging uncertain results

- Exon **uncertain** if mean depth < 20x, fewer than 95 % of bases ≥ 20x, or more than 10 % of reads
  with MAPQ < 20 (pseudogene BRCA1P1, repeats). Thresholds: `--min_exon_depth`, `--min_exon_pct20`,
  `--max_mapq0`.
- Every interval below `--low_depth` (10x) inside the genes is listed.
- Variant quality flags: FILTER, DP < 15, GQ < 20, allele fraction outside 25-75 % for heterozygous
  calls, MQ < 40, location in an uncertain exon.

## ACMG/AMP hints

Automated pre-checks of single criteria using ClinGen ENIGMA BRCA1/2 VCEP thresholds: PVS1
candidate (stop, frameshift, canonical splice site; restricted in the last exon), PM2_Supporting
(absent from gnomAD), BA1 (> 0.1 %), BS1 (> 0.01 %), open PP3/BP4 for missense, splice warnings,
BP4/BP7 candidates, and known exceptions such as BRCA2 p.Lys3326Ter. ClinVar significance and
review status (stars) are shown separately. **No classification is assigned.**

## Reproducibility

- Pinned versions in `envs/brca12.yml` and `envs/manta.yml` (Manta 1.6 needs Python 2.7; on Apple
  Silicon use `CONDA_SUBDIR=osx-64`).
- The report lists the tool versions that actually ran.
- `brca12_commands.sh` contains every executed command verbatim; `pipeline_info/` holds the
  Nextflow trace, timeline and report.

## Limitations

- Single sample: CNVs are candidates only (no cohort for GATK gCNV); confirm by MLPA.
- Mosaics below ~15 %, Alu insertions (BRCA2 exon 3) and complex rearrangements are not reliably
  detected by short-read WGS.
- DeepVariant on macOS is only available via Docker; the default is GATK HaplotypeCaller.

## Tests

```bash
python3 -m unittest discover -s tests   # report logic on synthetic data
bash tests/run_stub_pipeline.sh         # full Nextflow run with stand-in tools
```

Both run in GitHub Actions on every push. They check the pipeline wiring and the reporting logic,
not the accuracy of the external callers. Validation against GIAB or a characterised BRCA reference
sample is still to come.

## Deutsch (Kurzfassung)

BRCA1/BRCA2-Keimbahnauswertung aus WGS, nur GRCh38: BWA-MEM, samtools, GATK HaplotypeCaller,
Manta, CNVkit, VEP mit ClinVar und mosdepth je Exon, Bericht als HTML/Markdown/JSON
(`--report_lang de` für Deutsch). Nur für Forschung (RUO). Ergebnisse dürfen nicht ohne klinisch
validierte Bestätigung medizinisch interpretiert werden.

## Citation

Jungwirth P. brca12-wgs: BRCA1/BRCA2 germline analysis from whole-genome sequencing (GRCh38),
version 1.0.1. Zenodo, 2026. https://doi.org/10.5281/zenodo.24173638
(see also [CITATION.cff](CITATION.cff)).

## Licence

MIT, see [LICENSE](LICENSE). The pipeline calls third-party tools that keep their own licences
(GATK BSD-3, Manta GPL-3, CNVkit Apache-2.0, VEP Apache-2.0, BWA GPL-3, samtools/bcftools MIT,
mosdepth MIT); they are installed separately and not redistributed here.

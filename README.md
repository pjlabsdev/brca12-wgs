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

## Comparison with BRCA-analyzer

[BRCA-analyzer](https://github.com/aakechin/BRCA-analyzer) (Kechin et al. 2018, last commit February
2020) is the best-known open BRCA1/BRCA2 workflow. The table compares it, as described in its README
and code, with brca12-wgs 1.1.0. Feature by feature, brca12-wgs covers what BRCA-analyzer offers and
adds the rows it lacks; the one area where BRCA-analyzer is ahead is the size of its clinical
validation (last row).

| | brca12-wgs 1.1.0 | BRCA-analyzer |
|---|---|---|
| Data types | WGS, exome, hybrid-capture panels and amplicon panels (`--assay`), paired or single-end, Illumina or Ion Torrent | amplicon panels (MiSeq) |
| Input | FASTQ, BAM or CRAM | FASTQ |
| Reference | GRCh38 (checked from contig lengths); GRCh37 position given for every small variant | hg19 |
| Primer removal | `samtools ampliconclip` with a primer BED, share of clipped reads reported | cutPrimers, share of trimmed reads reported |
| Duplicates | marked (not for amplicons) | not marked |
| SNVs/indels | GATK 4 HaplotypeCaller or DeepVariant | Pisces (default), GATK 3 UnifiedGenotyper or FreeBayes, GATK 3 realignment/BQSR |
| Low allele fraction | mosaic candidates with GATK Mutect2 (`--mosaic`), validated by spike-ins | Pisces calls low fractions |
| CNVs | CNVkit (WGS, hybrid or amplicon mode, optional pooled normal reference) and per-exon dosage | none |
| Structural variants | Manta, plus a read-end signature for Alu insertions | none |
| Annotation | Ensembl VEP (MANE Select, HGVS, gnomAD, SIFT/PolyPhen), your ClinVar VCF, optional SpliceAI/BayesDel | ANNOVAR (licence required), snpEff, ClinVar 2016 bundled |
| Legacy names | BIC-style names (e.g. 185delAG, 5382insC, 6174delT) computed from the MANE CDS | BIC names from the bundled BIC database |
| ACMG/AMP | hints per ClinGen ENIGMA v1.2 criterion incl. ENIGMA tables 4 and 9, never a classification | class from ClinVar/BIC |
| Read QC | reads, mapping, pairing, duplicates, insert size, error rate | read and index counts |
| Coverage | per exon (10/20/30x, MAPQ < 20 share), per target, uniformity, fold-80, on-target rate, every interval below 10x | per amplicon, uniformity figure |
| Several samples | sample sheet, parallel runs, cohort table with carrier counts | yes, joint table |
| Contamination | hint when a low-fraction variant is clearly present in another sample of the batch | same idea |
| Figures | depth per exon and MLPA-like exon dosage plot per gene | amplicon uniformity |
| Output | HTML, Markdown, PDF, Excel and JSON (schema `genome.brca12.v1`) | Excel, Word |
| Report language | English, German, Russian (Russian not reviewed by a native speaker) | English, Russian |
| Reproducibility | Nextflow, pinned conda environments, Docker/Singularity image, every command and tool version in the report | Python script, Docker image |
| Tests | unit tests and three full pipeline runs in CI | none in the repository |
| Open validation | GIAB HG001/HG002 (WGS, 2 exomes, Ion AmpliSeq), real pathogenic carriers, spike-ins incl. mosaics, real SVs with long-read confirmation; all scripts in `validation/` | example data only |
| Clinical validation | none (research use only) | more than 900 patient samples, pathogenic variants confirmed by Sanger sequencing |

The last row cannot be matched with open data: it would need sequencing of characterised reference
DNA (e.g. Coriell) or access to a clinical cohort.

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

# Amplicon panel, single-end or paired, primers clipped, mosaic candidates, in the container
docker build -t brca12-wgs:1.1.0 .
nextflow run main.nf -profile docker --sample S1 --reads_r1 S1_R1.fastq.gz --reads_r2 S1_R2.fastq.gz \
  --assay amplicon --targets_bed amplicons.bed --primers_bed primers.bed --mosaic \
  --reference GRCh38.fa --vep_cache ~/.vep --clinvar_vcf clinvar.vcf.gz
```

Requirements: Nextflow ≥ 23.10, Java 17, and conda/mamba, Docker or Singularity (or the tools on
`PATH`). The container (`Dockerfile`, profiles `docker` and `singularity`, built in CI) holds every
tool at the pinned versions; VEP cache, reference and ClinVar are mounted from the host. The reference needs
`.fai`, `.dict` and the BWA index next to the FASTA. VEP needs an offline GRCh38 cache
(`vep_install -a cf -s homo_sapiens -y GRCh38`). ClinVar: `clinvar.vcf.gz` + `.tbi` from NCBI (GRCh38).

## Several samples (batch) and Excel

Every run also writes `brca12_report.xlsx` (overview with the research-use note, small variants with
GRCh38/GRCh37 positions, HGVS, BIC names, ClinVar, ACMG/AMP hints and quality flags, CNV, SV, exon
coverage, low-depth intervals, tool versions; no extra Python packages).

For a batch, list the samples in a sample sheet and start them together:

```bash
# samples.csv: sample,reads_r1,reads_r2,bam   (FASTQ pair or BAM/CRAM per line)
python3 bin/brca12_batch.py --samplesheet samples.csv --outdir batch --run-name RUN01 --jobs 2 -- \
  -profile conda --reference GRCh38.fa --vep_cache ~/.vep --clinvar_vcf clinvar.vcf.gz
```

Each sample runs as its own Nextflow run in `batch/<sample>/` (resumable alone with `--resume`).
`bin/brca12_cohort.py` then writes `batch/cohort/`: `cohort_samples.tsv` (one line per sample with
variant counts, uncertain exons, median exon depth), `cohort_variants.tsv` (all calls with the number
of carriers in the batch), `cohort_contamination.tsv` and `cohort.xlsx`. A **contamination hint**
is given when a variant has a low allele fraction (3 to 25 %) in one sample but is clearly
heterozygous or homozygous in another sample of the batch; the hint names the possible source. It is
a reason to check (identity SNPs, repeated library), not a finding.

## Exome, gene panels and amplicons

`--assay exome|panel|amplicon` with `--targets_bed` (the capture or amplicon intervals; contigs with
or without `chr`). Only targets that overlap BRCA1/BRCA2 are used.

- **Amplicons:** no duplicate marking; HaplotypeCaller without per-start downsampling;
  `--primers_bed` clips primers from both read ends with `samtools ampliconclip` (strand-aware when
  the BED has a strand column). Reads shorter than 30 bp after clipping are dropped.
- **CNVs:** CNVkit in hybrid (exome, panel) or amplicon mode. Without a reference from normal
  samples a single sample gives candidates only; build one with `cnvkit.py reference` from normals of
  the same assay and pass it with `--cnv_reference` (also possible for WGS). The per-exon dosage uses
  the median exon depth of both genes as baseline, so a whole-gene loss stays visible.
- **Structural variants:** Manta in exome mode; skipped for single-end reads.
- **Coverage:** every target with mean and minimum depth and share of bases ≥ 20x; targets below
  `--min_target_depth` (30x) are listed. Summary: uniformity (share of target bases with at least
  0.2 x the mean depth), fold-80 base penalty, on-target rate and, for amplicons, the share of reads
  with a clipped primer. Exons not covered by any target are marked uncertain.
- **Ion Torrent:** `--platform iontorrent` (single-end FASTQ with only `--reads_r1`). Indels are
  flagged as unreliable (homopolymer errors of the platform) and mosaic candidates are SNVs only.

## Read QC and mosaic candidates

Every run writes `qc/read_qc.json` and a report section with reads, mapping, proper pairs,
duplicates, MAPQ 0, insert size, error rate and base quality over the gene regions, with warnings
(for example more than 30 % duplicates, error rate above 1 %).

`--mosaic` adds GATK Mutect2 in tumour-only mode over the gene regions. Calls that pass
FilterMutectCalls, have an allele fraction between `--mosaic_min_af` (0.02) and `--mosaic_max_af`
(0.30), at least `--mosaic_min_alt_reads` (3) supporting reads and `--mosaic_min_depth` (20) reads,
and that the germline caller did not report are listed as **mosaic candidates**, annotated like
germline variants. They always need confirmation by an orthogonal method; detection depends on depth
(see Validation).

## What it does

| Step | Tool | Output |
|---|---|---|
| Input check (GRCh38 via chr13/chr17 lengths, BED inside the loci) | `bin/brca12_prepare.py` | `qc/reference_check.json` |
| Alignment | BWA-MEM, samtools fixmate/sort/markdup/index (single-end and amplicon aware) | `alignment/` |
| Primer clipping (amplicon) | samtools ampliconclip | `alignment/*.primerclip.bam` |
| Read and target QC | samtools stats/idxstats/depth | `qc/read_qc.json` |
| SNVs/indels | GATK HaplotypeCaller (default) or DeepVariant (`--snv_caller deepvariant`) | `calls/<sample>.small.vcf.gz` |
| Mosaic candidates (`--mosaic`) | GATK Mutect2 + FilterMutectCalls | `calls/<sample>.mosaic.vcf.gz` |
| Structural variants | Manta, breakpoints within ±1 Mb of the genes | `calls/<sample>.manta.diploidSV.vcf.gz` |
| CNVs | CNVkit (WGS, hybrid or amplicon; flat or `--cnv_reference`), plus per-exon dosage | `calls/cnvkit/` |
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

## Figures in the report

The HTML report shows two figures per gene: mean depth per exon with the 10x and minimum-depth lines
(uncertain exons in orange) and the exon dosage relative to the chromosome median with the no-hint
band 0.65 to 1.35, an MLPA-like view of possible exon losses or gains. Plain inline SVG, readable in
print; dosage outliers are candidates, not CNV calls.
## Report languages and PDF

`--report_lang en|de|ru` sets the language of the HTML and Markdown report. **The Russian version is
a machine-assisted translation that has not been reviewed by a native speaker**; the texts of the
individual ACMG/AMP hints stay in English there. The JSON keeps German values (stable schema) and
carries the disclaimer and limitations in all three languages.

`--report_pdf` also writes `brca12_report.pdf`, printed from the HTML report with headless
Chrome/Chromium (`--chrome /path/to/chrome`, otherwise searched on `PATH` and in the usual macOS
locations). Without a browser the step is skipped with a note in the log. The PDF has A4 pages,
a header with sample name and RUO note and page numbers.

## Flagging uncertain results

- Exon **uncertain** if mean depth < 20x, fewer than 95 % of bases ≥ 20x, or more than 10 % of reads
  with MAPQ < 20 (pseudogene BRCA1P1, repeats). Thresholds: `--min_exon_depth`, `--min_exon_pct20`,
  `--max_mapq0`.
- Every interval below `--low_depth` (10x) inside the genes is listed.
- Variant quality flags: FILTER, DP < 15, GQ < 20, allele fraction outside 25-75 % for heterozygous
  calls, MQ < 40, location in an uncertain exon.

## ACMG/AMP hints

Automated pre-checks of single criteria following the ClinGen ENIGMA BRCA1/BRCA2 VCEP specifications
v1.2 (Criteria Specification Registry GN092 and GN097). **No classification is assigned.** ClinVar
significance and review status (stars) are shown separately.

| Criterion | Rule as implemented | Needs |
|---|---|---|
| BA1 / BS1 / BS1_Supporting | filtering allele frequency above 0.1 % / above 0.01 % / above 0.002 % in a non-founder population of gnomAD non-cancer, exomes and genomes considered separately | `--gnomad_exomes_vcf`, `--gnomad_genomes_vcf` |
| PM2_Supporting | absent from gnomAD non-cancer exomes and genomes; not for insertions, deletions or delins | both gnomAD files |
| PP3 / BP4 (protein) | missense or in-frame inside a functional domain: BayesDel no-AF ≥ 0.28 (BRCA1) or ≥ 0.30 (BRCA2) for PP3; ≤ 0.15 or ≤ 0.18 together with SpliceAI ≤ 0.1 for BP4 | `--bayesdel_tsv`, `--spliceai_vcf` |
| PP3 (splicing) | SpliceAI ≥ 0.2 for silent, missense, in-frame and intronic variants outside the ±1,2 sites | `--spliceai_vcf` |
| BP1_Strong | silent, missense or in-frame outside the functional domains and SpliceAI ≤ 0.1 | `--spliceai_vcf` |
| BP4 / BP7 | silent inside a domain, or intronic outside ±1,2, with SpliceAI ≤ 0.1; BP7 for silent in a domain and for intronic at or beyond +7/−21 | `--spliceai_vcf` |
| PVS1, PM5_PTC | ENIGMA decision trees, Table 4: strength per exon for termination codons (including the NMD boundaries), per splice-site substitution (±1,2) and for the initiation codon; a candidate hint only where no row applies (e.g. indels at a splice site) | bundled `assets/enigma/pvs1.tsv` |
| PS3 / BS3 | calibrated functional assays, Table 9 (4,729 variants) | bundled `assets/enigma/functional.tsv` |
| PP4 / BP5 | published combined likelihood ratios (multifactorial analysis) | bundled `assets/enigma/multifactorial.tsv` |
| exon DEL / DUP | Table 4 code next to CNV candidates and exon dosage hints (tandem status unknown) | bundled `assets/enigma/pvs1.tsv` |

Functional domains: BRCA1 RING 2-101, coiled-coil 1391-1424, BRCT 1650-1857; BRCA2 PALB2 binding
10-40, DNA binding 2481-3186 (amino acids).

The four evidence files are optional, bgzip-compressed and tabix-indexed, and you download them
yourself: gnomAD is open; SpliceAI scores and BayesDel (for example from dbNSFP) are free for
non-commercial research only and are not redistributed here. `--bayesdel_tsv` expects the columns
chrom, pos, ref, alt, BayesDel no-AF score. Without the gnomAD files the frequency hints fall back
to the highest VEP population frequency and are labelled as approximations. Without SpliceAI the
predictor criteria are reported as open.

The bundled tables come from the [ENIGMA track hub](https://hgdownload.soe.ucsc.edu/hubs/enigma/enigma.html)
(specifications v1.2, ClinGen data, research use only) and are rebuilt with
`bin/brca12_make_enigma.py`; `assets/enigma/SOURCE.json` holds URLs, row counts and checksums. Check
the version against the ClinGen CSpec registry before use.

Not implemented: PS1, PM3, BS2, BS4 and RNA evidence of the sample itself. The gnomAD read-depth conditions (at least 20x for BA1/BS1, 25x
for PM2_Supporting) are not checked; each hint says so. Known exceptions such as BRCA2 p.Lys3326Ter
are noted.

## Legacy names and GRCh37 positions

For comparison with older laboratory reports every small variant also carries its GRCh37 position
and BIC-style legacy names (`grch37Pos`, `bicNames`; shown next to HGVS in the report). Both loci
map to GRCh37 as one contiguous block (Ensembl assembly mapping): BRCA1 −1,847,983 bp, BRCA2
+574,137 bp. BIC numbering: BRCA1 c. + 119, BRCA2 c. + 228 (5'UTR one more, as c. has no position
0); introns use the legacy exon numbers (BRCA1 has no legacy exon 4). BIC itself is inconsistent for
deletions and insertions, so both the HGVS-based and the 5'-shifted name are given
(c.68_69del → 187delAG / 185delAG, c.5266dup → 5385insC / 5382insC, c.5946del → 6174delT).
Against the 4,464 BRCA Exchange entries with a nucleotide BIC name (release of 2026-03-08) one of the
computed names was identical in 89.3 % and matched position and type in 96.9 %; the rest differ in
spelling (for example `343del4` instead of `343delAAAG`). HGVS remains the reference nomenclature.

## Reproducibility

- Pinned versions in `envs/brca12.yml` and `envs/manta.yml` (Manta 1.6 needs Python 2.7; on Apple
  Silicon use `CONDA_SUBDIR=osx-64`).
- The report lists the tool versions that actually ran.
- `brca12_commands.sh` contains every executed command verbatim; `pipeline_info/` holds the
  Nextflow trace, timeline and report.

## Limitations

- Single sample: CNVs are candidates only (no cohort for GATK gCNV); confirm by MLPA.
- Mosaics: at about 30x WGS, spiked variants at 5 % allele fraction were reported in 28 of 94 cases and
  at 20 % in 94 of 94; deeper data (panels, exomes) are needed for low fractions. Complex
  rearrangements are not reliably detected by short reads.
- Exome and panel data cover only what the capture covers; untargeted regions (for example the BRCA2
  3' UTR in both tested exomes) are listed as uncertain, not as negative.
- Mobile element insertions have no dedicated caller. Manta missed real Alu insertions in several
  1000 Genomes carriers; the read-end signature recognises Alu sequence only and gives hints, not
  calls (see `validation/pathogenic/README.md`).
- DeepVariant on macOS is only available via Docker; the default is GATK HaplotypeCaller.

## Tests

```bash
python3 -m unittest discover -s tests   # report logic on synthetic data
bash tests/run_stub_pipeline.sh         # full Nextflow run with stand-in tools
bash tests/run_stub_batch.sh            # two samples from a sample sheet, then the cohort summary
```

Both run in GitHub Actions on every push. They check the pipeline wiring and the reporting logic,
not the accuracy of the external callers.

## Validation (GIAB HG001, HG002)

Small-variant calls were compared with the GIAB v4.2.1 truth sets inside the GIAB high-confidence
part of BRCA1/BRCA2. Reads were taken from GIAB Illumina data at about 30-36x and run from FASTQ.

| Sample | SNV recall / precision | Indel recall / precision | Exons ±20 bp |
|---|---|---|---|
| HG001 | 269/269 · 0.996 | 48/49 · 0.980 | 22/22 found, 0 FP |
| HG002 | 95/95 · 1.000 | 19/20 · 0.950 | 12/12 found, 0 FP |

Both indel misses are genotype errors in intronic homopolymers at 12-14x, and the report flags both.
This is a targeted benchmark: reads were extracted around the genes and realigned to chr13/chr17
only. Annotation and CNVs were not part of it, and neither sample carries a pathogenic variant.
Method, the three discordant calls, limitations and a script to reproduce it are in
[validation/README.md](validation/README.md).

## Validation with pathogenic variants, large deletions and Alu insertions

Three further experiments, all for research use and none a clinical validation:

| Experiment | Result |
|---|---|
| 10 real carriers of pathogenic SNVs/indels (1000 Genomes) | 10 of 10 called heterozygous, GQ 99 |
| 185 small variants spiked into real HG002 reads (SNVs, indels up to 40 bp) | 185 of 185 called, no false positive in the exons |
| 12 large events spiked into real reads (120 bp to whole gene, 2 duplications, 1 Alu) | Manta 12 of 12; CNVkit 3 of 12 (0 of 12 before a fix in 1.1.0); exon dosage hint 9 of 11 |
| Real deletions in BRCA1 (1000 Genomes, not disease-causing) | 1,081 bp: 6 of 6 carriers; 271 bp: 2 of 3 |
| Real Alu insertions in BRCA2 (1000 Genomes, not disease-causing) | called in 4 of 11 carriers; read-end signature in 11 of 11 and in 0 of 47 non-carrier evaluations |

Spiked variants are simulated variants in real reads, not patient samples. Single-sample CNVkit
does not find events below about 10 kb; exon-level events currently rest on Manta. No real sample
with a pathogenic exon deletion, duplication or Alu insertion was available. Details, tables and
limits: [validation/pathogenic/README.md](validation/pathogenic/README.md).

## Validation on exome and amplicon data

Three public GIAB data sets, re-aligned from FASTQ to GRCh38 and compared with GIAB v4.2.1
(`validation/targeted/run_targeted.sh`). Counts are inside the targets and the GIAB confident regions.

| Data set | SNVs found | SNV false positives | Indels | Notes |
|---|---|---|---|---|
| HG002, Agilent SureSelect v5 exome (paired) | 10 of 12 | 0 | 0 of 1 | all misses in the BRCA2 3' UTR with 0-2 reads in the original data; flagged as uncertain |
| HG001, Illumina TruSeq exome (paired, vendor BED) | 15 of 17 | 0 | 1 of 2, 1 FP | same 3' UTR region; the indel FP is a genotype error at 2 reads |
| HG002, Ion AmpliSeq exome (single-end amplicons) | 8 of 8 | 3, all flagged for allele balance | 98 FP | Ion homopolymer errors: indels are flagged with `--platform iontorrent` |

With `--mosaic`, none of the three data sets produced a mosaic candidate. The minimum of 3
supporting reads and 20x depth and the SNV-only rule for Ion Torrent were set after a first run on
these data (which gave 4 candidates for the TruSeq exome, all with 2 supporting reads at 8x, and 23
indel candidates for Ion Torrent), so this is not an independent test of those thresholds. The same
GIAB HG002 WGS benchmark run inside the Docker image gave the same result as the conda run (SNVs
95 of 95, indels 19 of 20, exons ±20 bp 12 of 12, one indel false positive).

## Validation of mosaic candidates

94 small variants (SNVs and indels up to 40 bp) were spiked into real HG002 reads (about 30x) at
allele fractions of 5 to 30 %, plus a negative control (`validation/pathogenic/run_mosaic.sh`).

| Allele fraction | Mosaic candidates | Germline caller | Reported in either list | False positives |
|---|---|---|---|---|
| 5 % | 16 | 12 | 28 of 94 | 0 |
| 10 % | 19 | 46 | 65 of 94 | 0 |
| 15 % | 15 | 70 | 85 of 94 | 0 |
| 20 % | 4 | 90 | 94 of 94 | 0 |
| 30 % | 2 | 92 | 94 of 94 | 0 |
| negative control | 0 | 0 | none spiked | 0 |

Germline calls at these fractions carry the flag "unusual allele balance". At 30x, about 1.5 reads
carry a 5 % variant, which sets the limit; spiked variants are simulated, not patient samples.

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

# Output

```
brca12_results/
  brca12_report.html        standalone report (open in a browser)
  brca12_report.md          same content as Markdown
  brca12_results.json       machine-readable results, schema genome.brca12.v1
  brca12_commands.sh        every executed command, verbatim
  qc/reference_check.json   GRCh38 check, regions and exon counts
  alignment/                BAM + index + duplicate stats (FASTQ input only)
  calls/                    raw HaplotypeCaller/DeepVariant, Manta and CNVkit output
  variants/                 per-gene VEP-annotated VCFs (BRCA1, BRCA2)
  coverage/                 mosdepth tables and per-base depth inside the genes
  pipeline_info/            Nextflow trace, timeline and report
```

## brca12_results.json

Values are stable German keywords (e.g. `"status": "unsicher"`, `"type": "Deletion"`) so that
downstream parsers do not depend on the report language. Main fields:

| Field | Content |
|---|---|
| `schema` | `genome.brca12.v1` |
| `researchUseOnly`, `clinicalClassification` | always `true` / `null` |
| `disclaimer.de`, `disclaimer.en` | the mandatory confirmation notice |
| `genes.BRCA1`, `genes.BRCA2` | locus, strand, MANE transcript, coverage summary, `smallVariants[]` |
| `smallVariants[]` | `pos`, `ref`, `alt`, `genotype`, `depth`, `alleleFraction`, `consequence`, `transcript`, `hgvsc`, `hgvsp`, `grch37Pos`, `bicNames[]`, `exon`, `maxAF`, `gnomadFAF{dataset, population, faf, nonCancer}`, `spliceAI`, `bayesDel`, `clinvar{significance, reviewStatus, stars, condition}`, `acmgHints[{code, textDE, textEN, direction}]`, `qualityFlags[]` |
| `cnvs[]` | CNVkit segments with copy number != 2 overlapping a gene, `affectedExons[]` |
| `exonDosageHints[]` | exons whose depth ratio to the chromosome median is < 0.65 or > 1.35 |
| `insertionSignals[]` | insertion signatures from soft-clipped read ends inside the genes (two clip positions facing each other, overlap 0 to 30 bp, clipped sequence with Alu k-mers; `type` Insertionssignatur) and extra tandem repeat copies (clips at both borders of a repeat unit pointing to the same length, `duplicatedLength`, `type` Tandemkopie-Signatur; for the known 126 bp copy upstream of BRCA1 exon 1 also `markerPosition`, `markerReads`, `evidence` (2 of: clipped reads, marker base, haplotype), `haplotypeTags`, `haplotypeCopies`); hints, not calls |
| `svs[]` | Manta calls with a breakpoint in a gene (`type`, `pos`, `end`, mate, `PR`/`SR`, flags) |
| `exonCoverage[]` | per exon: mean depth, % ≥ 10/20/30x, MAPQ < 20 fraction, dosage ratio, `status`, `reasons[]` |
| `lowCoverageIntervals[]` | intervals below `--low_depth` inside the genes |
| `assay`, `platform` | `wgs`/`exome`/`panel`/`amplicon`; `illumina`/`iontorrent` |
| `input.singleEnd` | `true` for single-end FASTQ (Manta skipped) |
| `readQC` | `regionStats` (reads, mapped/paired/duplicate/MAPQ-0 %, insert size, error rate, read length, base quality), for targeted assays `targets[]` (mean/min depth, % ≥ 20x, `status`) and `targetSummary` (mean depth, % ≥ 20x, `uniformityPct`, `fold80BasePenalty`, `onTargetReadPct`, `weakTargetCount`), for amplicons `primerClipping`, and `flags[]` |
| `mosaicCalling`, `mosaicCandidates[]` | settings of `--mosaic` and Mutect2 calls between the allele-fraction limits that the germline caller did not report; same fields as `smallVariants[]`, `genotype` = `Mosaik-Kandidat` |
| `toolVersions` | versions that actually ran |
| `summary` | counts and `needsAttention` |
| `limitations`, `limitationsEN` | known limits of the method |

## Excel and cohort files

- `brca12_report.xlsx`: sheets Übersicht, SNVs Indels, CNV, SV, Exon-Coverage, Geringe Tiefe, Mosaik-Kandidaten (with `--mosaic`), Lese-QC, Ziel-Coverage (targeted assays), Tool-Versionen.
- Batch (`bin/brca12_batch.py`): `cohort/cohort_samples.tsv`, `cohort/cohort_variants.tsv` (with `cohort_carriers`, `cohort_fraction`, `contamination_hint`), `cohort/cohort_contamination.tsv` (`sample`, variant, `allele_fraction`, `possible_source`), `cohort/cohort.xlsx`.
## PDF and languages

- `brca12_report.pdf` (only with `--report_pdf` and an available Chrome/Chromium).
- JSON `disclaimer{de, en, ru}`, `limitations` (German), `limitationsEN`, `limitationsRU`.

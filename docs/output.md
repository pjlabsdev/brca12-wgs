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
| `insertionSignals[]` | insertion signatures from soft-clipped read ends inside the genes (two clip positions facing each other, overlap 0 to 30 bp, clipped sequence with Alu k-mers); hints, not calls |
| `svs[]` | Manta calls with a breakpoint in a gene (`type`, `pos`, `end`, mate, `PR`/`SR`, flags) |
| `exonCoverage[]` | per exon: mean depth, % ≥ 10/20/30x, MAPQ < 20 fraction, dosage ratio, `status`, `reasons[]` |
| `lowCoverageIntervals[]` | intervals below `--low_depth` inside the genes |
| `toolVersions` | versions that actually ran |
| `summary` | counts and `needsAttention` |
| `limitations`, `limitationsEN` | known limits of the method |

## Excel and cohort files

- `brca12_report.xlsx`: sheets Übersicht, SNVs Indels, CNV, SV, Exon-Coverage, Geringe Tiefe, Tool-Versionen.
- Batch (`bin/brca12_batch.py`): `cohort/cohort_samples.tsv`, `cohort/cohort_variants.tsv` (with `cohort_carriers`, `cohort_fraction`, `contamination_hint`), `cohort/cohort_contamination.tsv` (`sample`, variant, `allele_fraction`, `possible_source`), `cohort/cohort.xlsx`.

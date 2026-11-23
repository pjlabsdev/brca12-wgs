# Output

```
brca12_results/
  SAMPLE.tsv                one line per gene
  SAMPLE.findings.tsv       one line per finding with its certainty
  SAMPLE.exons.tsv          one line per exon with its certainty
  SAMPLE.json               machine-readable results, schema genome.brca12.v1
  SAMPLE.md                 the report as Markdown, English and German
  SAMPLE.xlsx               the tables as an Excel workbook
  brca12_report.html        standalone report (open in a browser): German and English, tabs Result/Details
  brca12_commands.sh        every executed command, verbatim
  qc/reference_check.json   GRCh38 check, regions and exon counts
  alignment/                BAM + index + duplicate stats (FASTQ input only)
  calls/                    raw HaplotypeCaller/DeepVariant, Delly or Manta and CNVkit output
  variants/                 per-gene VEP-annotated VCFs (BRCA1, BRCA2)
  coverage/                 mosdepth tables and per-base depth inside the genes
  pipeline_info/            Nextflow trace, timeline and report
```

`SAMPLE` is the value of `--sample` (the start file takes the input file name). The three tables
follow the form of hla-wgs (`PREFIX.tsv`, `PREFIX.conf.tsv`): tab-separated, one header line, English
values, `-` for an empty cell. Every value comes from `SAMPLE.json`. Up to 1.1 the JSON was called
`brca12_results.json`, the Markdown `brca12_report.md` and the workbook `brca12_report.xlsx`.

## SAMPLE.tsv

One line per gene (BRCA1, BRCA2).

| column | meaning |
|---|---|
| `gene` | gene |
| `result` | `attention`, `unclear` or `none`, as the overall statement of the result page but for this gene |
| `attention`, `unclear` | number of entries of each level in `SAMPLE.findings.tsv` |
| `small_variants` | small variants reported for the gene (all, also benign ones) |
| `copy_number` | CNVkit CNVs and exon copy number calls |
| `structural` | SV calls with a breakpoint in the gene |
| `insertion_signals` | insertion and tandem copy signals, including known common alleles |
| `mosaic` | mosaic candidates (with `--mosaic`) |
| `coverage` | `limited` when an exon is uncertain or the mean exon depth is below `--min_exon_depth`, else `ok`; with `limited` a `none` does not exclude a finding in the weak parts |
| `exons`, `exons_uncertain` | exons with a depth value, and those marked uncertain in `SAMPLE.exons.tsv` |
| `mean_depth`, `min_exon_depth` | mean of the exon mean depths, and the lowest exon mean depth (MAPQ ≥ 20) |
| `low_depth_bp` | bases inside the gene below `--low_depth` |
| `transcript` | MANE Select transcript |

## SAMPLE.findings.tsv

One line per finding, in the order of the report: small variants, mosaic candidates, copy number
(CNVkit, exon copy number calls and hints), structural variants, insertion signals, exon dosage hints.

| column | meaning |
|---|---|
| `gene` | gene |
| `kind` | `small`, `mosaic`, `cnv`, `exon_copy_number`, `sv`, `insertion_signal`, `exon_dosage` |
| `level` | `attention` (may point to a predisposition: ClinVar pathogenic or likely pathogenic, loss of function, CNV, exon copy number call, SV with FILTER PASS), `unclear` (ClinVar VUS, mosaic candidate, SV without PASS, insertion signal), `common` (a known common allele, reported but not counted), `hint` (a lead that is not counted: exon copy number hint, dosage ratio) or `none` (listed without a finding, for example benign); `attention` and `unclear` are the entries the result page counts |
| `reason` | `clinvar_pathogenic`, `loss_of_function`, `clinvar_vus`, `mosaic`, `cnv`, `exon_copy_number`, `sv`, `insertion_signal`, `common_allele`, `exon_copy_number_hint`, `dosage_ratio`, `ba1` (a small variant with a BA1 code) or `listed` |
| `chrom`, `start`, `end` | GRCh38, 1-based; for small variants `end` is the last reference base |
| `change` | `REF>ALT`, or the type with copy number, caller or label |
| `hgvsc`, `hgvsp` | on the MANE transcript (small variants) |
| `exon` | exon or intron of a small variant (`10/23`), affected exons of a larger event |
| `genotype` | `heterozygous`, `homozygous`; VCF genotype for SVs; zygosity for known alleles |
| `reads` | small variants `alt/depth`; SVs `PR ref,alt SR ref,alt`; CNVkit `log2`; exon copy number `quality`; insertion signals the read counts that support them |
| `allele_fraction` | alternative allele fraction (small variants) |
| `clinvar`, `clinvar_stars` | ClinVar significance and review stars |
| `acmg` | the ACMG/AMP codes given as hints (ENIGMA specification), for example `PVS1 (candidate)`, `BA1 (approximation)`, `PP3/BP4 open` |
| `certainty` | `ok`: the call carries no quality flag; `limited`: it carries one (low depth, strand bias, known weak spot, CNV flags), it is a mosaic candidate or an SV without FILTER PASS; `hint`: a read signal or depth ratio, not a call |
| `flags` | the quality flags and notes, English, separated by `; ` |

Every entry needs confirmation with a clinically validated method before any medical use.

## SAMPLE.exons.tsv

One line per coding exon of the MANE transcripts.

| column | meaning |
|---|---|
| `gene`, `exon`, `label` | gene, exon number, and the name with the legacy BRCA1 numbering |
| `chrom`, `start`, `end` | GRCh38, 1-based |
| `mean_depth` | mean depth (MAPQ ≥ 20) |
| `pct10x`, `pct20x`, `pct30x` | percent of the bases with at least 10, 20 and 30 reads |
| `low_mapq_fraction` | fraction of reads with MAPQ below 20 (pseudogene or repeat) |
| `dosage_ratio` | exon depth relative to the chromosome median, scaled to the exons; `-` when the baseline is not plausible |
| `status` | `ok` or `uncertain` (depth below `--min_exon_depth`, less than `--min_exon_pct20` % at 20x, or too many reads with low MAPQ) |
| `reasons` | why the exon is uncertain |

A variant in an uncertain exon can be missed; an empty finding list says nothing about those exons.

## SAMPLE.json

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
| `cnvSegmentsNotSupported[]` | CNVkit segments with copy number != 2 that reach beyond the gene while the bins inside the gene lie in CNVkit's two-copy range (mean log2 above -0.25 and at most 0.2): `segmentType`, `overlapBinMeanLog2`, `supportedInGene` false; not candidates, each adds a limitation |
| `exonDosageHints[]` | exons whose depth ratio to the chromosome median is < 0.65 or > 1.35 |
| `exonDosageBaseline` | `plausible` (false when the chromosome median lies below a quarter of the median exon depth, e.g. when only the gene regions were read; dosage ratios and hints are then switched off and a limitation is added), `chromosome`, `chromosomeMedian`, `medianExonDepth` |
| `exonCopyNumber` | exon copy number against a reference cohort: `enabled` (false without a panel for the assay, with `reason`), `panel`, `panelSize`, `pad`, `minQuality`, `normaliserDepth`, `exons[]` (`log2Deviation`, `z`, `copyNumber`, `pTwoCopies`) and `calls[]` (`type` Deletion or Duplikation, `copyNumber`, `exons`, `affectedExons`, `quality`, `meanLog2Deviation`, `flags`, `enigmaCodes`), `hints[]` (same fields, quality 5 to 9), `hintQuality`, `alleleBalance`, `hetSitesUsed`; calls and hints also appear in the CNV table of the report |
| `insertionSignals[]` | insertion signatures from soft-clipped read ends inside the genes (two clip positions facing each other, overlap 0 to 30 bp, clipped sequence with Alu k-mers; `type` Insertionssignatur) and extra tandem repeat copies (clips at both borders of a repeat unit pointing to the same length, `duplicatedLength`, `type` Tandemkopie-Signatur; for the known 126 bp copy upstream of BRCA1 exon 1 also `markerPosition`, `markerReads`, `evidence` (2 of: clipped reads, marker base, haplotype), `haplotypeTags`, `haplotypeCopies`, `realignedReads` from allele realignment; the known BRCA2 Alu allele also has `zygosity`; known alleles from the catalog have `label`, `realignedReads`, `zygosity` (`unknown` for a tandem duplication), `clinicalSignificance` (`pathogenic` for the founder rearrangements, `type` Bekannte pathogene Insertion, Deletion or Duplikation) and for duplications `discordantPairs`); tandem duplications from outward-facing read pairs (`signature` tandem duplication, `type` Tandemduplikation (Lesepaare), `outwardPairs`, `splitReads`, `duplicatedLength`, `breakpointPrecision` approximate); hints, not calls |
| `svs[]` | SV calls with a breakpoint in a gene (`caller` Delly or Manta, `type`, `pos`, `end`, mate, `pairedReads` and `splitReads` as reference,alternative: Manta PR/SR, Delly DR,DV and RR,RV; `length`, flags) |
| `exonCoverage[]` | per exon: mean depth, % ≥ 10/20/30x, MAPQ < 20 fraction, dosage ratio, `status`, `reasons[]` |
| `lowCoverageIntervals[]` | intervals below `--low_depth` inside the genes |
| `assay`, `platform` | `wgs`/`exome`/`panel`/`amplicon`; `illumina`/`iontorrent` |
| `input.singleEnd` | `true` for single-end FASTQ (SV caller skipped) |
| `readQC` | `regionStats` (reads, mapped/paired/duplicate/MAPQ-0 %, insert size, error rate, read length, base quality), for targeted assays `targets[]` (mean/min depth, % ≥ 20x, `status`) and `targetSummary` (mean depth, % ≥ 20x, `uniformityPct`, `fold80BasePenalty`, `onTargetReadPct`, `weakTargetCount`), for amplicons `primerClipping`, and `flags[]` |
| `mosaicCalling`, `mosaicCandidates[]` | settings of `--mosaic` and Mutect2 calls between the allele-fraction limits that the germline caller did not report; same fields as `smallVariants[]`, `genotype` = `Mosaik-Kandidat` |
| `toolVersions` | versions that actually ran |
| `summary` | counts and `needsAttention` |
| `limitations`, `limitationsEN` | known limits of the method |

## Excel and cohort files

- `SAMPLE.xlsx`: sheets Übersicht, SNVs Indels, CNV, SV, Exon-Coverage, Geringe Tiefe, Mosaik-Kandidaten (with `--mosaic`), Lese-QC, Ziel-Coverage (targeted assays), Tool-Versionen.
- Batch (`bin/brca12_batch.py`): `cohort/cohort_samples.tsv`, `cohort/cohort_variants.tsv` (with `cohort_carriers`, `cohort_fraction`, `contamination_hint`), `cohort/cohort_contamination.tsv` (`sample`, variant, `allele_fraction`, `possible_source`), `cohort/cohort.xlsx`.
## Languages

- JSON `disclaimer{de, en}`, `limitations` (German), `limitationsEN`.

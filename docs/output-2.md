# Output (brca12-wgs 2.0, in development)

`brca12-wgs INPUT.bam|INPUT.cram OUTPREFIX [--reference REF.fa]` writes three files, in the form of
the hla-wgs output. The HTML result page is made from `OUTPREFIX.json` at the end. Research use only.

A column that this version does not determine yet is `-` in the TSV files and `null` in the JSON;
`analysed` and `not_yet_analysed` in the JSON list the categories. Development state: coverage and
read QC are done (equal to 1.1); small variants, annotation, copy number and structural variants
follow.

## OUTPREFIX.tsv

One line per gene.

| column | meaning |
|---|---|
| `gene` | BRCA1 or BRCA2 |
| `transcript` | MANE Select transcript (BRCA1 NM_007294.4, BRCA2 NM_000059.4) |
| `result` | `pathogenic` (ClinVar pathogenic or likely pathogenic), `candidate` (finding that needs confirmation) or `none` |
| `small_variants` | SNVs and indels found in the gene region |
| `clinvar_plp` | of these, ClinVar pathogenic or likely pathogenic (not conflicting) |
| `loss_of_function` | loss-of-function candidates (ENIGMA PVS1 hint) |
| `exon_copy_number` | exon deletions or duplications called against the reference cohort |
| `cnv_candidates` | CNV segments supported by the gene's own bins |
| `structural_variants` | structural variants with a breakpoint in the gene |
| `coding_depth` | mean depth over the exons (MAPQ >= 20), weighted by exon length |
| `pct20x` | share of exon bases with depth >= 20, weighted by exon length |
| `uncertain_exons` | exons marked uncertain (see `OUTPREFIX.exons.tsv`), comma-separated |
| `low_coverage_bp` | bases of the gene region (introns included) with depth < 10 |
| `statement` | the overall statement for the gene in one sentence: result (pathogenic variant found / finding that needs confirmation / no pathogenic variant found), exon coverage, and which exons are uncertain or that all exons were checked; in the JSON as `statement` with `en` and `de` |

## OUTPREFIX.exons.tsv

One line per exon of the MANE Select transcripts (BRCA1 23, BRCA2 27).

| column | meaning |
|---|---|
| `gene`, `exon`, `label` | exon number; BRCA1 exons 4 to 23 carry the legacy number in brackets |
| `chrom`, `start`, `end`, `length` | 1-based, inclusive (GRCh38) |
| `mean_depth` | mean depth with MAPQ >= 20, records with flags 0x704 skipped, the overlap of two mates of a properly paired fragment counted once (as mosdepth) |
| `pct10x`, `pct20x`, `pct30x` | share of exon bases with at least that depth |
| `low_mapq_fraction` | share of depth lost by the MAPQ >= 20 filter (pseudogene or repeat) |
| `dosage_ratio` | mean depth relative to the median depth of the chromosome (100 kb windows) |
| `status` | `ok` or `uncertain` |
| `reasons` | why uncertain: mean depth < 20x, less than 95 % of bases >= 20x, more than 10 % of reads with MAPQ < 20 |

## OUTPREFIX.json

Schema `genome.brca12.v2`.

```json
{
  "tool": "brca12-wgs", "version": "2.0.0-dev", "schema": "genome.brca12.v2",
  "sample": "PJ", "input": {"files": ["PJ.bam"], "reference": null},
  "analysis_time": 1795426512, "research_use_only": true,
  "analysed": ["coverage", "read_qc"],
  "not_yet_analysed": ["small_variants", "annotation", "exon_copy_number", "cnv", "structural_variants", "insertion_signals"],
  "thresholds": {"min_mapq": 20, "min_exon_depth": 20.0, "min_exon_pct20x": 95.0, "max_low_mapq_fraction": 0.1,
                 "low_coverage_depth": 10.0, "calling_padding": 200, "chromosome_window": 100000},
  "read_qc": {"reads_in_regions": 48298, "mapped_pct": 99.98, "properly_paired_pct": 98.42, "duplicate_pct": 0.0,
              "mapq0_pct": 0.37, "median_insert_size": 290, "mean_insert_size": 291.9, "error_rate": 0.00304,
              "average_read_length": 150.0, "average_base_quality": 35.7},
  "chromosome_median": {"chr13": 33.335, "chr17": 35.05},
  "genes": [{"gene": "BRCA1", "transcript": {"refseq": "NM_007294.4", "ensembl": "ENST00000357654"}, "result": null,
             "statement": {"en": "BRCA1: variants not analysed yet (coverage only); exons 34.1x, 97.0 % of exon bases >= 20x; uncertain: Exon 6 (Legacy 7), Exon 23 (Legacy 24)", "de": "..."},
             "small_variants": null, "clinvar_plp": null, "loss_of_function": null, "exon_copy_number": null,
             "cnv_candidates": null, "structural_variants": null,
             "coverage": {"exon_count": 23, "coding_mean_depth": 34.1, "pct20x": 97.0,
                          "uncertain_exons": ["Exon 6 (Legacy 7)", "Exon 23 (Legacy 24)"], "low_coverage_bp": 2473,
                          "low_coverage_intervals": [{"chrom": "chr17", "start": 43044805, "end": 43044805, "length": 1, "min_depth": 1.0}, ...]}}],
  "exons": [{"gene": "BRCA1", "exon": 6, "label": "Exon 6 (Legacy 7)", "chrom": "chr17", "start": 43104122, "end": 43104261,
             "length": 140, "mean_depth": 22.7, "pct10x": 100.0, "pct20x": 82.9, "pct30x": 0.0, "low_mapq_fraction": 0.002,
             "dosage_ratio": 0.65, "status": "uncertain",
             "reasons": {"en": ["only 82.9 % of bases >= 20x"], "de": ["nur 82.9 % der Basen >= 20x"]}}]
}
```

`read_qc` follows samtools stats over the gene regions +/- 200 bp: records overlapping the region,
per-read counts for primary records, insert sizes from the main 99 % of pairs. Records without an
NM tag (CRAM) get NM computed from the reference as htslib does.

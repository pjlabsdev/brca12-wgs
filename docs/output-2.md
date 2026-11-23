# Output (brca12-wgs 2.0, in development)

`brca12-wgs INPUT.bam|INPUT.cram OUTPREFIX --reference REF.fa` writes four files, in the form of
the hla-wgs output. The HTML result page is made from `OUTPREFIX.json` at the end. Research use only.

Exon copy number is the method of 1.1 (bin/brca12lib/exon_cnv.py), ported one to one: mean depth
over exon +/- 200 bp (as samtools depth -a -Q 20) against eight single-copy control windows, the
panel of 80 genomes, a hidden Markov model over the exons of a gene with allele balance at the own
heterozygous SNV calls. A segment with quality >= 10 is a call (gene result `candidate`, confirmation
e.g. by MLPA), 5-9 a hint in the JSON (`exon_copy_number.hints`).

Without reads in the gene regions (truncated file, other reference or chromosome names) the program
stops with an error instead of writing a result.

A column that this version does not determine yet is `-` in the TSV files and `null` in the JSON;
`analysed` and `not_yet_analysed` in the JSON list the categories. Development state: coverage,
read QC, small variants and their annotation are done; exon copy number, CNV and structural
variants follow.

## OUTPREFIX.tsv

One line per gene.

| column | meaning |
|---|---|
| `gene` | BRCA1 or BRCA2 |
| `transcript` | MANE Select transcript (BRCA1 NM_007294.4, BRCA2 NM_000059.4) |
| `result` | `pathogenic` (a PASS variant that ClinVar lists as pathogenic or likely pathogenic, not conflicting), `candidate` (a loss-of-function variant that is not benign in ClinVar and rarer than 0.1 % in gnomAD, or a ClinVar pathogenic variant that did not pass the filters: needs confirmation), `none`, or `not_assessable` (no finding, but the exons are covered too little to say none: coding depth below 15x or less than 80 % of exon bases >= 20x) |
| `small_variants` | SNVs and indels with FILTER PASS in the gene region |
| `clinvar_plp` | of these, ClinVar pathogenic or likely pathogenic (not conflicting) |
| `loss_of_function` | of these, frameshift, stop gained, start lost, splice donor or acceptor (+/-1, 2), not benign |
| `exon_copy_number` | exon deletions or duplications against the reference cohort (quality >= 10), `none`, or `-` when not determined (a region BAM/CRAM without the panel's control windows elsewhere on chr13/chr17) |
| `cnv_candidates` | CNV segments supported by the gene's own bins |
| `structural_variants` | structural variants with a breakpoint in the gene |
| `coding_depth` | mean depth over the exons (MAPQ >= 20), weighted by exon length |
| `pct20x` | share of exon bases with depth >= 20, weighted by exon length |
| `uncertain_exons` | exons marked uncertain (see `OUTPREFIX.exons.tsv`), comma-separated |
| `low_coverage_bp` | bases of the gene region (introns included) with depth < 10 |
| `statement` | the overall statement for the gene in one sentence: result (pathogenic variant found / finding that needs confirmation / no pathogenic variant found), exon coverage, and which exons are uncertain or that all exons were checked; in the JSON as `statement` with `en` and `de` |

## OUTPREFIX.variants.tsv

One line per called allele in the gene regions +/- 200 bp (a site with two different alleles gives
two lines). All calls are listed, also those that did not pass the filters (column `filter`).

| column | meaning |
|---|---|
| `gene` | BRCA1 or BRCA2 |
| `chrom`, `pos`, `ref`, `alt` | the allele as in a VCF (GRCh38, 1-based, left-aligned, full REF and ALT sequences) |
| `genotype`, `zygosity` | `0/1` heterozygous, `1/1` homozygous, `1/2` two different alleles |
| `filter` | `PASS`, or `QD` (QUAL / alt reads < 2), `FS` (strand bias, Fisher > 60 SNV / > 200 indel) |
| `qual`, `gq` | variant quality (phred probability that there is no variant) and genotype quality |
| `depth`, `ref_reads`, `alt_reads`, `allele_fraction` | reads over the site and the reads assigned to REF and ALT |
| `hgvs_c` | HGVS c. on the MANE Select transcript (NM_007294.4, NM_000059.4), 3'-shifted |
| `hgvs_p` | HGVS p. (three-letter amino acids, `Ter`) on NP_009225.1 / NP_000050.3 |
| `consequence` | Sequence Ontology terms, most severe first (as Ensembl VEP) |
| `exon`, `intron` | number / total on the MANE transcript |
| `clinvar_id`, `clinvar_significance`, `clinvar_review`, `clinvar_stars` | ClinVar (version in the JSON) |
| `gnomad_present`, `gnomad_faf95`, `gnomad_population`, `gnomad_dataset` | gnomAD v2.1.1 exomes and v3.1.2 genomes: highest filtering allele frequency (FAF95) over the populations, founder populations skipped, non-cancer subset preferred (as 1.1) |
| `class` | the tool's hint (no ACMG classification): `pathogenic` (ClinVar P/LP), `benign` (ClinVar B/LB or FAF95 >= 0.1 %, ENIGMA BA1), `candidate` (loss of function, not benign), `other` |
| `in_gene_region` | `yes` inside the gene region, `no` in the 200 bp padding |

The HGVS names and consequences come from the program itself. Checked against Ensembl VEP 116
(`--hgvs --shift_hgvs 1`) on all 37,125 ClinVar alleles in BRCA1/BRCA2: c. equal for 99.97 %, p.
for 99.96 %, most severe consequence for 99.70 %. The differences: a few very large insertions and
deletions across the transcript ends; VEP's `p.LeuLeu1521=` form; and deletions that reach the
coding exon only after the 3' shift (for example c.518del), where VEP gives no p. but this program
does, because the protein change is the same.

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
| `cnv_log2_deviation`, `cnv_copy_number`, `cnv_p_two_copies` | exon copy number against the panel of 80 1000 Genomes 30x genomes: log2 deviation from the panel centre, copy number of the hidden Markov model, posterior probability of two copies |
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

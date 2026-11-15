# GIAB HG003 to HG007 (before release 1.1.0)

Plan frozen in `PREREGISTRATION.md` (commit bc270b2) before the first run. Code: `main` ce825de
(Delly 1.3.3 default). Arm A, targeted method identical to HG001/HG002 (`run.sh`): reads over
BRCA1/BRCA2 ±10 kb from GIAB NHGRI Illumina BAMs, ~30x, realigned from FASTQ. VEP and CNVkit are
stand-ins, so annotation and CNVs are not tested here. Arm B (whole genomes on the Mac with real
VEP) follows when the downloads are done.

## Result, arm A

GIAB confident regions inside BRCA1/BRCA2 (`results/summary.tsv`, errors in `results/errors.tsv`).

| Sample | Truth | SNVs TP/truth, FP | SNV recall · precision | Indels TP/truth, FP | Exons ±20 bp |
|---|---|---|---|---|---|
| HG003 | v4.2.1 | 263/263, 1 FP | 1.000 · 0.996 | 47/47, 0 FP | 22/22, 0 FP |
| HG004 | v4.2.1 | 86/86, 2 FP | 1.000 · **0.977** | 21/21, 1 FP | 12/12, 0 FP |
| HG005 | v4.2.1 | 252/255, 0 FP | **0.988** · 1.000 | 31/32, 0 FP | 23/23, 0 FP |
| HG006 | v3.3.2 | 184/184, 0 FP | 1.000 · 1.000 | 38/38, 0 FP | 24/24, 0 FP |
| HG007 | v3.3.2 | 106/106, 0 FP | 1.000 · 1.000 | 21/22, 0 FP | 13/13, 0 FP |
| All | | 891/894, 3 FP | 0.997 · 0.997 | 158/160, 1 FP | 94/94, 0 FP |

**The preregistered test is not passed in full.** The indel and exon thresholds hold for every
sample, but two samples miss an SNV threshold of 0.99:

- HG004, precision 0.977: the two SNV false positives and the indel false positive are one
  haplotype at chr17:43126851-43126856 (intron), built by HaplotypeCaller from 4 of 17 reads. Two
  of the three records carry the report flag for an unusual allele fraction.
- HG005, recall 0.988: three intronic heterozygous SNVs (chr17:43096467, 43102469, 43125170) were
  not called. Local depth there is 8 to 19 reads with 0, 5 and 3 alternative reads, far below the
  35x gene mean.
- HG003's single false positive (chr17:43084901, 4 of 29 reads) is flagged as well.

Every exonic variant of the five genomes (±20 bp) was found, and there was no false positive in
the exons. No structural variant or copy-number change was called in any of the five genes.
Following the plan, nothing was changed after these results; any fix comes as a separate change and
is reported before and after.

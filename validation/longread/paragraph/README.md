# Head-to-head with Paragraph (preregistered samples)

Research use only. Not part of the preregistration; run afterwards on the same 200 samples and the
same short reads (windows of the 126 bp copy, the BRCA2 Alu insertion and the 980 bp BRCA1 deletion),
with the same allele sequences (`variants.vcf`). Paragraph 2.3 (Chen et al., Genome Biol 2019)
from bioconda, `multigrmpy.py` with depth 30 and read length 150 (`run_paragraph.sh`). Per sample:
`paragraph_vs_ours.tsv`.

| Allele | Paragraph, any call | Paragraph, PASS | This pipeline | Non-carrier calls (both) |
|---|---:|---:|---:|---:|
| 126 bp tandem copy | 114/116 | 102/116 | 116/116 | 4/84 |
| BRCA2 Alu insertion | 23/25 | 12/25 | 25/25 | 0/175 |
| 980 bp BRCA1 deletion | 2/2 | 2/2 | 2/2 | 0/198 |

Paragraph calls the same 4 long-read non-carriers as this pipeline (HG01527, HG01566, NA18940,
NA19789), which supports that they carry the insertion and are missed by the thin long-read
coverage. Zygosity agreement with the long reads among called carriers: tandem copy 113/114
(this pipeline, from the haplotype: 113/116), Alu 23/23 (this pipeline 25/25).

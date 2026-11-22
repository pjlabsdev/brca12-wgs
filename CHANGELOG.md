# Changelog

## Unreleased (1.1.0)

- Read QC: the sequencing platform is recognised from the read-name format (MGI/DNBSEQ or Illumina,
  >= 90 % of up to 1,000 reads of the first calling region). MGI/DNBSEQ reads give a note that the
  pipeline is validated on Illumina data, and a second one when the read group says ILLUMINA. A BAM
  without a duplicate-marking step in its header (samtools/sambamba markdup, Picard/GATK
  MarkDuplicates, bammarkduplicates, umi_tools dedup; a file name such as x.markdup.bam does not
  count) and with 0 % duplicates gives a note too. Seen on a DNBSEQ genome aligned with bwa and
  sorted only, with ILLUMINA in its read group.
- A task killed by a signal (exit > 128, for example a bus error 138 or kill 137 when an external USB
  disk drops out briefly) is run again, at most `--signal_retries` times (2); other errors stop the
  run as before.
- Start file (`start/brca12-wgs-start.command`): removes the macOS download quarantine from its own
  folder so the shipped programs run on other Macs, and after a disk drop-out waits until the input
  and the reference are readable again and continues with `-resume` (at most `BRCA12_RESUME`
  times, 3). The output folder is created only after the start is confirmed, and an unreadable
  confirmation stops the start. Both points follow the hla-wgs start file.
- HTML report rebuilt after the hla-wgs result page: one file with German and English (switch at the
  top, `--report_lang` sets the language shown first), two tabs (Result with report header, counts
  and the entries to review; Details with all tables, figures, versions, commands and limitations),
  PJ Labs design with the embedded font, explanations on hover or tap, DOI and link to the research
  article. No script is used. A notice appears when the mean exon depth of a gene is below
  `--min_exon_depth`, and the header states the quality as passed or limited with the reason. The
  PDF prints both tabs with its own footer (sample, RUO note, page n / N). New JSON field
  `pipeline` (name, version, DOI).
- Paths with spaces: the pipeline scripts and the reference path are quoted in the process
  commands. Before, a run from a folder such as `brca12-wgs Start 2` failed in the first step with
  `can't open file '/…/brca12-wgs'`, and a reference path with a space failed in the report step.
  The stub test now runs once with a space in the pipeline, input and output folders.
- CNVkit candidates: a segment with copy number != 2 that reaches beyond the gene counts as a
  candidate only if the bins inside the gene show the change too. If their mean log2 lies in
  CNVkit's own two-copy range (above -0.25, at most 0.2) the segment is listed under
  `cnvSegmentsNotSupported` with a limitation instead. In a 30x DNBSEQ genome two segments of 96 Mb
  and 24 Mb (log2 -0.30 and -0.59 against the flat reference) had been reported as whole-gene
  deletion candidates of BRCA2 and BRCA1 although the gene bins were at -0.05 and -0.10 and both
  genes carried heterozygous variants. Segments inside a gene and the simulated deletions of the
  spike-in tests are unchanged.
- Exon coverage: the reads overlapping the analysed intervals are fetched once through the index
  (`samtools view -M -L`) into a small BAM, and mosdepth, `samtools depth` and the insertion
  signals read that file. Before, `mosdepth --by` and `samtools view -L` read the whole alignment
  twice or more, which took hours on a 70 GB genome. Depths are identical (checked on GIAB HG003).
- CNVkit: `Rscript` is resolved explicitly (next to `--cnvkit`, else `PATH`, or `--rscript`) and
  passed with `--rscript-path`; before the run it must load DNAcopy, after `batch` the segments and
  the .cnr must exist, otherwise the step stops with a clear message. Before, an `Rscript` without
  DNAcopy earlier on `PATH` (a conda base) made `batch -p` end without segments and without an
  error message (whole-genome NA18949 run). The DNAcopy version is listed in the tool versions.
- Per-exon dosage: when the chromosome median lies below a quarter of the median exon depth (for
  example when only the gene regions were read), the ratios and hints are switched off with a
  limitation instead of one hint per exon (`exonDosageBaseline`).
- Small variants: before HaplotypeCaller and Mutect2 the reads of the calling regions are written to
  a BAM in which soft-clipped poly-G read-through tails (two-colour chemistry) are hard-clipped
  (`bin/brca12_clip_polyg.py`). GATK assembles soft-clipped bases and made a false three-record
  haplotype from such tails in GIAB HG004. The same step fixes CRAM input on which htsjdk fails
  while converting the .crai ("Attempt to create a bai entry for an unmapped slice", for example
  the 1000 Genomes 30x CRAMs).
- GIAB HG003 to HG007, preregistered (`validation/giab_more/`): every exonic variant found (94 of
  94), no false positive in the exons. Before the poly-G fix 891/894 SNVs with 3 false positives
  and 158/160 indels with 1; two samples missed the 0.99 SNV threshold (HG004 precision 0.977,
  HG005 recall 0.988). After it 1 SNV false positive and no indel false positive; HG005 still
  misses three intronic SNVs in 2x250 bp reads.
- Structural variants: Delly 1.3.3 is the default SV caller (`--sv_caller delly|manta`); it runs
  natively on Apple Silicon. PASS calls with a non-reference genotype; DEL/DUP/INV need both ends in
  the search windows (genes +/- `--sv_padding`) and at least 50 bp, as Manta's call regions
  (`bin/brca12_filter_sv.py`). Manta 1.6 remains available. Report reads Delly DR/DV and RR/RV.
- Delly 1.3.3 versus Manta 1.6, preregistered head-to-head (`validation/sv_delly/`): 40 simulated
  deletions (50 bp-20 kb), 20 inversions and 40 tandem duplications at three depths: Delly 108/120,
  60/60 and 118/120, Manta 105/120, 60/60 and 109/120, no false call for either in 309 simulated analyses;
  200 real 1000 Genomes genomes: 1 call not confirmed by ONT long reads; all 9 ONT-confirmed carriers of
  the 980 bp BRCA1 deletion found. `spike.py` simulates inversions (LINV). Read-pair deletion/inversion
  scan added as an opt-in research arm (`--read-pair-svs`).
- Report languages: English and German only; the Russian report language (`--report_lang ru`,
  added earlier in this release cycle) has been removed.
- Tandem duplications from outward-facing read pairs (`insertionSignals[]`, type
  Tandemduplikation (Lesepaare)): at least 4 outward pairs or joining split reads, MAPQ >= 20,
  span 1 to 200 kb, independent of read depth. Preregistered test with 40 simulated exon-spanning
  duplications (1-20 kb) in real HG002 reads at three depths: 119/120 found (Manta 109/120), no false
  call in 147 simulated analyses, one call in 200 real 1000 Genomes genomes that the ONT reads do not
  confirm. `validation/tandem_dup/`.
- Exon copy number against a reference cohort (`exonCopyNumber`, `--exon_cnv_panel`): depth over
  each exon +/- 200 bp relative to eight single-copy control windows, compared with a panel of 80
  1000 Genomes 30x genomes (`assets/exon_cnv_panel_1kg30x.json`, WGS default), hidden Markov model
  over the exons of each gene, segments with quality >= 10. Rules frozen before a test on 60 new
  1000 Genomes samples with simulated events: 569/600 deletions, 209/360 duplications, BRCA1 exon
  1-2 via BRCA1P1 60/60 and 51/60, no call in the 60 unmodified genomes (current dosage rule: 87
  hints). `validation/exon_cnv/`.
- Exon copy number, step 2: segments with quality 5 to 9 are reported as hints (`exonCopyNumber.hints`),
  and allele balance at heterozygous SNVs of the small-variant calls enters the model. Panel rebuilt
  from reads of the full window. Frozen test on 58 further samples, events simulated on real
  haplotypes: 567/580 deletions and 244/348 duplications at quality >= 5, no false call.
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
  Clip positions a few bases apart are pooled, clips from 10 bp count, and the common 126 bp copy
  upstream of BRCA1 exon 1 is also recognised by its marker base (chr17:43128365 C) and by tag SNPs
  of its haplotype; 2 of the 3 signals are required, and the haplotype gives the likely zygosity.
  It is found in 19 of 19 HPRC and 6 of 6 ONT carriers (with allele realignment, below), none in 25 non-carriers, and in an
  independent set of 120 1KG-ONT samples in 67 of 67 carriers and 1 of 53 non-carriers.
- Known-allele genotyping by allele-aware realignment (`brca12_insertion_signals.py`): reads near
  the 126 bp copy and the common intronic BRCA2 Alu insertion are aligned to the reference and the
  allele haplotype (HPRC sequences, Myers bit-vector edit distance, no new dependency); fields
  `realignedReads`, `zygosity`. Preregistered test on 200 new 1KG-ONT samples
  (`validation/longread/frozen/`): tandem copy 116/116 carriers (4/84 non-carriers with calls),
  Alu 25/25 and 0/175 with 25/25 correct genotypes. Deletions in the catalogue: 980 bp BRCA1
  intronic deletion (registered addendum: 10/10 carriers, 0/198 non-carriers, genotype 10/10).
  Head-to-head on the same reads: ExpansionHunter v5 14/19 HPRC and 54/64 independent carriers of
  the tandem copy; Paragraph 114/116 and 23/25 (PASS 102 and 12), same 4 non-carrier calls.
- Pathogenic founder rearrangements in the known-allele catalog, flagged as pathogenic in the
  report (`clinicalSignificance`): BRCA2 c.156_157insAlu (AluYa5 consensus with the published 14 bp
  target site duplication), BRCA1 exon 13 duplication (c.4186-1787_4357+4122dup; realignment across
  the whole segment plus outward-facing read pairs, `discordantPairs`, zygosity unknown) and the
  Dutch BRCA1 exon 13 deletion of 3,835 bp. Preregistered tests with simulated carriers and 200 real
  non-carriers (`validation/pathogenic/founder/`): 89/90, 29/30, duplication 9/30 frozen and 58/60
  in the registered confirmation; Manta 5/90, 27/30, 37/60; no false calls. The 510 bp exon 22
  deletion is left out (published coordinates fail the Alu check). `spike.py --alu-sequence`.
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

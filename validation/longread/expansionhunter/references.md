# References: genotyping tandem copies, VNTRs and SVs from short reads

Research use only. Collected 2026-11-14 via PubMed; one-line findings are from the abstracts.

1. Dolzhenko E, Deshpande V, Schlesinger F, et al. ExpansionHunter: a sequence-graph-based tool to analyze variation in short tandem repeat regions. Bioinformatics 2019;35(22):4754-4756. doi:10.1093/bioinformatics/btz431
   Sequence-graph genotyping of repeats, including loci with adjacent variants or imperfect repeats; targeted catalogue of loci.
2. Mousavi N, Shleizer-Burko S, Yanicky R, Gymrek M. Profiling the genome-wide landscape of tandem repeat expansions. Nucleic Acids Res 2019;47(15):e90. doi:10.1093/nar/gkz501
   GangSTR: maximum-likelihood TR lengths from paired-end information, genome-wide, including alleles longer than the read.
3. Bakhtiari M, Shleizer-Burko S, Gymrek M, Bansal V, Bafna V. Targeted genotyping of variable number tandem repeats with adVNTR. Genome Res 2018;28(11):1709-1719. doi:10.1101/gr.235119.118
   HMM per VNTR (units 6-100 bp) to count repeat units and detect sequence variation from Illumina or PacBio reads.
4. Lu TY, Human Genome Structural Variation Consortium, Chaisson MJP. Profiling variable-number tandem repeat variation across populations using repeat-pangenome graphs. Nat Commun 2021;12(1):4250. doi:10.1038/s41467-021-24378-0
   danbing-tk: repeat-pangenome graph built from haplotype-resolved assemblies; estimates VNTR composition from short reads; population-stratified VNTR lengths and eQTLs.
5. Chen S, Krusche P, Dolzhenko E, et al. Paragraph: a graph-based structural variant genotyper for short-read sequence data. Genome Biol 2019;20(1):291. doi:10.1186/s13059-019-1909-7
   Genotypes known SVs by aligning reads to a local sequence graph of ref and alt alleles; validated against long-read SV calls.
6. Ebler J, Ebert P, Clarke WE, et al. Pangenome-based genome inference allows efficient and accurate genotyping across a wide spectrum of variant classes. Nat Genet 2022;54(4):518-525. doi:10.1038/s41588-022-01043-w
   PanGenie: k-mer counts plus a haplotype-resolved panel (HMM over panel haplotypes); best gains for insertions >=50 bp and repeats; uses unique k-mers.
7. Eggertsson HP, Kristmundsdottir S, Beyter D, et al. GraphTyper2 enables population-scale genotyping of structural variation using pangenome graphs. Nat Commun 2019;10(1):5402. doi:10.1038/s41467-019-13341-9
   Pangenome-graph SV genotyping at population scale; adding assembly data improves sensitivity for large insertions.
8. Prodanov T, Bansal V. Robust and accurate estimation of paralog-specific copy number for duplicated genes using whole-genome sequencing. Nat Commun 2022;13(1):3221. doi:10.1038/s41467-022-30930-3
   Parascopy: joint analysis of reads across repeat copies; uses paralogous sequence variants (PSVs) to separate copies and multiple samples to reduce bias.
9. Chen X, Shen F, Gonzaludo N, et al. Cyrius: accurate CYP2D6 genotyping using whole-genome sequencing data. Pharmacogenomics J 2021;21(2):251-261. doi:10.1038/s41397-020-00205-5
   Copy number and SV star alleles of CYP2D6 vs its paralog CYP2D7 from depth and copy-specific diagnostic positions; 96.5 % then 99.3 % concordance.
10. Mukamel RE, Handsaker RE, Sherman MA, et al. Protein-coding repeat polymorphisms strongly shape diverse human phenotypes. Science 2021;373(6562):1499-1505. doi:10.1126/science.abg8289
    VNTR lengths estimated from exome reads and imputed into SNP haplotypes; 118 VNTRs in 415,280 UK Biobank participants.
11. Beyter D, Ingimundardottir H, Oddsson A, et al. Long-read sequencing of 3,622 Icelanders provides insight into the role of structural variants in human diseases and other traits. Nat Genet 2021;53(6):779-786. doi:10.1038/s41588-021-00865-4
    133,886 SV alleles from long reads imputed into 166,281 short-read/chip individuals via haplotypes; multiallelic 57 bp repeat in ACAN.
12. Liao WW, Asri M, Ebler J, et al. A draft human pangenome reference. Nature 2023;617(7960):312-324. doi:10.1038/s41586-023-05896-x
    47 phased diploid assemblies; pangenome-based short-read analysis doubles SVs per haplotype vs GRCh38.

No publication on the specific 126 bp tandem copy upstream of BRCA1 exon 1 (chr17:43128247-43128372) was found in PubMed (search on BRCA1 promoter polymorphism/insertion/tandem repeat, 22 hits, none on this allele).

Position of the approaches relative to brca12-wgs:
- Local ref/alt realignment of known alleles: Paragraph (graph alignment), here edit distance of each read against the reference and the allele haplotype (bin/brca12_insertion_signals.py).
- Copy-specific diagnostic bases (PSVs): Parascopy, Cyrius, SMN copy-number callers; here chr17:43128365 C.
- Haplotype prior / imputation of SVs from SNP haplotypes: Mukamel 2021, Beyter 2021, PanGenie panel HMM; here 5 tag SNPs.
- Repeat-unit counting: ExpansionHunter, GangSTR, adVNTR, danbing-tk; here ExpansionHunter v5 as comparator.

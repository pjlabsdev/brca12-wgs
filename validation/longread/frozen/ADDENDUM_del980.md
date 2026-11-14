# Addendum: 980 bp BRCA1 intronic deletion (written before the data were fetched)

Research use only. Written on 2026-11-14 after the main frozen test, before any short reads of the
deletion window were downloaded.

Method: `bin/brca12_insertion_signals.py` at commit 30c053b9fe6b5775c2e336b7c02cb13cfa9489c0, allele
`del980_brca1` (anchor chr17:43071727, 980 bases deleted; breakpoints from 38 of 49 long reads of 5
1KG-ONT carriers), same realignment rules as the main test; zygosity homozygous when reference reads
are at most 10 % of the informative reads.

Samples: the 200 samples of `samples.tsv` (2 of them are long-read carriers) and all other 1KG-ONT
carriers of the deletion with 1000 Genomes 30x data that are not among the 72 samples used before
(`del980_extra_carriers.tsv`, 8 samples; 1 carrier has no 30x data). Short reads: original 1000
Genomes alignments of chr17:43070800-43073700. Truth: carriers of ONT_BRCA1_DEL_43071728 in
`../results/all_disc.json`.

Endpoints: carriers found, non-carriers with a call (exact 95 % confidence intervals), genotype
agreement with the long reads (homozygous when at least 90 % of spanning long reads carry the
deletion). The 6 carriers of the 1000 Genomes short-read panel among the 72 earlier samples (all
found, HG02322 homozygous) were seen before this addendum and are reported separately.

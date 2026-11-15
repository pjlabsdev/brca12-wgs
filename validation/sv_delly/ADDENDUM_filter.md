# Addendum: Delly window filter (written before any test result was looked at)

Committed on 2026-11-15 while the frozen test of `PREREGISTRATION.md` was running and before any of its
output was opened. The pipeline (branch sv-delly) adds one filter after the frozen Delly rule, found on
the real NA18949 reads, which are not a test sample: `bin/brca12_filter_sv.py` applies Manta's call-region
logic to the SV windows (genes ±1 Mb by default):

- POS lies inside a window;
- for DEL, DUP and INV, END also lies inside a window and END − POS is at least 50 bp;
- BND records are kept whatever their partner.

On NA18949 this removed a 31 bp deletion (also called by HaplotypeCaller as a small indel) and an
inversion whose END lay in a distant control window. The evaluation reports Delly with this filter as
the primary Delly arm (it is what the pipeline reports) and the unfiltered frozen rule next to it. The
success criterion is unchanged.

## Harness fix for the real genomes (after the first run)

In the first run Delly wrote no output for any real genome: the region reads fetched from the 1000 Genomes
CRAMs carry @SQ lines for all GRCh38 contigs, while the test reference holds only chr13 and chr17, and
Delly stops on the mismatch. The simulated sets were not affected. `realneg.sh` now keeps only the chr13 and
chr17 @SQ lines before sorting; the real genomes were fetched and analysed again. The method and the
evaluation are unchanged. The pipeline itself uses the full reference and is not affected.

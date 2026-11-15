# Real-data screen: rule for read-level confirmation (fixed before any screen result was seen)

Written 2026-11-14 ~21:56 UTC, before the screen of `screen1.tsv` (340 samples) was started.
Research use only.

A screen hit is a cohort-model segment (frozen step-2 rules: call Q >= 10, hint 5 <= Q < 10) in an
unmodified 1000 Genomes 30x sample. Each hit is checked on the original 1000 Genomes CRAM alignment
in a window from 20 kb upstream of the first to 20 kb downstream of the last exon of the segment
(reads MAPQ >= 20, no duplicates/secondary/supplementary):

1. Depth: mean depth in 500 bp bins over the segment's exons +/- 200 bp, divided by the median of the
   bins in the window outside the segment.
   Deletion: ratio <= 0.70. Duplication: ratio >= 1.30.
2. Breakpoint evidence (at least one):
   a. Deletion: >= 3 read pairs in normal (inward) orientation whose insert exceeds the sample's
      median insert by >= 500 bp, with one mate upstream and one downstream of the depth change.
      Duplication: >= 3 read pairs in outward orientation (first-in-pair maps after its mate on
      reverse/forward strands, i.e. RF) inside the window.
   b. >= 3 soft-clipped reads (clip >= 15 bp) at the same position (+/- 5 bp) within 2 kb of an
      edge of the depth change.

Confirmed = 1 AND 2. Depth only (1 without 2) = "depth-supported, breakpoint not seen" (e.g.
breakpoints in repeats). Neither = "not confirmed" (counted as false call of the screen).
Coordinates are reported as GRCh38; the exact junction is given when clipped reads allow it.

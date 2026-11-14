#!/usr/bin/env python3
"""Insertion signatures from soft-clipped reads inside BRCA1/BRCA2.

A mobile element insertion (for example Alu) leaves a characteristic pattern in short reads: reads
coming from the left are clipped where the foreign sequence starts, reads coming from the right
are clipped at the other junction, and the two clip positions overlap by the target-site
duplication (0 to about 30 bp). SV callers have to assemble the inserted sequence and can fail in
poly-A and repeat context; the pattern itself is simple to count. A deletion gives the opposite
order of clip positions (far apart) and is not reported here.

Input: SAM records on stdin (samtools view of the gene regions). Output: JSON list of candidates.
By default only candidates whose clipped sequence carries Alu consensus k-mers are written, because
paired clip positions without such a signature also occur at alignment artefacts (every tested
sample shows a few). --all keeps the others. Other element families (L1, SVA) are not recognised.
These are hints, not calls.
"""
import argparse
import json
import re
import sys
from collections import Counter, defaultdict

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
COMP = str.maketrans("ACGTN", "TGCAN")
# 12-mers from both ends of the Alu consensus (left monomer start, right monomer end)
ALU_KMERS = ("GGCCGGGCGCGG", "GGCTCACGCCTG", "CTGTAATCCCAG", "GCACTTTGGGAG", "GAGGCGGAGGTT", "GCGAGACTCCGT",
             "CAGCCTGGGCGA", "GGGAGGCTGAGG", "TACTCGGGAGGC", "GGGCGTGGTGGC")
ALU_KMERS = ALU_KMERS + tuple(k.translate(COMP)[::-1] for k in ALU_KMERS)


def scan(lines, min_clip=10, min_mapq=20):
    """Count clip positions. right[pos]: reads aligned up to pos and clipped after it;
    left[pos]: reads clipped before pos and aligned from pos on. Also depth of aligned reads."""
    right, left = defaultdict(list), defaultdict(list)
    starts = defaultdict(list)
    for line in lines:
        if line.startswith("@"):
            continue
        f = line.split("\t", 11)
        if len(f) < 11:
            continue
        flag, chrom, pos, mapq, cigar, seq = int(f[1]), f[2], int(f[3]), int(f[4]), f[5], f[9]
        if flag & 0xF04 or mapq < min_mapq or cigar == "*":      # unmapped, secondary, dup, supplementary, QC fail
            continue
        ops = CIGAR.findall(cigar)
        ref_len = sum(int(n) for n, op in ops if op in "MDN=X")
        end = pos + ref_len - 1
        starts[chrom].append((pos, end))
        if ops[0][1] == "S" and int(ops[0][0]) >= min_clip:
            left[(chrom, pos)].append(seq[:int(ops[0][0])])
        if ops[-1][1] == "S" and int(ops[-1][0]) >= min_clip:
            right[(chrom, end)].append(seq[-int(ops[-1][0]):])
    return right, left, starts


def describe(clips):
    alu = sum(1 for clip in clips if any(k in clip for k in ALU_KMERS))
    poly = sum(1 for clip in clips if re.search(r"A{10,}|T{10,}", clip))
    return alu, poly


def candidates(right, left, starts, regions, min_reads=3, max_tsd=30, min_fraction=0.15, keep_all=False):
    out = []
    left_by_chrom = defaultdict(list)
    for (chrom, pos), clips in left.items():
        if len(clips) >= min_reads:
            left_by_chrom[chrom].append(pos)
    used = set()
    for (chrom, r_pos), r_clips in sorted(right.items()):
        if len(r_clips) < min_reads:
            continue
        gene = next((g for c, s, e, g in regions if c == chrom and s < r_pos <= e), None)
        if not gene:
            continue
        # left-clipped reads start at most max_tsd before, or directly after, the right-clip position
        partners = [p for p in left_by_chrom[chrom] if -1 <= r_pos - p <= max_tsd]
        if not partners:
            continue
        l_pos = max(partners, key=lambda p: len(left[(chrom, p)]))
        if (chrom, l_pos) in used:
            continue
        l_clips = left[(chrom, l_pos)]
        lo, hi = min(l_pos, r_pos), max(l_pos, r_pos)
        depth = sum(1 for s, e in starts[chrom] if s < lo and e > hi) + len(r_clips) + len(l_clips)
        fraction = (len(r_clips) + len(l_clips)) / depth if depth else 0.0
        if fraction < min_fraction:
            continue
        used.add((chrom, l_pos))
        alu, poly = describe(r_clips + l_clips)
        signature = []
        if alu >= 2:
            signature.append("Alu")
        if poly >= 2:
            signature.append("poly-A")
        if "Alu" not in signature and not keep_all:
            continue
        out.append({"gene": gene, "chrom": chrom, "pos": l_pos if l_pos <= r_pos else r_pos, "end": hi,
                    "targetSiteDuplication": max(0, r_pos - l_pos + 1),
                    "clippedReadsLeftJunction": len(r_clips), "clippedReadsRightJunction": len(l_clips),
                    "spanningReads": depth - len(r_clips) - len(l_clips),
                    "clippedFraction": round(fraction, 2),
                    "signature": signature, "clipsWithAluSequence": alu, "clipsWithPolyA": poly})
    return out


def read_regions(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                f = line.rstrip("\n").split("\t")
                rows.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else ""))
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--regions", required=True, help="BED with the two gene regions (contig names as in the BAM)")
    p.add_argument("--out", required=True)
    p.add_argument("--min-reads", type=int, default=3, help="clipped reads needed at each junction")
    p.add_argument("--min-clip", type=int, default=10)
    p.add_argument("--min-fraction", type=float, default=0.15)
    p.add_argument("--all", action="store_true", help="also write candidates without Alu sequence in the clips")
    args = p.parse_args(argv)
    regions = read_regions(args.regions)
    right, left, starts = scan(sys.stdin, args.min_clip)
    found = candidates(right, left, starts, regions, args.min_reads, min_fraction=args.min_fraction, keep_all=args.all)
    with open(args.out, "w") as handle:
        json.dump(found, handle, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())

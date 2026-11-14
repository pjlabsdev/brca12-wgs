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

With --gene-fasta (the reference sequence of the gene regions) a second pattern is searched: an
extra copy of a tandem repeat unit. Reads from the left are clipped at the end of a unit and their
clipped part continues with the start of the same unit (found L bases upstream in the reference);
reads from the right are clipped at the start of a unit and their clipped part matches the end of a
unit L bases downstream. Both junctions with the same length L (50 to 1,000 bp) are required. Such
insertions align to the existing copy with soft-clipped ends only and are missed by Manta and by
the small-variant caller.

Known tandem-copy alleles (KNOWN_COPIES) are also checked with a marker base: the common extra
126 bp copy upstream of BRCA1 exon 1 copies unit 1, which differs from unit 2 at chr17:43128365
(T in unit 2, C in unit 1), so reads of the extra copy placed on unit 2 show C there. Reads of any
mapping quality count (the array is duplicated), base quality at least 13.
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


# chrom, marker position, reference base, base of the extra copy, unit start, unit end, unit length, label
KNOWN_COPIES = (("chr17", 43128365, "T", "C", 43128247, 43128372, 126,
                 "zusätzliche 126-bp-Kopie vor BRCA1 Exon 1 (häufiges Allel)"),)


def count_markers(lines, markers=KNOWN_COPIES, min_baseq=13):
    """Reads with the reference and with the extra-copy base at each marker position."""
    want = {(m[0].replace("chr", ""), m[1]): m for m in markers}
    counts = {m: [0, 0] for m in markers}
    for line in lines:
        if line.startswith("@"):
            continue
        f = line.split("\t", 11)
        if len(f) < 11 or int(f[1]) & 0xF04 or f[5] == "*":
            continue
        chrom, pos = f[2].replace("chr", ""), int(f[3])
        hits = [k for k in want if k[0] == chrom and k[1] >= pos]
        if not hits:
            continue
        ref, rd = pos, 0
        for n, op in CIGAR.findall(f[5]):
            n = int(n)
            if op in "M=X":
                for k in hits:
                    if ref <= k[1] < ref + n:
                        i = rd + k[1] - ref
                        if len(f[10]) <= i or ord(f[10][i]) - 33 >= min_baseq:
                            m = want[k]
                            base = f[9][i].upper()
                            if base == m[2]:
                                counts[m][0] += 1
                            elif base == m[3]:
                                counts[m][1] += 1
                ref += n
                rd += n
            elif op in "DN":
                ref += n
            elif op in "IS":
                rd += n
    return counts


def known_copies(counts, regions, min_reads=2, min_fraction=0.05):
    out = []
    for m, (n_ref, n_alt) in counts.items():
        depth = n_ref + n_alt
        if n_alt < min_reads or not depth or n_alt / depth < min_fraction:
            continue
        chrom = next((c for c, s, e, g in regions if c.replace("chr", "") == m[0].replace("chr", "")), m[0])
        gene = next((g for c, s, e, g in regions if c == chrom and s < m[1] <= e), "")
        out.append({"gene": gene, "chrom": chrom, "pos": m[4], "end": m[5], "duplicatedLength": m[6],
                    "targetSiteDuplication": 0, "clippedReadsLeftJunction": 0, "clippedReadsRightJunction": 0,
                    "spanningReads": depth, "clippedFraction": 0.0, "markerPosition": m[1],
                    "markerReads": f"{n_alt}/{depth}", "label": m[7],
                    "signature": ["tandem copy", "known allele"], "clipsWithAluSequence": 0, "clipsWithPolyA": 0})
    return out


def merge_nearby(clips_by_pos, distance=10):
    """Clip positions within a few bases of a better supported one are counted with it."""
    out = {}
    for key in sorted(clips_by_pos, key=lambda k: -len(clips_by_pos[k])):
        home = next((k for k in out if k[0] == key[0] and abs(k[1] - key[1]) <= distance), None)
        if home:
            out[home] = out[home] + clips_by_pos[key]
        else:
            out[key] = list(clips_by_pos[key])
    return out


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


def read_fasta(path):
    """FASTA of samtools faidx region output (names chrom:start-end) -> {chrom: (start, sequence)}."""
    seqs, name = {}, None
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip().upper())
    out = {}
    for name, parts in seqs.items():
        chrom, _, span = name.rpartition(":")
        if not chrom:
            chrom, span = name, "1-"
        out[chrom] = (int(span.split("-")[0]), "".join(parts))
    return out


def _identity(a, b):
    return sum(x == y for x, y in zip(a, b)) / len(a)


def _place(clip, start, ref, lo, hi, upstream_side, min_identity=0.9):
    """Best reference position p in [lo, hi) where the clipped sequence matches: for a right clip
    ref[p:p+k] (copy start), for a left clip ref[p-k+1:p+1] (copy end)."""
    k, best = len(clip), None
    # nearest copy first, so that identical units give the length of one unit
    for p in (range(hi - 1, lo - 1, -1) if upstream_side else range(lo, hi)):
        i = p - start
        seg = ref[i:i + k] if upstream_side else ref[i - k + 1:i + 1]
        if i - (0 if upstream_side else k - 1) < 0 or len(seg) < k:
            continue
        ident = _identity(clip, seg)
        if ident >= min_identity and (best is None or ident > best[1]):
            best = (p, ident)
    return best


def _lengths(clips_by_pos, fasta, right_side, min_reads, min_clip, min_len, max_len, seed=40):
    """Per clip position: the duplication length most of its clipped reads point to."""
    out = []
    for (chrom, pos), clips in clips_by_pos.items():
        clips = [c for c in clips if len(c) >= min_clip]
        if len(clips) < min_reads or chrom not in fasta:
            continue
        start, ref = fasta[chrom]
        lengths = []
        for clip in clips:
            if right_side:       # clipped part follows pos: it starts the next copy, L bases upstream
                hit = _place(clip[:seed], start, ref, pos - max_len + 1, pos - min_len + 2, True)
                if hit:
                    lengths.append(pos - hit[0] + 1)
            else:                # clipped part precedes pos: it ends the previous copy, L bases downstream
                hit = _place(clip[-seed:], start, ref, pos + min_len - 1, pos + max_len, False)
                if hit:
                    lengths.append(hit[0] - pos + 1)
        if len(lengths) >= min_reads:
            length, n = Counter(lengths).most_common(1)[0]
            support = sum(1 for x in lengths if abs(x - length) <= 0.1 * length)
            if support >= min_reads:
                out.append((chrom, pos, length, support))
    return out


def duplications(right, left, starts, regions, fasta, min_reads=2, min_clip=10, min_len=50, max_len=1000):
    """Tandem copy signatures: a right-clip and a left-clip position whose clips point to the same length."""
    rights = _lengths(merge_nearby(right), fasta, True, min_reads, min_clip, min_len, max_len)
    lefts = _lengths(merge_nearby(left), fasta, False, min_reads, min_clip, min_len, max_len)
    out, used = [], set()
    for chrom, b, length, n_r in sorted(rights, key=lambda x: -x[3]):
        gene = next((g for c, s, e, g in regions if c == chrom and s < b <= e), None)
        if not gene:
            continue
        partners = [x for x in lefts if x[0] == chrom and abs(x[2] - length) <= 0.1 * length
                    and abs(b - x[1]) <= 2 * length and (x[0], x[1]) not in used]
        if not partners:
            continue
        _, a, _, n_l = max(partners, key=lambda x: x[3])
        used.add((chrom, a))
        lo, hi = min(a, b - length + 1), max(a + length - 1, b)
        spanning = sum(1 for s, e in starts[chrom] if s < b - 10 and e > b + 10)   # reads across the junction
        out.append({"gene": gene, "chrom": chrom, "pos": lo, "end": hi, "duplicatedLength": length,
                    "targetSiteDuplication": 0, "clippedReadsLeftJunction": n_r, "clippedReadsRightJunction": n_l,
                    "spanningReads": spanning,
                    "clippedFraction": round((n_r + n_l) / (spanning + n_r + n_l), 2) if spanning + n_r + n_l else 0.0,
                    "signature": ["tandem copy"], "clipsWithAluSequence": 0, "clipsWithPolyA": 0})
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
    p.add_argument("--gene-fasta", help="reference sequence of the gene regions (samtools faidx output); "
                                        "enables the search for extra tandem repeat copies")
    args = p.parse_args(argv)
    regions = read_regions(args.regions)
    lines = sys.stdin.readlines()
    right, left, starts = scan(lines, args.min_clip)
    found = candidates(right, left, starts, regions, args.min_reads, min_fraction=args.min_fraction, keep_all=args.all)
    if args.gene_fasta:
        dups = duplications(right, left, starts, regions, read_fasta(args.gene_fasta))
        for known in known_copies(count_markers(lines), regions):
            same = [d for d in dups if d["chrom"] == known["chrom"] and abs(d["end"] - known["end"]) <= 20
                    and abs(d["duplicatedLength"] - known["duplicatedLength"]) <= 10]
            if same:
                same[0].update(markerPosition=known["markerPosition"], markerReads=known["markerReads"],
                               label=known["label"], signature=["tandem copy", "known allele"])
            else:
                dups.append(known)
        found += dups
    with open(args.out, "w") as handle:
        json.dump(found, handle, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())

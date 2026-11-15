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

A third, independent signal is the haplotype the extra copy sits on (HAPLOTYPE_TAGS): SNPs 3 to 21 kb
downstream that are almost always inherited together with it (linkage in the HPRC assemblies and the
1KG-ONT long reads). The known allele is reported when at least 2 of the 3 signals agree: clipped
reads at both unit borders, the marker base, and at least 3 of the tag SNPs (alternative base in at
least 3 reads and 20 % of the reads, mapping quality at least 20). When at least 3 tag SNPs are
homozygous (80 %), two copies of the haplotype are likely. A tandem-copy signature from clipped
reads alone is still reported.

Known alleles from long-read assemblies (KNOWN_ALLELES: the 126 bp copy and the common intronic BRCA2
Alu insertion, sequences from the HPRC pangenome) are also genotyped by allele-aware realignment:
every read near the site (unmapped mates placed there included) is aligned to the reference
haplotype and to the haplotype carrying the allele (edit distance, Myers bit-vector algorithm,
free start and end in the haplotype). A read supports the allele when it fits that haplotype with at
least 2 fewer edits, and the reference when the opposite holds; reads with more than 8 % edits to
both are ignored. The allele is present when at least 2 reads support it. For an insertion that is
not a tandem copy, reference-supporting reads across the site give the zygosity; for the tandem
copy no short read can tell the alleles apart in that way, and the zygosity comes from the haplotype.
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
# id, chrom, anchor position (VCF POS), inserted sequence, tandem copy (no reference-specific reads), label,
# deleted bases after the anchor (optional, 0 for insertions)
KNOWN_ALLELES = (
    ("tandem126", "chr17", 43128264,
     "GATGGGAGGAACGAACAAACTCCAGACATGCCGCCTTAAGAGCTGTAACACTCACCGCAAAGGTCTGCAGCTTCACTCCTGAAGCCAGCGAGACCACGAACCCACC"
     "AGGAGGAATGAACAACTCCA", True, "zusätzliche 126-bp-Kopie vor BRCA1 Exon 1 (häufiges Allel)"),
    ("alu_brca2", "chr13", 32348295,
     "AAAAGTTTTGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCCGGCTAAAAC"
     "GGTGAAACCCCGTCTCTACTAAAAATACAAAAAAATTAGCCGGGCGTAGTGGCGGGCGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGTGTGAAC"
     "CCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTC" + "A" * 27, False,
     "bekannte intronische Alu-Insertion in BRCA2 (häufiges Allel)"),
    # 980 bp deletion in BRCA1 intron (1KG-ONT: 13 of 993 samples; breakpoints from 38 of 49 long reads
    # of 5 carriers; listed as 1,081 bp in the 1000 Genomes short-read panel)
    ("del980_brca1", "chr17", 43071727, "", False, "bekannte intronische 980-bp-Deletion in BRCA1 (seltenes Allel)", 980),
    # BRCA2 c.156_157insAlu, pathogenic founder allele (Portugal; Teugels et al. 2005, Hum Mutat,
    # Mutation in Brief 837; Peixoto et al. 2009): AluYa5 in the same orientation as BRCA2, flanked by the 14 bp target site
    # duplication AAGAATCTGAACAT (chr13:32319152-32319165, ending at c.156). The exact element sequence is not
    # published; the Dfam AluYa5 consensus with its poly-A tail stands in for it.
    ("alu_brca2_c156", "chr13", 32319165,
     "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCCGGCTAAAACGG"
     "TGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTAGTGGCGGGCGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCGT"
     "GAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAA"
     "AAAAAAAAAAA"
     "AAGAATCTGAACAT", False, "BRCA2 c.156_157insAlu (Gründermutation, Exon 3)"),
)
# known alleles that are pathogenic (the report flags them for laboratory confirmation)
# BRCA1 exon 12 (legacy exon 13) duplication, pathogenic founder allele (Northern Europe and elsewhere;
# The BRCA1 Exon 13 Duplication Screening Group 2000; Edwinsdotter Ardnor et al. 2019:
# NM_007294 c.4186-1787_4357+4122dup). On GRCh38 c.4186 = chr17:43082575 and c.4357 = chr17:43082404
# (exon boundaries, ../validation/pathogenic/cds.py), so the duplicated segment is chr17:43078282-43084362
# (6,081 bp). Both ends lie in antisense Alu elements with long homology, as expected for the published
# Alu-Alu recombination. Inserted sequence: ("dup", start, end) = the reference segment, after its end.
KNOWN_ALLELES = KNOWN_ALLELES + (
    ("dup13_brca1", "chr17", 43084362, ("dup", 43078282, 43084362), False,
     "BRCA1 Duplikation Exon 13 (alte Zählung; Gründermutation)"),
)
# BRCA1 exon 12 (legacy exon 13) deletion of 3,835 bp, Dutch founder allele (Petrij-Bosch et al. 1997;
# "IVS12-1643del3835" in Zeegers et al. 2004, = c.4186-1643_4357+2020del). On GRCh38 the deleted segment
# is chr17:43080384-43084218; both ends lie in antisense Alu elements, as described for this
# Alu-mediated deletion.
KNOWN_ALLELES = KNOWN_ALLELES + (
    ("del3835_brca1", "chr17", 43080383, "", False, "BRCA1 Deletion Exon 13 (alte Zählung; 3.835 bp, Gründermutation)", 3835),
)
PATHOGENIC = frozenset({"alu_brca2_c156", "dup13_brca1", "del3835_brca1"})
# realigned tandem-copy allele -> marker position of the same copy in KNOWN_COPIES
ALLELE_MARKER = {"tandem126": 43128365}
# marker position of a known copy -> SNPs on the same haplotype (chrom, position, reference, alternative)
HAPLOTYPE_TAGS = {43128365: (("chr17", 43131837, "G", "A"), ("chr17", 43135863, "A", "G"),
                             ("chr17", 43137747, "A", "G"), ("chr17", 43138596, "G", "T"),
                             ("chr17", 43148782, "C", "A"))}


def count_markers(lines, markers=KNOWN_COPIES, min_baseq=13, min_mapq=0):
    """Reads with the reference and with the extra-copy base at each marker position."""
    want = {(m[0].replace("chr", ""), m[1]): m for m in markers}
    counts = {m: [0, 0] for m in markers}
    for line in lines:
        if line.startswith("@"):
            continue
        f = line.split("\t", 11)
        if len(f) < 11 or int(f[1]) & 0xF04 or f[5] == "*" or int(f[4]) < min_mapq:
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


def haplotype(tag_counts, min_reads=3, min_fraction=0.2, hom_fraction=0.8, min_tags=3):
    """(tags carrying the alternative base, tags tested, copies of the haplotype: 0, 1 or 2)."""
    alt = [k for k, (r, a) in tag_counts.items() if a >= min_reads and a / (r + a) >= min_fraction]
    hom = [k for k in alt if tag_counts[k][1] / sum(tag_counts[k]) >= hom_fraction]
    copies = 2 if len(hom) >= min_tags else (1 if len(alt) >= min_tags else 0)
    return len(alt), len(tag_counts), copies


def known_copies(counts, regions, min_reads=2, min_fraction=0.05, haplotypes=None, clipped=(), realigned=None):
    """Known alleles found by allele realignment (at least min_reads supporting reads) or supported by
    at least 2 of: marker base, haplotype tags, clipped reads. haplotypes: {marker position: haplotype()
    result}; clipped: marker positions with a clip signature; realigned: {marker position: (ref, allele)}."""
    haplotypes, realigned = haplotypes or {}, realigned or {}
    out = []
    for m, (n_ref, n_alt) in counts.items():
        depth = n_ref + n_alt
        marker = n_alt >= min_reads and depth and n_alt / depth >= min_fraction
        hap = haplotypes.get(m[1])
        evidence = [name for name, ok in (("clipped reads", m[1] in clipped), ("marker base", marker),
                                          ("haplotype", bool(hap and hap[2]))) if ok]
        n_realigned = realigned.get(m[1], (0, 0))[1]
        if len(evidence) < 2 and n_realigned < min_reads:
            continue
        if n_realigned >= min_reads:
            evidence.insert(0, "allele realignment")
        chrom = next((c for c, s, e, g in regions if c.replace("chr", "") == m[0].replace("chr", "")), m[0])
        gene = next((g for c, s, e, g in regions if c == chrom and s < m[1] <= e), "")
        out.append({"gene": gene, "chrom": chrom, "pos": m[4], "end": m[5], "duplicatedLength": m[6],
                    "targetSiteDuplication": 0, "clippedReadsLeftJunction": 0, "clippedReadsRightJunction": 0,
                    "spanningReads": depth, "clippedFraction": 0.0, "markerPosition": m[1],
                    "markerReads": f"{n_alt}/{depth}", "label": m[7],
                    "signature": ["tandem copy", "known allele"], "clipsWithAluSequence": 0, "clipsWithPolyA": 0,
                    "evidence": evidence})
        if m[1] in realigned:
            out[-1]["realignedReads"] = f"{n_realigned} allele"
        if hap:
            out[-1].update(haplotypeTags=f"{hap[0]}/{hap[1]}", haplotypeCopies=hap[2])
    return out


def infix_distance(pattern, text):
    """Smallest edit distance of pattern to any substring of text (Myers 1999 bit-vector algorithm)."""
    m = len(pattern)
    if not m:
        return 0
    peq = {}
    for i, c in enumerate(pattern):
        peq[c] = peq.get(c, 0) | (1 << i)
    mask, high = (1 << m) - 1, 1 << (m - 1)
    pv, mv, score, best = mask, 0, m, m
    for c in text:
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (~(xh | pv) & mask)
        mh = pv & xh
        if ph & high:
            score += 1
        elif mh & high:
            score -= 1
        ph = (ph << 1) & mask          # free start in the text: no carry-in
        mh = (mh << 1) & mask
        pv = mh | (~(xv | ph) & mask)
        mv = ph & xv
        best = min(best, score)
    return best


def realign_known(lines, fasta, alleles=KNOWN_ALLELES, flank=700, before=400, after=300, min_diff=2,
                  max_error=0.08, min_len=50):
    """{allele id: (reads supporting the reference, reads supporting the allele)} by local realignment."""
    comp = str.maketrans("ACGTN", "TGCAN")
    sites = []
    for a in alleles:
        chrom = next((c for c in fasta if c.replace("chr", "") == a[1].replace("chr", "")), None)
        if chrom is None:
            continue
        start, seq = fasta[chrom]
        i = a[2] - start                          # index of the base after the anchor
        deleted = a[6] if len(a) > 6 else 0
        if i - flank < 0 or i + deleted + flank > len(seq):
            continue
        if isinstance(a[3], str):
            inserted, left = a[3], i - flank
        elif a[3][1] - start - flank < 0:
            continue
        else:                                     # ("dup", start, end): copy of a reference segment
            inserted = seq[a[3][1] - start:a[3][2] - start + 1]
            left = a[3][1] - start - flank        # both haplotypes span the whole segment, so every
                                                  # reference junction is in both and only the new one differs
        ref = seq[left:i + deleted + flank]
        sites.append((a, chrom.replace("chr", ""), ref, seq[left:i + 1] + inserted + seq[i + 1 + deleted:i + deleted + flank]))
    counts = {a[0]: [0, 0] for a, *_ in sites}
    seen = set()
    for line in lines:
        if line.startswith("@"):
            continue
        f = line.split("\t", 11)
        if len(f) < 11 or int(f[1]) & 0xF00 or len(f[9]) < min_len:
            continue
        chrom, pos = f[2].replace("chr", ""), int(f[3])
        for a, c, ref, alt in sites:
            deleted = a[6] if len(a) > 6 else 0
            # a duplication junction read can be placed at either end of the duplicated segment
            dup_start = a[3][1] if not isinstance(a[3], str) else 0
            if c != chrom or not (a[2] - before < pos < a[2] + after or
                                  deleted and a[2] + deleted - before < pos < a[2] + deleted + after or
                                  dup_start and dup_start - before < pos < dup_start + after):
                continue
            key = (a[0], f[0], int(f[1]) & 0xC0)
            if key in seen:
                continue
            seen.add(key)
            read = f[9].upper()
            reads = (read, read.translate(comp)[::-1]) if int(f[1]) & 4 else (read,)
            d_ref, d_alt = min(((infix_distance(r, ref), infix_distance(r, alt)) for r in reads), key=min)
            if min(d_ref, d_alt) > max_error * len(read):
                continue
            if d_alt + min_diff <= d_ref:
                counts[a[0]][1] += 1
            elif d_ref + min_diff <= d_alt:
                counts[a[0]][0] += 1
    return counts


def duplication_pairs(lines, alleles=KNOWN_ALLELES, reach=700, min_mapq=20):
    """{allele id: read pairs that face outward across a known tandem duplication}.

    A fragment across the new junction (segment end joined to segment start) has its reverse read near
    the segment start and its forward mate near the segment end, the opposite of a normal pair."""
    sites = [(a[0], a[1].replace("chr", ""), a[3][1], a[3][2]) for a in alleles if not isinstance(a[3], str)]
    names = {s[0]: set() for s in sites}
    for line in lines:
        if line.startswith("@"):
            continue
        f = line.split("\t", 10)
        if len(f) < 10:
            continue
        flag = int(f[1])
        if flag & 0xF0C or int(f[4]) < min_mapq or f[6] != "=" or not flag & 16 or flag & 32:
            continue
        chrom, pos, mpos = f[2].replace("chr", ""), int(f[3]), int(f[7])
        for allele, c, start, end in sites:
            if c == chrom and start - 50 <= pos <= start + reach and end - reach - 150 <= mpos <= end:
                names[allele].add(f[0])
    return {a: len(n) for a, n in names.items()}


def outward_pairs(lines, regions, min_mapq=20, window=600):
    """Clusters of read pairs that face outward (reverse read left of its forward mate), the read-pair
    signature of a tandem duplication: [(chrom, first reverse-read start, last forward-mate end, names)].
    Clipped reads whose supplementary alignment (SA tag) joins the two ends are added to the cluster."""
    pairs, clips = [], []
    seen = set()
    for line in lines:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        flag = int(f[1])
        if flag & 0xF0C or int(f[4]) < min_mapq:
            continue
        chrom, pos = f[2], int(f[3])
        if not any(c == chrom and s < pos <= e for c, s, e, g in regions):
            continue
        if flag & 1 and f[6] == "=" and flag & 16 and not flag & 32 and int(f[7]) > pos and f[0] not in seen:
            seen.add(f[0])
            pairs.append((chrom, pos, int(f[7]) + len(f[9]), f[0]))
        sa = next((x[5:] for x in f[11:] if x.startswith("SA:Z:")), None)
        if sa:
            cig = CIGAR.findall(f[5])
            aligned = sum(int(n) for n, op in cig if op in "MD=X")
            for entry in sa.rstrip(";").split(";"):
                parts = entry.split(",")
                if len(parts) >= 5 and parts[0] == chrom and int(parts[4]) >= min_mapq:
                    other = int(parts[1])
                    # left part aligned at the segment end, right part at the segment start (or the reverse)
                    if cig and cig[-1][1] == "S" and other < pos:
                        clips.append((chrom, other, pos + aligned, f[0]))
                    elif cig and cig[0][1] == "S" and other > pos:
                        clips.append((chrom, pos, other + aligned, f[0]))
    clusters = []
    for chrom, start, end, name in sorted(pairs):
        home = next((c for c in clusters if c[0] == chrom and abs(c[1] - start) <= window and
                     abs(c[2] - end) <= window), None)
        if home:
            home[1], home[2] = min(home[1], start), max(home[2], end)
            home[3].add(name)
        else:
            clusters.append([chrom, start, end, {name}, set()])
    for chrom, start, end, name in clips:
        home = next((c for c in clusters if c[0] == chrom and abs(c[1] - start) <= window and
                     abs(c[2] - end) <= window), None)
        if home:
            home[4].add(name)
    return clusters


def tandem_duplications(lines, regions, exons=None, min_support=4, min_span=1000, max_span=200000):
    """Report entries for tandem duplications found by outward-facing read pairs plus joining split reads."""
    out = []
    for chrom, start, end, names, split in outward_pairs(lines, regions):
        support = len(names | split)
        if support < min_support or not min_span <= end - start <= max_span:
            continue
        gene = next((g for c, s, e, g in regions if c == chrom and s < start <= e), "")
        out.append({"gene": gene, "chrom": chrom, "pos": start, "end": end, "duplicatedLength": end - start,
                    "signature": ["tandem duplication", "read pairs"], "evidence": ["read pairs"] +
                    (["split reads"] if split else []), "outwardPairs": len(names), "splitReads": len(split),
                    "breakpointPrecision": "approximate (read-pair span)",
                    "clippedReadsLeftJunction": 0, "clippedReadsRightJunction": 0})
    return out


def _insert_stats(tlens):
    """Median and robust spread of the insert size of properly paired reads."""
    v = sorted(tlens)
    if len(v) < 100:
        return 450.0, 100.0
    med = v[len(v) // 2]
    sample = v[:: max(1, len(v) // 5000)]
    mad = sorted(abs(x - med) for x in sample)[len(sample) // 2]
    return float(med), max(1.4826 * mad, 20.0)


def sv_signals(lines, regions, min_mapq=20):
    """Read pairs and split reads that point to deletions or inversions inside the gene regions.
    Deletion pairs: forward read, reverse mate further right, insert above median + 6 spreads.
    Inversion pairs: read and mate on the same strand, mate at least 300 bp further right.
    Split reads: supplementary alignment (SA tag) on the same chromosome; same strand with a gap in
    reference order is a deletion junction, opposite strand an inversion junction."""
    reads, tlens = [], []
    for line in lines:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        flag = int(f[1])
        if flag & 0xF0C or int(f[4]) < min_mapq:
            continue
        chrom, pos = f[2], int(f[3])
        if not any(c == chrom and s < pos <= e for c, s, e, g in regions):
            continue
        cig = CIGAR.findall(f[5])
        span = sum(int(n) for n, op in cig if op in "MD=X")
        reads.append((f, flag, chrom, pos, cig, span))
        if flag & 2 and f[6] == "=" and 0 < int(f[8]) < 2000:
            tlens.append(int(f[8]))
    med, sd = _insert_stats(tlens)
    cut = med + 6 * sd
    dpairs, ipairs, dsplit, isplit = [], [], [], []
    for f, flag, chrom, pos, cig, span in reads:
        name = f[0]
        if flag & 1 and f[6] == "=" and not flag & 8:
            mpos = int(f[7])
            rev, mrev = bool(flag & 16), bool(flag & 32)
            if mpos > pos:
                if not rev and mrev and mpos + len(f[9]) - pos > cut:
                    dpairs.append((chrom, pos + span, mpos, name))
                elif rev == mrev and mpos - pos > 300:
                    ipairs.append((chrom, pos, mpos, name, "--" if rev else "++"))
        sa = next((x[5:] for x in f[11:] if x.startswith("SA:Z:")), None)
        if not sa or not cig:
            continue
        strand = "-" if flag & 16 else "+"
        for entry in sa.rstrip(";").split(";"):
            p = entry.split(",")
            if len(p) < 6 or p[0] != chrom or int(p[4]) < min_mapq:
                continue
            opos, ocig = int(p[1]), CIGAR.findall(p[3])
            ospan = sum(int(n) for n, op in ocig if op in "MD=X")
            if not ocig:
                continue
            if p[2] == strand:
                if cig[-1][1] in "SH" and opos > pos + span + 30 and ocig[0][1] in "SH":
                    dsplit.append((chrom, pos + span, opos, name))
                elif cig[0][1] in "SH" and opos + ospan + 30 < pos and ocig[-1][1] in "SH":
                    dsplit.append((chrom, opos + ospan, pos, name))
            else:
                def junction(start, c, length):
                    right = c[-1][1] in "SH" and not (c[0][1] in "SH" and int(c[0][0]) > int(c[-1][0]))
                    return start + length if right else start
                a, b = sorted((junction(pos, cig, span), junction(opos, ocig, ospan)))
                if b - a > 300:
                    isplit.append((chrom, a, b, name))
    return dpairs, ipairs, dsplit, isplit, (med, sd)


def _clusters(items, window):
    out = []
    for it in sorted(items):
        home = next((c for c in out if c[0] == it[0] and abs(c[1] - it[1]) <= window and
                     abs(c[2] - it[2]) <= window), None)
        if home:
            home[3].append(it)
        else:
            out.append([it[0], it[1], it[2], [it]])
    return out


def _sv_entry(regions, chrom, start, end, signature, pairs, split, precise, kind_len):
    gene = next((g for c, s, e, g in regions if c == chrom and s < start <= e), "")
    return {"gene": gene, "chrom": chrom, "pos": start, "end": end, kind_len: end - start + 1,
            "signature": signature, "evidence": (["read pairs"] if pairs else []) + (["split reads"] if split else []),
            "discordantPairs": pairs, "splitReads": split,
            "breakpointPrecision": "exact (split reads)" if precise else "approximate (read-pair span)",
            "clippedReadsLeftJunction": 0, "clippedReadsRightJunction": 0}


def read_pair_svs(lines, regions, min_support=4, min_del=50, min_inv=300, max_size=200000):
    """Deletions and inversions from discordant read pairs and split reads (distinct read names >= min_support)."""
    dpairs, ipairs, dsplit, isplit, (med, sd) = sv_signals(lines, regions)
    w = int(med + 6 * sd)
    out, used = [], set()
    for chrom, a, b, its in _clusters(dsplit, 20):
        (s, e), _ = Counter((x[1], x[2]) for x in its).most_common(1)[0]
        names = {x[3] for x in its}
        pnames = {x[3] for x in dpairs if x[0] == chrom and s - w <= x[1] <= s + 20 and e - 20 <= x[2] <= e + w}
        if len(names | pnames) >= min_support and min_del <= e - s <= max_size:
            out.append(_sv_entry(regions, chrom, s, e - 1, ["read-pair deletion"], len(pnames), len(names), True,
                                 "deletedLength"))
            used |= pnames
    for chrom, a, b, its in _clusters([x for x in dpairs if x[3] not in used], w):
        names = {x[3] for x in its}
        s, e = max(x[1] for x in its), min(x[2] for x in its)
        if len(names) >= min_support and min_del <= e - s <= max_size and not any(
                o["chrom"] == chrom and abs(o["pos"] - s) < w and abs(o["end"] - e) < w for o in out):
            out.append(_sv_entry(regions, chrom, s, e - 1, ["read-pair deletion"], len(names), 0, False,
                                 "deletedLength"))
    wi = int(med + 5 * sd)
    invs = []
    for chrom, a, b, its in _clusters(ipairs, wi):
        names = {x[3] for x in its}
        pp, mm = [x for x in its if x[4] == "++"], [x for x in its if x[4] == "--"]
        sp = [x for x in isplit if x[0] == chrom and abs(x[1] - a) <= wi and abs(x[2] - b) <= wi]
        A, B = a, b
        if pp and mm:
            A = (max(x[1] for x in pp) + 150 + min(x[1] for x in mm)) // 2
            B = (max(x[2] for x in pp) + 150 + min(x[2] for x in mm)) // 2
        if sp:
            A, B = sorted(x[1] for x in sp)[len(sp) // 2], sorted(x[2] for x in sp)[len(sp) // 2]
        if any(o["chrom"] == chrom and abs(o["pos"] - A) <= wi and abs(o["end"] - B) <= wi for o in invs):
            continue
        spn = {x[3] for x in sp}
        if len(names | spn) >= min_support and min_inv <= B - A <= max_size:
            invs.append(_sv_entry(regions, chrom, A, B - 1, ["inversion", "read pairs"], len(names), len(spn),
                                  bool(sp), "invertedLength"))
    return out + invs


def realigned_alleles(counts, regions, haplotypes=None, min_reads=2, pairs=None, min_pairs=3):
    """Report entries for known alleles found by realignment (non-tandem insertions and deletions)."""
    out = []
    for a in KNOWN_ALLELES:
        if a[0] not in counts or a[4]:
            continue
        n_ref, n_alt = counts[a[0]]
        n_pairs = (pairs or {}).get(a[0], 0)
        # duplication: the Alu homology at the junction leaves few informative reads, so outward-facing
        # pairs plus one junction read also count
        if n_alt < min_reads and not (not isinstance(a[3], str) and n_pairs >= min_pairs and n_alt >= 1):
            continue
        chrom = next((c for c, s, e, g in regions if c.replace("chr", "") == a[1].replace("chr", "")), a[1])
        gene = next((g for c, s, e, g in regions if c == chrom and s < a[2] <= e), "")
        deleted = a[6] if len(a) > 6 else 0
        kind = ["Alu", "known allele"] if a[0].startswith("alu") else ["known allele"]
        if deleted:
            kind.append("deletion")
        duplicated = 0 if isinstance(a[3], str) else a[3][2] - a[3][1] + 1
        if duplicated:
            kind.append("duplication")
        out.append({"gene": gene, "chrom": chrom, "pos": a[3][1] if duplicated else a[2] + (1 if deleted else 0),
                    "end": a[3][2] if duplicated else a[2] + max(deleted, 1),
                    "insertedLength": duplicated or len(a[3]), "deletedLength": deleted,
                    "targetSiteDuplication": 0, "clippedReadsLeftJunction": 0, "clippedReadsRightJunction": 0,
                    "spanningReads": n_ref, "clippedFraction": 0.0, "label": a[5],
                    "signature": kind,
                    "clipsWithAluSequence": 0, "clipsWithPolyA": 0,
                    "evidence": ["allele realignment"] + (["read pairs"] if n_pairs >= min_pairs else []),
                    "realignedReads": f"{n_alt} allele / {n_ref} reference",
                    **({"discordantPairs": n_pairs} if not isinstance(a[3], str) else {}),
                    "clinicalSignificance": "pathogenic" if a[0] in PATHOGENIC else "",
                    # a tandem duplication keeps every reference junction, so no read is reference-specific
                    "zygosity": "unknown" if duplicated else
                    "homozygous" if n_ref <= 0.1 * (n_ref + n_alt) and n_alt >= 2 * min_reads
                    else "heterozygous"})
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
    p.add_argument("--read-pair-svs", action="store_true",
                   help="also report deletions and inversions from discordant read pairs (research arm, off by default)")
    p.add_argument("--gene-fasta", help="reference sequence of the gene regions (samtools faidx output); "
                                        "enables the search for extra tandem repeat copies")
    args = p.parse_args(argv)
    regions = read_regions(args.regions)
    lines = sys.stdin.readlines()
    right, left, starts = scan(lines, args.min_clip)
    found = candidates(right, left, starts, regions, args.min_reads, min_fraction=args.min_fraction, keep_all=args.all)
    if args.gene_fasta:
        fasta = read_fasta(args.gene_fasta)
        dups = duplications(right, left, starts, regions, fasta)

        def same_copy(d, m):
            return (d["chrom"].replace("chr", "") == m[0].replace("chr", "") and abs(d["end"] - m[5]) <= 20
                    and abs(d["duplicatedLength"] - m[6]) <= 10)

        clipped = {m[1] for m in KNOWN_COPIES if any(same_copy(d, m) for d in dups)}
        haps = {pos: haplotype(count_markers(lines, tags, min_mapq=20)) for pos, tags in HAPLOTYPE_TAGS.items()}
        realigned = realign_known(lines, fasta)
        by_marker = {pos: tuple(realigned[a]) for a, pos in ALLELE_MARKER.items() if a in realigned}
        for known in known_copies(count_markers(lines), regions, haplotypes=haps, clipped=clipped,
                                  realigned=by_marker):
            m = next(m for m in KNOWN_COPIES if m[1] == known["markerPosition"])
            same = [d for d in dups if same_copy(d, m)]
            if same:
                same[0].update({k: v for k, v in known.items() if k in (
                    "markerPosition", "markerReads", "label", "signature", "evidence", "haplotypeTags",
                    "haplotypeCopies", "realignedReads")})
            else:
                dups.append(known)
        found += dups
        # known non-tandem alleles: replace a matching Alu read-end signature, else add the entry
        for known in realigned_alleles(realigned, regions, pairs=duplication_pairs(lines)):
            same = [d for d in found if d["chrom"] == known["chrom"] and "tandem copy" not in d["signature"]
                    and "deletion" not in known["signature"] and "duplication" not in known["signature"]
                    and abs(d["pos"] - known["pos"]) <= 30]
            if same:
                same[0].update({k: known[k] for k in ("label", "evidence", "realignedReads", "zygosity",
                                                      "insertedLength", "clinicalSignificance")})
                same[0]["evidence"] = ["allele realignment", "clipped reads"]
                same[0]["signature"] = sorted(set(same[0]["signature"]) | {"known allele"})
            else:
                found.append(known)
        # tandem duplications from outward-facing read pairs, unless a known allele already covers them
        for dup in tandem_duplications(lines, regions):
            if not any(d["chrom"] == dup["chrom"] and "duplication" in d.get("signature", []) and
                       abs(d["pos"] - dup["pos"]) <= 600 and abs(d["end"] - dup["end"]) <= 600 for d in found):
                found.append(dup)
    # research arm, off by default (Delly is the structural-variant caller): deletions and inversions from
    # discordant pairs and split reads, unless a known allele covers them
    if args.read_pair_svs:
        for sv in read_pair_svs(lines, regions):
            if not any(d["chrom"] == sv["chrom"] and "deletion" in d.get("signature", []) and
                       abs(d["pos"] - sv["pos"]) <= 600 and abs(d["end"] - sv["end"]) <= 600 for d in found):
                found.append(sv)
    with open(args.out, "w") as handle:
        json.dump(found, handle, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())

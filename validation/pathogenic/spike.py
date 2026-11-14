#!/usr/bin/env python3
"""Spike-in simulator: known heterozygous variants in real WGS background reads.

Idea: inside the two gene windows (gene +/- 20 kb) every real read pair of the background sample
is replaced by a pair with the same position, fragment length, read lengths and base qualities,
but with its sequence read from one of two haplotypes: A (the reference) or B (the reference with
the spiked events). Coverage profile, fragment sizes and base-quality profile therefore stay those
of the real run; sequencing errors are drawn from the real base qualities. Outside the windows the
real reads are kept untouched (they give CNVkit its chromosome-wide baseline).

The spiked events are synthetic (some are modelled on well-known pathogenic variants); the
background sample's own variants inside the windows are not carried over. This measures whether a
variant that is present in the reads is called and reported. It does not replace real samples.

Sub-commands
  design    write the small-variant spike lists (deterministic, seeded)
  prepare   collect the window fragments and write the background BAM without them
  simulate  write FASTQ + truth files for one scenario sample
"""
import argparse
import gzip
import hashlib
import json
import os
import random
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cds import Gene, fetch, revcomp  # noqa: E402

WINDOWS = {"BRCA1": ("chr17", 43024295, 43190327), "BRCA2": ("chr13", 32295077, 32420268)}
# Alu element copied from the reference (BRCA1 intron, forward strand, found by the Alu 5' consensus
# motif and its A-rich tail); inserted with a 14 bp target-site duplication.
ALU_SOURCE = ("chr17", 43152573, 43152885)
TSD = 14
SMALL_TYPES = ["SNV", "DEL1", "INS1", "DEL2", "SNV", "DEL4", "INS3", "DUP5", "DEL11", "SNV", "INS8",
               "DEL19", "INS25", "DEL40"]
# Modelled on well-known pathogenic variants (HGVS on the MANE Select transcript); the reference
# bases are checked against the CDS model before use.
KNOWN = {
    "SPIKE_S1": [("BRCA1", "c.68_69del", 68, "AG", ""), ("BRCA1", "c.181T>G", 181, "T", "G"),
                 ("BRCA1", "c.5266dup", 5266, "C", "CC"), ("BRCA2", "c.771_775del", 771, "TCAAA", ""),
                 ("BRCA2", "c.5946del", 5946, "T", "")],
    "SPIKE_S2": [],
}


def stable_random(*parts):
    return int(hashlib.sha1("|".join(map(str, parts)).encode()).hexdigest()[:12], 16) / float(1 << 48)


# ---------------------------------------------------------------- design
def known_event(gene, name, c, ref, alt, fasta):
    """HGVS-like (c position, transcript ref, transcript alt) -> VCF-style genomic event."""
    if gene.cds[c - 1:c - 1 + len(ref)] != ref:
        raise ValueError(f"{gene.name} {name}: reference bases do not match the CDS")
    positions = [gene.genomic(c + i) for i in range(len(ref))]
    lo, hi = min(positions), max(positions)
    if hi - lo + 1 != len(ref):
        raise ValueError(f"{gene.name} {name}: spans an intron")
    g_ref = ref if gene.strand == "+" else revcomp(ref)
    g_alt = alt if gene.strand == "+" else revcomp(alt)
    if len(g_ref) == 1 and len(g_alt) == 1:
        return lo, g_ref, g_alt
    anchor = fetch(fasta, gene.chrom, lo - 1, lo - 1)
    return lo - 1, anchor + g_ref, anchor + g_alt


def design(args):
    genes = {n: Gene(n, args.exons, args.reference) for n in ("BRCA1", "BRCA2")}
    rows = []
    for sample, seed in (("SPIKE_S1", 11), ("SPIKE_S2", 23)):
        rng = random.Random(seed)
        counter = 0
        for name, gene in genes.items():
            chrom, lo, hi = WINDOWS[name]
            seq = fetch(args.reference, chrom, lo, hi)
            base = lambda g: seq[g - lo]                       # noqa: E731
            taken = []
            events = []
            for kname, hgvs, c, ref, alt in KNOWN[sample]:
                if kname != name:
                    continue
                pos, g_ref, g_alt = known_event(gene, hgvs, c, ref, alt, args.reference)
                events.append((pos, g_ref, g_alt, "known", hgvs))
                taken.append(pos)
            sites = []
            coding_numbers = {n for n, _, _ in gene.coding_exons()}
            for idx, (number, start, end) in enumerate(gene.exons):
                if number not in coding_numbers:
                    continue
                # canonical splice sites of every second intron end, alternating donor/acceptor
                if idx % 2 == (0 if sample == "SPIKE_S1" else 1):
                    continue
                for g, which in ((start - rng.choice((1, 2)), "left"), (end + rng.choice((1, 2)), "right")):
                    first, last = idx == 0, idx == len(gene.exons) - 1
                    is_acceptor = (which == "left") == (gene.strand == "+")
                    if (first and is_acceptor) or (last and not is_acceptor):
                        continue
                    sites.append((number, g, g, "splice"))
            for number, start, end in gene.coding_exons():
                step = 260
                n = max(1, (end - start + 1) // step)
                for i in range(n):
                    a = start + i * (end - start + 1) // n
                    b = start + (i + 1) * (end - start + 1) // n - 1
                    sites.append((number, a, b, "coding"))
            for number, a, b, where in sites:
                kind = "SNV" if where == "splice" else SMALL_TYPES[counter % len(SMALL_TYPES)]
                size = int(kind[3:]) if kind[:3] in ("DEL", "INS", "DUP") else 1
                for _ in range(30):
                    pos = rng.randint(a, max(a, b - size)) if where == "coding" else a
                    if all(abs(pos - t) >= 150 for t in taken):
                        break
                else:
                    continue
                if where == "coding":
                    counter += 1
                if kind == "SNV":
                    ref = base(pos)
                    alt = rng.choice([x for x in "ACGT" if x != ref])
                elif kind.startswith("DEL"):
                    ref, alt = seq[pos - lo:pos - lo + size + 1], base(pos)
                elif kind.startswith("INS"):
                    ref = base(pos)
                    alt = ref + "".join(rng.choice("ACGT") for _ in range(size))
                else:                                           # tandem duplication of `size` bases
                    ref = base(pos)
                    alt = ref + seq[pos - lo + 1 - size:pos - lo + 1]
                events.append((pos, ref, alt, where, f"{name} exon {number}"))
                taken.append(pos)
            for pos, ref, alt, where, note in sorted(events):
                kind = "SNV" if len(ref) == len(alt) == 1 else ("DEL" if len(ref) > len(alt) else "INS")
                rows.append([sample, name, kind, chrom, pos, ref, alt, abs(len(ref) - len(alt)) or 1, where, note])
    with open(args.out, "w") as handle:
        handle.write("#sample\tgene\tkind\tchrom\tpos\tref\talt\tsize\tcontext\tnote\n")
        for r in rows:
            handle.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)} small spike-in variants written to {args.out}")


# ---------------------------------------------------------------- prepare
def in_window(chrom, lo, hi):
    for wc, ws, we in WINDOWS.values():
        if chrom == wc and lo >= ws and hi <= we:
            return True
    return False


def prepare(args):
    """Pass 1: fragments fully inside a window. Pass 2: background BAM without those pairs."""
    names, pairs = set(), {}
    for chrom, ws, we in WINDOWS.values():
        view = subprocess.Popen([args.samtools, "view", args.bam, f"{chrom}:{ws}-{we}"], stdout=subprocess.PIPE, text=True)
        for line in view.stdout:
            f = line.split("\t", 11)
            flag, pos, pnext, tlen = int(f[1]), int(f[3]), int(f[7]), int(f[8])
            if f[6] != "=" or tlen == 0 or abs(tlen) > 2000 or flag & 0xC:
                continue
            lo = min(pos, pnext)
            if not in_window(chrom, lo, lo + abs(tlen) - 1):
                continue
            names.add(f[0])
            if flag & 0x900:
                continue
            pairs.setdefault(f[0], []).append((flag, pos, tlen, len(f[9]), f[10]))
        view.wait()
    kept = 0
    with gzip.open(args.fragments, "wt") as out:
        for name, recs in pairs.items():
            if len(recs) != 2:
                continue
            fwd = [r for r in recs if not r[0] & 0x10]
            rev = [r for r in recs if r[0] & 0x10]
            if len(fwd) != 1 or len(rev) != 1 or fwd[0][0] & 0x400 or rev[0][0] & 0x400:
                continue
            fwd, rev = fwd[0], rev[0]
            flen = fwd[2]
            if not 50 <= flen <= 1000 or fwd[1] > rev[1]:
                continue
            out.write("\t".join(map(str, [name, fwd[1], flen, 1 if fwd[0] & 0x40 else 2, fwd[4], rev[4]])) + "\n")
            kept += 1
    # (the chromosome of a fragment follows from its position: the two windows do not share coordinates)
    view = subprocess.Popen([args.samtools, "view", "-h", args.bam], stdout=subprocess.PIPE, text=True)
    sink = subprocess.Popen([args.samtools, "view", "-@", "3", "-b", "-o", args.out, "-"], stdin=subprocess.PIPE, text=True)
    keep_contigs = {c for c, _, _ in WINDOWS.values()}
    rg = args.read_group
    tag = f"\tRG:Z:{rg}\n"
    dropped = other = 0
    header = True
    write = sink.stdin.write
    for line in view.stdout:
        if line[0] == "@":
            if line.startswith("@SQ"):
                if line.split("\t")[1][3:] in keep_contigs:
                    write(line)
            elif line.startswith("@HD"):
                write("@HD\tVN:1.6\tSO:coordinate\n")
            elif not line.startswith("@RG"):
                write(line)
            continue
        if header:
            write(f"@RG\tID:{rg}\tSM:{rg}\tPL:ILLUMINA\tLB:{rg}\n")
            header = False
        f = line.split("\t", 7)
        if f[0] in names:
            dropped += 1
            continue
        if f[6] != "=" and f[6] not in keep_contigs:
            other += 1
            continue
        if "\tSA:Z:" in line or "\tRG:Z:" in line:
            g = line.rstrip("\n").split("\t")
            line = "\t".join(g[:11] + [t for t in g[11:] if not t.startswith(("RG:Z:", "SA:Z:"))]) + "\n"
        write(line[:-1] + tag)
    sink.stdin.close()
    sink.wait()
    view.wait()
    subprocess.run([args.samtools, "index", args.out], check=True)
    print(json.dumps({"windowPairsReplaced": kept, "windowRecordsRemoved": dropped,
                      "recordsWithMateOnOtherContigRemoved": other}))


# ---------------------------------------------------------------- simulate
ERR = [10 ** (-q / 10) for q in range(100)]


def with_errors(seq, qual, rng):
    out = None
    for i, ch in enumerate(qual):
        if rng.random() < ERR[ord(ch) - 33]:
            if out is None:
                out = list(seq)
            out[i] = rng.choice([b for b in "ACGT" if b != seq[i]])
    return seq if out is None else "".join(out)


class Haplotype:
    """Reference window with edits (ref offset, deleted length, inserted sequence)."""

    def __init__(self, ref, start, edits):
        self.start = start
        self.edits = sorted(edits)
        parts, self.blocks, cursor, hap_len = [], [], 0, 0
        self.inserts = []                                  # (ref offset of insertion, hap offset, length)
        for off, dele, ins in self.edits:
            if off < cursor:
                raise ValueError("overlapping events")
            parts.append(ref[cursor:off])
            self.blocks.append((cursor, off, hap_len - cursor))   # ref [cursor, off) -> hap shift
            hap_len += off - cursor
            if ins:
                self.inserts.append((off, hap_len, len(ins)))
                parts.append(ins)
                hap_len += len(ins)
            cursor = off + dele
        parts.append(ref[cursor:])
        self.blocks.append((cursor, len(ref), hap_len - cursor))
        self.seq = "".join(parts)

    def lift(self, pos):
        off = pos - self.start
        for a, b, shift in self.blocks:
            if a <= off < b:
                return off + shift
        return None


def load_events(args, sample):
    """Return {gene: [(ref offset 0-based in window, del_len, ins_seq)]} plus truth records."""
    edits = {g: [] for g in WINDOWS}
    small, large = [], []
    for path in args.events:
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != sample:
                continue
            gene, kind, chrom = f[1], f[2], f[3]
            _, ws, _ = WINDOWS[gene]
            if kind in ("SNV", "DEL", "INS"):
                pos, ref, alt = int(f[4]), f[5], f[6]
                if fetch(args.reference, chrom, pos, pos + len(ref) - 1) != ref:
                    raise ValueError(f"reference mismatch at {chrom}:{pos}")
                edits[gene].append((pos - ws, len(ref), alt))
                small.append((chrom, pos, ref, alt))
            else:
                start, end = int(f[4]), int(f[5])
                record = {"sample": sample, "gene": gene, "kind": kind, "chrom": chrom, "start": start, "end": end,
                          "size": end - start + 1, "note": f[7] if len(f) > 7 else ""}
                if kind == "LDEL":
                    edits[gene].append((start - ws, end - start + 1, ""))
                elif kind == "LDUP":
                    edits[gene].append((end + 1 - ws, 0, fetch(args.reference, chrom, start, end)))
                elif kind == "ALU":
                    alu = fetch(args.reference, *ALU_SOURCE)
                    tsd = fetch(args.reference, chrom, start - TSD + 1, start)
                    edits[gene].append((start + 1 - ws, 0, alu + tsd))
                    record.update({"end": start, "size": len(alu) + TSD,
                                   "aluSource": "%s:%d-%d" % ALU_SOURCE, "targetSiteDuplication": TSD})
                else:
                    raise ValueError(kind)
                large.append(record)
    return edits, small, large


def simulate(args):
    edits, small, large = load_events(args, args.sample)
    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    stats = {"sample": args.sample, "pairs": 0, "hapA": 0, "hapB": 0, "skippedInDeletion": 0, "extraForInsertions": 0}
    haps = {}
    for gene, (chrom, ws, we) in WINDOWS.items():
        ref = fetch(args.reference, chrom, ws, we)
        haps[gene] = (chrom, ws, we, ref, Haplotype(ref, ws, edits[gene]))
    r1 = gzip.open(os.path.join(args.outdir, f"{args.sample}_R1.fastq.gz"), "wt", compresslevel=3)
    r2 = gzip.open(os.path.join(args.outdir, f"{args.sample}_R2.fastq.gz"), "wt", compresslevel=3)

    def emit(name, seq, hstart, flen, first, qf, qr):
        if hstart < 0 or hstart + flen > len(seq):
            return
        nf, nr = min(len(qf), flen), min(len(qr), flen)
        qf2, qr2 = qf[:nf], qr[len(qr) - nr:]
        fwd = with_errors(seq[hstart:hstart + nf], qf2, rng)
        rev = with_errors(seq[hstart + flen - nr:hstart + flen], qr2, rng)
        a = (fwd, qf2)
        b = (revcomp(rev), qr2[::-1])
        if first == 2:
            a, b = b, a
        r1.write(f"@{name}\n{a[0]}\n+\n{a[1]}\n")
        r2.write(f"@{name}\n{b[0]}\n+\n{b[1]}\n")
        stats["pairs"] += 1

    for line in gzip.open(args.fragments, "rt"):
        name, pos, flen, first, qf, qr = line.rstrip("\n").split("\t")
        pos, flen, first = int(pos), int(flen), int(first)
        gene = next(g for g, (c, ws, we) in WINDOWS.items() if ws <= pos <= we)
        chrom, ws, we, ref, hap = haps[gene]
        if stable_random(name, args.seed) < 0.5:
            emit(name, ref, pos - ws, flen, first, qf, qr)
            stats["hapA"] += 1
            # positions of haplotype-A fragments are reused for sequence that only exists on B
            for off, hoff, length in hap.inserts:
                if off - length <= pos - ws < off:
                    emit(name + ":ins", hap.seq, hoff + (pos - ws) - (off - length), flen, first, qf, qr)
                    stats["extraForInsertions"] += 1
        else:
            h = hap.lift(pos)
            if h is None:
                stats["skippedInDeletion"] += 1
                continue
            emit(name, hap.seq, h, flen, first, qf, qr)
            stats["hapB"] += 1
    r1.close()
    r2.close()
    with open(os.path.join(args.outdir, f"{args.sample}.truth.vcf"), "w") as handle:
        handle.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        for chrom, length in (("chr13", 114364328), ("chr17", 83257441)):
            handle.write(f"##contig=<ID={chrom},length={length}>\n")
        handle.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n")
        for chrom, pos, ref, alt in sorted(small):
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0/1\n")
    json.dump({"stats": stats, "largeEvents": large, "smallVariants": len(small)},
              open(os.path.join(args.outdir, f"{args.sample}.truth.json"), "w"), indent=2)
    print(json.dumps(stats))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("design")
    d.add_argument("--exons", required=True)
    d.add_argument("--reference", required=True)
    d.add_argument("--out", required=True)
    q = sub.add_parser("prepare")
    q.add_argument("--bam", required=True)
    q.add_argument("--fragments", required=True)
    q.add_argument("--out", required=True)
    q.add_argument("--read-group", default="spike")
    q.add_argument("--samtools", default="samtools")
    s = sub.add_parser("simulate")
    s.add_argument("--sample", required=True)
    s.add_argument("--events", nargs="+", required=True)
    s.add_argument("--fragments", required=True)
    s.add_argument("--reference", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--seed", type=int, default=1)
    args = p.parse_args(argv)
    return {"design": design, "prepare": prepare, "simulate": simulate}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())

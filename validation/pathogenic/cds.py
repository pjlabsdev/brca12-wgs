"""MANE Select coding sequence model of BRCA1 (NM_007294.4) and BRCA2 (NM_000059.4) on GRCh38.

Built from assets/brca12_exons.bed and the reference FASTA. The position of the start codon is
the only constant; it is checked on load: the CDS must start with ATG, contain no internal stop
codon, end with a stop codon and have the known protein length.
"""
import subprocess

ATG = {"BRCA1": ("chr17", 43124096, "-"), "BRCA2": ("chr13", 32316461, "+")}   # genomic position of c.1
PROTEIN_LENGTH = {"BRCA1": 1863, "BRCA2": 3418}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODON = {}
_bases = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_bases):
    for _j, _b in enumerate(_bases):
        for _k, _c in enumerate(_bases):
            CODON[_a + _b + _c] = _aa[16 * _i + 4 * _j + _k]


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def fetch(fasta, chrom, start, end, samtools="samtools"):
    """1-based inclusive."""
    out = subprocess.run([samtools, "faidx", fasta, f"{chrom}:{start}-{end}"], capture_output=True, text=True, check=True)
    return "".join(out.stdout.splitlines()[1:]).upper()


class Gene:
    def __init__(self, name, exon_bed, fasta, samtools="samtools"):
        self.name = name
        self.chrom, atg, self.strand = ATG[name]
        self.exons = []                      # (number, start1, end1) in transcript order
        with open(exon_bed) as handle:
            for line in handle:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4 and f[3].startswith(name + "_ex"):
                    self.exons.append((int(f[3].split("_ex")[1]), int(f[1]) + 1, int(f[2])))
        self.exons.sort()
        self.cds_pos = []                    # genomic position of every CDS base, c.1 first
        started = False
        for _, start, end in self.exons:
            rng = range(start, end + 1) if self.strand == "+" else range(end, start - 1, -1)
            for g in rng:
                started = started or g == atg
                if started:
                    self.cds_pos.append(g)
        lo, hi = min(self.cds_pos), max(self.cds_pos)
        seq = fetch(fasta, self.chrom, lo, hi, samtools)
        base = {g: seq[g - lo] for g in self.cds_pos}
        tx = "".join(base[g] if self.strand == "+" else base[g].translate(COMP) for g in self.cds_pos)
        protein = "".join(CODON.get(tx[i:i + 3], "X") for i in range(0, len(tx) - 2, 3))
        stop = protein.find("*")
        if not tx.startswith("ATG") or stop != PROTEIN_LENGTH[name]:
            raise ValueError(f"{name}: CDS check failed (start {tx[:3]}, first stop after {stop} codons)")
        self.cds_pos = self.cds_pos[:3 * (stop + 1)]
        self.cds = tx[:3 * (stop + 1)]
        self.protein = protein[:stop]
        self.index = {g: i + 1 for i, g in enumerate(self.cds_pos)}

    def genomic(self, c):
        return self.cds_pos[c - 1]

    def coding(self, g):
        return self.index.get(g)

    def exon_of(self, g):
        for number, start, end in self.exons:
            if start <= g <= end:
                return number
        return None

    def coding_exons(self):
        """(number, genomic start, genomic end) of the coding part of each exon."""
        out = []
        for number, start, end in self.exons:
            inside = [g for g in (start, end) if g in self.index]
            pos = [g for g in range(start, end + 1) if g in self.index]
            if pos:
                out.append((number, pos[0], pos[-1]))
        return out

"""GRCh37 coordinates and BIC-style legacy names for BRCA1/BRCA2 variants.

GRCh37: both analysed loci map to GRCh37 as one contiguous block (Ensembl map endpoint,
GRCh38 -> GRCh37): chr17:43044295-43170327 -> 41196312-41322344 (offset -1,847,983) and
chr13:32315077-32400268 -> 32889214-32974405 (offset +574,137). Positions outside these loci are
not converted.

BIC names follow the Breast Cancer Information Core nucleotide numbering used in older reports:
BRCA1 nt = c. + 119 (U14680), BRCA2 nt = c. + 228 (U43746). Deletions and duplications are
shifted to their 5'-most position in the coding sequence, as BIC did (HGVS uses the 3'-most one):
c.68_69del -> 185delAG, c.5266dup -> 5382insC, c.5946del -> 6174delT. Intronic variants use the
legacy intron number (IVS); BRCA1 legacy exon numbering skips exon 4. The names are computed and
meant for comparison with older reports only; HGVS stays the reference nomenclature.
"""
import csv
import os
import re

ASSETS = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "assets")

GRCH37 = {
    "chr17": (43044295, 43170327, -1847983),
    "chr13": (32315077, 32400268, 574137),
}
BIC_OFFSET = {"BRCA1": 119, "BRCA2": 228}
SPELL_MAX = 4          # deleted/inserted bases up to this length are written out (185delAG, 6486delCAGA)

_CDS = {}
_EXONS = {}


def grch37(chrom, pos):
    """GRCh37 position for a GRCh38 position inside the BRCA1/BRCA2 loci, else None."""
    block = GRCH37.get(chrom if chrom.startswith("chr") else "chr" + chrom)
    if not block or pos is None or not block[0] <= pos <= block[1]:
        return None
    return pos + block[2]


def _cds(gene):
    if not _CDS:
        name, seq = None, []
        with open(os.path.join(ASSETS, "legacy", "mane_cds.fa")) as handle:
            for line in handle:
                line = line.strip()
                if line.startswith(">"):
                    if name:
                        _CDS[name] = "".join(seq)
                    name, seq = line[1:].split()[0], []
                elif line:
                    seq.append(line)
        if name:
            _CDS[name] = "".join(seq)
    return _CDS.get(gene)


def _exons(gene):
    """MANE exons as (c_start, c_end, legacy_number) from the ENIGMA exon table."""
    if gene not in _EXONS:
        rows = []
        path = os.path.join(ASSETS, "enigma", "exon_weights.tsv")
        if os.path.exists(path):
            with open(path, newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if row["gene"] != gene:
                        continue
                    match = re.match(r"(\d+)(?: \((\d+)\))?", row["exon"])
                    if not match:
                        continue
                    legacy = int(match.group(2) or match.group(1))
                    rows.append((_num(row["c_start"]), _num(row["c_end"]), legacy))
        _EXONS[gene] = rows
    return _EXONS[gene]


def _num(text):
    try:
        return int(str(text).replace("*", ""))
    except ValueError:
        return None


def _left_shift(seq, start, end):
    """Shift a 1-based CDS interval [start, end] to its 5'-most equivalent position."""
    while start > 1 and seq[start - 2] == seq[end - 1]:
        start -= 1
        end -= 1
    return start, end


def _nt(pos, off):
    """BIC nucleotide number; c. has no position 0, so 5'UTR positions are shifted by one more."""
    return pos + off + 1 if pos < 0 else pos + off


def _bases(text):
    return text if len(text) <= SPELL_MAX else str(len(text))


def _intron_name(gene, anchor, sign, offset, change):
    """IVS number = legacy number of the exon before the intron (BRCA1 has no legacy exon 4)."""
    exons = sorted((e for e in _exons(gene) if e[0] is not None and e[1] is not None), key=lambda e: e[0])
    for i, (c_start, c_end, legacy) in enumerate(exons):
        if sign == "+" and anchor == c_end:
            return f"IVS{legacy}+{offset}{change}"
        if sign == "-" and anchor == c_start and i > 0:
            return f"IVS{exons[i - 1][2]}-{offset}{change}"
    return None


def bic_names(gene, hgvsc):
    """BIC-style legacy names for a c. notation of the MANE transcript (list, may be empty).

    BIC itself is not consistent for deletions and insertions: most entries use the HGVS (3'-most)
    position, some well-known ones the 5'-most position (185delAG, 5382insC). Both are returned,
    the HGVS-based name first.
    """
    if gene not in BIC_OFFSET or not hgvsc:
        return []
    c = hgvsc.split(":", 1)[1] if ":" in hgvsc else hgvsc
    off = BIC_OFFSET[gene]
    match = re.match(r"c\.(-?\d+)([+-])(\d+)([ACGT]>[ACGT])$", c)
    if match:
        name = _intron_name(gene, int(match.group(1)), match.group(2), int(match.group(3)), match.group(4))
        return [name] if name else []
    match = re.match(r"c\.(-?\d+)([ACGT])>([ACGT])$", c)
    if match:
        return [f"{_nt(int(match.group(1)), off)}{match.group(2)}>{match.group(3)}"]
    seq = _cds(gene)
    match = re.match(r"c\.(\d+)(?:_(\d+))?(del|dup)$", c)
    if match and seq:
        start = int(match.group(1))
        end = int(match.group(2) or start)
        if not 1 <= start <= end <= len(seq):
            return []
        names = []
        for s0, e0 in ((start, end), _left_shift(seq, start, end)):
            bases = seq[s0 - 1:e0]
            if match.group(3) == "del":
                names.append(f"{s0 + off}del{_bases(bases)}")
            else:
                names.append(f"{(e0 if (s0, e0) == (start, end) else s0 - 1) + off}ins{_bases(bases)}")
        return list(dict.fromkeys(names))
    match = re.match(r"c\.(\d+)_(\d+)ins([ACGT]+)$", c)
    if match:
        return [f"{int(match.group(1)) + off}ins{_bases(match.group(3))}"]
    return []


def bic_name(gene, hgvsc):
    names = bic_names(gene, hgvsc)
    return names[0] if names else None

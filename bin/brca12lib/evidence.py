"""Optional evidence sources for the ENIGMA hints: gnomAD non-cancer, SpliceAI and BayesDel.

All files are supplied by the user (licences: gnomAD open, SpliceAI and dbNSFP/BayesDel free for
non-commercial research only). Each file is either bgzip-compressed with a tabix index (queried with
`tabix`) or plain text (scanned, intended for tests and small extracts).

gnomAD   VCF (v2.1 exomes, v3.1/v4 genomes). Filtering allele frequencies are read from INFO keys
         that contain "faf95"; keys of the non-cancer subset ("non_cancer") are preferred.
         Founder or mixed populations (asj, fin, ami, oth, remaining, mid) are skipped.
SpliceAI VCF with INFO SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL.
BayesDel TSV: chrom, pos, ref, alt, BayesDel_noAF_score (e.g. extracted from dbNSFP).
"""
import gzip
import re
import subprocess

SKIP_POPULATIONS = {"asj", "fin", "ami", "oth", "remaining", "mid"}
_NOT_A_POPULATION = {"faf95", "non", "cancer", "popmax", "max", "adj", "grpmax"}
_SEX = {"xx", "xy", "male", "female"}


class TabixFile:
    def __init__(self, path, tabix="tabix"):
        self.path = path
        self.tabix = tabix
        self.indexed = path.endswith(".gz")
        self._lines = None

    def _all(self):
        if self._lines is None:
            opener = gzip.open if self.path.endswith(".gz") else open
            with opener(self.path, "rt") as handle:
                self._lines = [line.rstrip("\n").split("\t") for line in handle if not line.startswith("#")]
        return self._lines

    def query(self, chrom, pos):
        if self.indexed:
            names = [chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom]
            for name in names:
                out = subprocess.run([self.tabix, self.path, f"{name}:{pos}-{pos}"],
                                     capture_output=True, text=True)
                if out.returncode == 0 and out.stdout.strip():
                    return [line.split("\t") for line in out.stdout.splitlines()]
            return []
        bare = chrom[3:] if chrom.startswith("chr") else chrom
        return [f for f in self._all() if len(f) > 1 and f[0] in (chrom, bare, "chr" + bare) and f[1] == str(pos)]


def _info(field):
    out = {}
    for item in field.split(";"):
        key, _, value = item.partition("=")
        out[key] = value if value else True
    return out


def _float(value):
    try:
        return float(str(value).split(",")[0])
    except (TypeError, ValueError):
        return None


def _match(rows, ref, alt, ref_col=3, alt_col=4):
    for f in rows:
        if len(f) > alt_col and f[ref_col] == ref and alt in f[alt_col].split(","):
            return f
    return None


def gnomad_lookup(files, chrom, pos, ref, alt, tabix="tabix"):
    """Highest non-founder FAF per dataset. Returns None when no gnomAD file is configured."""
    if not files:
        return None
    result = {"datasets": [], "present": False, "best": None}
    for label, path in files:
        result["datasets"].append(label)
        hit = _match(TabixFile(path, tabix).query(chrom, pos), ref, alt)
        if not hit:
            continue
        info = _info(hit[7]) if len(hit) > 7 else {}
        ac_keys = [k for k in info if re.fullmatch(r"(non_cancer_)?AC|AC_non_cancer", k)]
        if any((_float(info[k]) or 0) > 0 for k in ac_keys) or not ac_keys:
            result["present"] = True
        candidates = []
        for key, value in info.items():
            tokens = key.lower().replace("-", "_").split("_")
            if "faf95" not in tokens or "joint" in tokens or _SEX & set(tokens):
                continue
            populations = [t for t in tokens if t not in _NOT_A_POPULATION]
            population = populations[-1] if populations else "max"
            if population in SKIP_POPULATIONS:
                continue
            faf = _float(value)
            if faf is None:
                continue
            candidates.append(("cancer" in tokens, faf, population))
        if not candidates:
            continue
        non_cancer = [c for c in candidates if c[0]]
        pool = non_cancer or candidates
        nc, faf, population = max(pool, key=lambda c: c[1])
        if result["best"] is None or faf > result["best"]["faf"]:
            result["best"] = {"dataset": label, "population": population, "faf": faf, "nonCancer": nc}
    return result


def spliceai_lookup(path, chrom, pos, ref, alt, gene, tabix="tabix"):
    if not path:
        return None
    hit = _match(TabixFile(path, tabix).query(chrom, pos), ref, alt)
    if not hit or len(hit) < 8:
        return None
    scores = []
    for entry in str(_info(hit[7]).get("SpliceAI", "")).split(","):
        parts = entry.split("|")
        if len(parts) >= 6 and parts[0] == alt and (not gene or parts[1] == gene):
            scores += [s for s in (_float(p) for p in parts[2:6]) if s is not None]
    return max(scores) if scores else None


def bayesdel_lookup(path, chrom, pos, ref, alt, tabix="tabix"):
    if not path:
        return None
    hit = _match(TabixFile(path, tabix).query(chrom, pos), ref, alt, ref_col=2, alt_col=3)
    return _float(hit[4]) if hit and len(hit) > 4 else None


class Evidence:
    def __init__(self, gnomad_exomes=None, gnomad_genomes=None, spliceai=None, bayesdel=None, tabix="tabix"):
        self.gnomad = [(label, p) for label, p in (("gnomAD exomes", gnomad_exomes),
                                                    ("gnomAD genomes", gnomad_genomes)) if p]
        self.spliceai = spliceai
        self.bayesdel = bayesdel
        self.tabix = tabix

    def configured(self):
        return {"gnomad": [label for label, _ in self.gnomad], "spliceai": bool(self.spliceai),
                "bayesdel": bool(self.bayesdel)}

    def for_variant(self, chrom, pos, ref, alt, gene):
        return {"gnomad": gnomad_lookup(self.gnomad, chrom, pos, ref, alt, self.tabix),
                "spliceai": spliceai_lookup(self.spliceai, chrom, pos, ref, alt, gene, self.tabix),
                "bayesdel": bayesdel_lookup(self.bayesdel, chrom, pos, ref, alt, self.tabix)}

#!/usr/bin/env python3
"""Evaluate the frozen Delly-versus-Manta test (PREREGISTRATION.md).
Usage: eval.py OUTDIR  -> results/sv_test.tsv, results/sv_real.tsv and results/sv_summary.txt"""
import glob
import gzip
import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pathogenic", "founder"))
from eval import clopper_pearson  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
VAL = os.path.dirname(HERE)
ROOT = os.path.dirname(VAL)
TYPES = {"LDEL": "DEL", "LDUP": "DUP", "LINV": "INV", "ALU": "INS"}


def read_events(path, prefix=None):
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if prefix and not f[0].startswith(prefix):
            continue
        out.setdefault(f[0], []).append((TYPES[f[2]], f[3], int(f[4]), int(f[5])))
    return out


NEW = read_events(os.path.join(HERE, "events.tsv"))
TDUP = read_events(os.path.join(VAL, "tandem_dup", "events.tsv"))
SPIKE = read_events(os.path.join(VAL, "pathogenic", "large_events.tsv"))
GENES = [(l.split("\t")[0], int(l.split("\t")[1]), int(l.split("\t")[2])) for l in open(os.path.join(ROOT, "assets", "brca12_regions.bed"))
         if not l.startswith("#")]


def tol(size):
    return max(20, min(1000, size // 2))


def in_gene(chrom, a, b):
    return any(c == chrom and a <= e and b > s for c, s, e in GENES)


def vcf_records(path):
    """(type, chrom, a, b) for PASS records of a Manta or Delly VCF; BND as a pair of positions."""
    out = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[6] != "PASS":
                continue
            info = dict(x.split("=", 1) if "=" in x else (x, True) for x in f[7].split(";"))
            t = info.get("SVTYPE", "")
            pos = int(f[1])
            if t == "BND":
                m = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", f[4])
                if not m or m.group(1) != f[0]:
                    out.append(("BND", f[0], pos, pos))
                    continue
                a, b = sorted((pos, int(m.group(2))))
                out.append(("BND", f[0], a, b))
            elif t in ("DEL", "DUP", "INV", "INS"):
                out.append((t, f[0], pos, int(info.get("END", pos))))
    return out


def locus(records):
    return [r for r in records if in_gene(r[1], r[2], r[3])]


WINDOWS = [(c, s - 1000000, e + 1000000) for c, s, e in GENES]


def window_filter(records):
    """Pipeline filter after Delly (bin/brca12_filter_sv.py, ADDENDUM_filter.md)."""
    def inside(c, p):
        return any(c == wc and ws < p <= we for wc, ws, we in WINDOWS)
    return [r for r in records if inside(r[1], r[2]) and
            (r[0] == "BND" or (inside(r[1], r[3]) and r[3] - r[2] >= 50))]


def ours(path):
    out = []
    for r in json.load(open(path)):
        sig = r.get("signature", [])
        t = "DEL" if "read-pair deletion" in sig else "INV" if "inversion" in sig else \
            "DUP" if "tandem duplication" in sig else None
        if t:
            out.append((t, r["chrom"], r["pos"], r["end"]))
    return [r for r in out if in_gene(r[1], r[2], r[3])]


def matches(call, event):
    t, c, a, b = call
    et, ec, es, ee = event
    if c != ec or et == "INS":
        return et == "INS" and c == ec and abs(a - es) <= 1000
    ok_type = t == et or (et == "INV" and t == "BND") or (et == "DUP" and t == "BND")
    return ok_type and abs(a - es) <= tol(ee - es + 1) + 1 and abs(b - ee) <= tol(ee - es + 1) + 1


def score(calls, events):
    found = [any(matches(c, e) for c in calls) for e in events]
    false = [c for c in calls if not any(matches(c, e) for e in events)]
    return found, false


def main():
    out = sys.argv[1]
    rows = []
    for path in sorted(glob.glob(os.path.join(out, "res", "*.json"))):
        t = os.path.basename(path)[:-5]
        name, f = t.rsplit(".", 1)
        events = NEW.get(name) or TDUP.get(name) or SPIKE.get(name) or []
        group = "new test" if name in NEW else "tandem duplications (../tandem_dup)" if name in TDUP else \
            "earlier spike-ins (seen in tuning)"
        raw = vcf_records(os.path.join(out, "res", t + ".delly.vcf"))
        arms = {"delly": locus(window_filter(raw)), "dellyUnfiltered": locus(raw),
                "manta": locus(vcf_records(os.path.join(out, "res", t + ".manta.vcf.gz"))),
                "readpairs": ours(path)}
        res = {k: score(v, events) for k, v in arms.items()}
        for i, e in enumerate(events or [None]):
            row = {"test": t, "sample": name, "group": group, "reads": {"f1": "all", "f2": "1/2", "f3": "1/3"}[f],
                   "type": e[0] if e else "none", "chrom": e[1] if e else "", "start": e[2] if e else 0,
                   "end": e[3] if e else 0, "size": e[3] - e[2] + 1 if e else 0}
            for k in arms:
                row[k] = ("yes" if res[k][0][i] else "no") if e else "n/a"
                row[k + "False"] = len(res[k][1]) if i == 0 else 0
            rows.append(row)
    cols = list(rows[0])
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "sv_test.tsv"), "w") as h:
        h.write("\t".join(cols) + "\n")
        for r in rows:
            h.write("\t".join(str(r[c]) for c in cols) + "\n")

    # real genomes: calls versus ONT long-read events of 50 bp or more (validation/longread)
    sites = {}
    for l in open(os.path.join(VAL, "longread", "results", "all_disc.sites.tsv")):
        if not l.startswith("#"):
            f = l.split("\t")
            sites[f[0]] = (f[5], f[2], int(f[3]), int(f[4]) if f[5] == "DEL" else int(f[3]) + int(f[6]))
    ont = {}
    for l in open(os.path.join(VAL, "longread", "results", "all_disc.carriers.tsv")):
        if not l.startswith("#"):
            f = l.split("\t")
            if f[1] in sites:
                ont.setdefault(f[0], []).append(sites[f[1]])
    extra = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2] for l in open(os.path.join(HERE, "real_carriers.tsv"))
             if not l.startswith("#")}
    real = []
    for path in sorted(glob.glob(os.path.join(out, "real", "*.json"))):
        s = os.path.basename(path)[:-5]
        truth = [e for e in ont.get(s, []) if e[0] in ("DEL", "INV")]
        raw = vcf_records(os.path.join(out, "real", s + ".delly.vcf"))
        arms = {"delly": locus(window_filter(raw)), "dellyUnfiltered": locus(raw), "readpairs": ours(path)}
        row = {"sample": s, "set": "real carrier list" if s in extra else "200 frozen genomes",
               "ontDeletions": ";".join(f"{c}:{a}-{b}" for _, c, a, b in truth) or "none",
               "inOntCohort": "yes" if s in ont or s not in extra else "unknown"}
        for k, calls in arms.items():
            found, false = score([c for c in calls if c[0] in ("DEL", "INV", "BND")], truth)
            row[k + "Found"] = sum(found)
            row[k + "Unconfirmed"] = ";".join(f"{t}:{c}:{a}-{b}" for t, c, a, b in false)
        real.append(row)
    if real:
        with open(os.path.join(HERE, "results", "sv_real.tsv"), "w") as h:
            h.write("\t".join(real[0]) + "\n")
            for r in real:
                h.write("\t".join(str(v) for v in r.values()) + "\n")

    def rate(sub, key):
        k = sum(1 for r in sub if r[key] == "yes")
        lo, hi = clopper_pearson(k, len(sub))
        return f"{k}/{len(sub)} ({100 * lo:.1f}-{100 * hi:.1f} %)"
    lines = ["Delly 1.3.3 versus Manta 1.6, frozen test. Simulated sets are SIMULATED reads (real HG002 positions,",
             "fragment sizes and base qualities); the 1000 Genomes samples are real reads. Read-pair = research arm.", ""]
    for group in ("new test", "tandem duplications (../tandem_dup)", "earlier spike-ins (seen in tuning)"):
        g = [r for r in rows if r["group"] == group and r["type"] != "none"]
        if not g:
            continue
        lines.append(group)
        for typ in ("DEL", "INV", "DUP", "INS"):
            sub = [r for r in g if r["type"] == typ]
            if not sub:
                continue
            lines.append(f"  {typ}: Delly {rate(sub, 'delly')} (unfiltered {rate(sub, 'dellyUnfiltered')}), Manta {rate(sub, 'manta')}, read pairs {rate(sub, 'readpairs')}")
            for reads in ("all", "1/2", "1/3"):
                s2 = [r for r in sub if r["reads"] == reads]
                lines.append(f"    reads {reads}: Delly {rate(s2, 'delly')}, Manta {rate(s2, 'manta')}, read pairs {rate(s2, 'readpairs')}")
            if typ == "DEL":
                for lo_, hi_ in ((50, 300), (300, 1000), (1000, 20001)):
                    s2 = [r for r in sub if lo_ <= r["size"] < hi_]
                    lines.append(f"    size {lo_}-{hi_ - 1} bp: Delly {rate(s2, 'delly')}, Manta {rate(s2, 'manta')}, read pairs {rate(s2, 'readpairs')}")
    fsim = [r for r in rows if r["group"] != "earlier spike-ins (seen in tuning)" or r["sample"] in ("SPIKE_NEG", "SPIKE_S1", "SPIKE_S2")]
    n_an = len({r["test"] for r in fsim})
    lines.append("")
    lines.append(f"false calls in {n_an} simulated analyses (new test, tandem duplications, SPIKE_NEG/S1/S2): "
                 f"Delly {sum(r['dellyFalse'] for r in fsim)} (unfiltered {sum(r['dellyUnfilteredFalse'] for r in fsim)}), Manta {sum(r['mantaFalse'] for r in fsim)}, "
                 f"read pairs {sum(r['readpairsFalse'] for r in fsim)}")
    if real:
        for st in ("200 frozen genomes", "real carrier list"):
            sub = [r for r in real if r["set"] == st]
            n_true = sum(len(r["ontDeletions"].split(";")) for r in sub if r["ontDeletions"] != "none")
            lines.append(f"{st} ({len(sub)}): ONT deletions {n_true}; found Delly {sum(r['dellyFound'] for r in sub)}, "
                         f"read pairs {sum(r['readpairsFound'] for r in sub)}; samples with unconfirmed calls Delly "
                         f"{sum(1 for r in sub if r['dellyUnconfirmed'])} (unfiltered {sum(1 for r in sub if r['dellyUnfilteredUnconfirmed'])}), read pairs {sum(1 for r in sub if r['readpairsUnconfirmed'])}")
    text = "\n".join(lines) + "\n"
    open(os.path.join(HERE, "results", "sv_summary.txt"), "w").write(text)
    print(text)


if __name__ == "__main__":
    main()

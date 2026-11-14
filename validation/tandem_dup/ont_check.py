#!/usr/bin/env python3
"""Check a short-read tandem duplication call against the sample's ONT long reads (1000 Genomes ONT, Vienna).
A carrier read either holds an insertion of about the duplicated length or is split, with a supplementary
alignment on the same chromosome that starts near the duplication start while the primary ends near its end
(or the reverse). Reads that span the whole segment plus 500 bp on each side without either count as reference.
Usage: ont_check.py SAMPLE CHROM START END REFERENCE   (needs samtools)"""
import re
import subprocess
import sys

URL = "https://1kg-ont-vienna.s3.amazonaws.com/hg38/{}.hg38.cram"
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def main():
    sample, chrom, start, end, ref = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    size, tol = end - start, 1000
    sam = subprocess.run(["samtools", "view", "-T", ref, URL.format(sample), f"{chrom}:{start - 5000}-{end + 5000}"],
                         capture_output=True, text=True, check=True).stdout
    alt, refreads = set(), set()
    for line in sam.splitlines():
        f = line.split("\t")
        flag, pos = int(f[1]), int(f[3])
        if flag & 0x904 or int(f[4]) < 20:
            continue
        cig = CIGAR.findall(f[5])
        span = sum(int(n) for n, op in cig if op in "MD=X")
        ins = [int(n) for n, op in cig if op == "I" and 0.5 * size <= int(n) <= 1.5 * size]
        sa = [e.split(",") for x in f[11:] if x.startswith("SA:Z:") for e in x[5:].rstrip(";").split(";")]
        split = any(s[0] == chrom and (abs(int(s[1]) - start) <= tol and abs(pos + span - end) <= tol or
                                       abs(pos - start) <= tol and int(s[1]) > start + 0.5 * size) for s in sa)
        if ins or split:
            alt.add(f[0])
        elif pos <= start - 500 and pos + span >= end + 500:
            refreads.add(f[0])
    verdict = "confirmed" if len(alt) >= 2 else "not confirmed"
    print("\t".join(map(str, [sample, f"{chrom}:{start}-{end}", len(alt), len(refreads), verdict])))


if __name__ == "__main__":
    main()

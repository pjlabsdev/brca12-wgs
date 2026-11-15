#!/usr/bin/env python3
"""Hard-clip soft-clipped poly-G read-through tails (SAM stdin -> SAM stdout).

Two-colour Illumina chemistry reads 'no signal' as G, so a read that runs past its insert ends in
a G run. Aligners soft-clip it, but HaplotypeCaller still assembles soft-clipped bases and can turn
such tails into false SNVs/indels. A soft clip at the sequencing 3' end (right end of a forward
read, left end of a reverse read, where the run shows as C) whose first WINDOW bases next to the
alignment are at least MIN_FRAC G (C) is turned into a hard clip. Nothing else is changed.
"""
import re, sys
WINDOW, MIN_FRAC = 10, 0.8
CIG = re.compile(r"(\d+)([MIDNSHP=X])")

def main():
    n = 0
    for line in sys.stdin:
        if line.startswith("@"):
            sys.stdout.write(line); continue
        f = line.rstrip("\n").split("\t")
        if f[5] == "*" or f[9] == "*":
            sys.stdout.write(line); continue
        ops = [(int(a), b) for a, b in CIG.findall(f[5])]
        rev = int(f[1]) & 16
        seq, qual = f[9], f[10]
        if not rev and ops[-1][1] == "S":
            k = ops[-1][0]; tail = seq[-k:][:WINDOW]
            if len(tail) >= min(WINDOW, k) and k >= 5 and tail.count("G") >= MIN_FRAC * len(tail):
                seq, qual = seq[:-k], (qual[:-k] if qual != "*" else qual)
                ops[-1] = (k, "H"); n += 1
        elif rev and ops[0][1] == "S":
            k = ops[0][0]; head = seq[:k][-WINDOW:]
            if k >= 5 and head.count("C") >= MIN_FRAC * len(head):
                seq, qual = seq[k:], (qual[k:] if qual != "*" else qual)
                ops[0] = (k, "H"); n += 1
        else:
            sys.stdout.write(line); continue
        f[5] = "".join(f"{a}{b}" for a, b in ops); f[9] = seq; f[10] = qual
        sys.stdout.write("\t".join(f) + "\n")
    sys.stderr.write(f"poly-G tails hard-clipped: {n}\n")

if __name__ == "__main__":
    main()

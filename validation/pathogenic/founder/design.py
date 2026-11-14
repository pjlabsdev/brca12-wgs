#!/usr/bin/env python3
"""Design of the simulated carriers for the BRCA2 c.156_157insAlu test (deterministic, seeded).

The real founder element has not been published base by base (Teugels et al. 2005 give subfamily
AluYa5, orientation and the 14 bp target site duplication). The catalog therefore carries the Dfam
AluYa5 consensus. To test that a carrier is still found when the real element differs from the
consensus, every simulated carrier gets its own element: the AluYa5 consensus body with random
substitutions (80 %) and 1 bp insertions or deletions (20 %) at a set divergence, and a poly-A tail
of a set length. The element is written to elements/<sample>.fa and the event to events.tsv
(spike.py format; the target site duplication is added by spike.py).

Usage: design.py ALUYA5_CONSENSUS_FILE OUTDIR
"""
import os
import random
import sys

DIVERGENCE = (0.0, 0.02, 0.05, 0.08, 0.10)
POLYA = (15, 30)
REPLICATES = 3
SEED = 20261115


def mutate(body, divergence, rng):
    out = []
    for base in body:
        r = rng.random()
        if r < divergence * 0.8:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        elif r < divergence * 0.9:
            continue                                  # 1 bp deletion
        elif r < divergence:
            out.append(base + rng.choice("ACGT"))     # 1 bp insertion
        else:
            out.append(base)
    return "".join(out)


def main():
    consensus, outdir = sys.argv[1], sys.argv[2]
    with open(consensus) as handle:
        seq = "".join(l.strip() for l in handle if not l.startswith(">")).upper()
    body = seq.rstrip("A")
    rng = random.Random(SEED)
    os.makedirs(os.path.join(outdir, "elements"), exist_ok=True)
    with open(os.path.join(outdir, "events.tsv"), "w") as events:
        events.write("#sample\tgene\tkind\tchrom\tstart\tend\tsize\tnote\n")
        for d in DIVERGENCE:
            for a in POLYA:
                for r in range(1, REPLICATES + 1):
                    name = f"FOUNDER_d{int(d * 100):02d}_a{a}_r{r}"
                    element = mutate(body, d, rng) + "A" * a
                    with open(os.path.join(outdir, "elements", name + ".fa"), "w") as fa:
                        fa.write(f">{name} AluYa5 consensus, divergence {d:.2f}, poly-A {a}\n{element}\n")
                    events.write(f"{name}\tBRCA2\tALU\tchr13\t32319165\t32319165\t{len(element) + 14}\t"
                                 f"simulated c.156_157insAlu, element {d:.0%} from AluYa5 consensus, poly-A {a}\n")


if __name__ == "__main__":
    main()

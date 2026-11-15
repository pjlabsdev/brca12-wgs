"""Poly-G read-through tails in soft clips are hard-clipped before small-variant calling."""
import os
import subprocess
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(ROOT, "bin", "brca12_clip_polyg.py")


def run(records):
    sam = "@HD\tVN:1.6\n" + "".join("\t".join(r) + "\n" for r in records)
    out = subprocess.run([sys.executable, TOOL], input=sam, capture_output=True, text=True, check=True)
    return [l.split("\t") for l in out.stdout.splitlines() if not l.startswith("@")], out.stderr


def rec(flag, cigar, seq):
    return ["r", str(flag), "chr17", "100", "60", cigar, "=", "300", "200", seq, "I" * len(seq), "NM:i:0"]


class ClipPolyG(unittest.TestCase):
    def test_forward_polyg_tail_becomes_hard_clip(self):
        seq = "ACGT" * 5 + "GGGGGGGGCTGGGGG"
        (r,), err = run([rec(99, "20M15S", seq)])
        self.assertEqual(r[5], "20M15H")
        self.assertEqual(r[9], "ACGT" * 5)
        self.assertEqual(len(r[10]), 20)
        self.assertIn("hard-clipped: 1", err)

    def test_reverse_read_shows_the_run_as_leading_c(self):
        seq = "TTCCCCCCCCCC" + "ACGT" * 5
        (r,), _ = run([rec(83, "12S20M", seq)])
        self.assertEqual(r[5], "12H20M")
        self.assertEqual(r[9], "ACGT" * 5)
        self.assertEqual(r[3], "100")

    def test_other_soft_clips_and_plain_reads_are_untouched(self):
        records = [
            rec(99, "20M15S", "ACGT" * 5 + "ACGTTGCAAGCTTAG"),      # real sequence in the clip (e.g. a breakpoint)
            rec(99, "15S20M", "GGGGGGGGGGGGGGG" + "ACGT" * 5),      # G run at the 5' end of a forward read
            rec(83, "20M15S", "ACGT" * 5 + "GGGGGGGGGGGGGGG"),      # G run at the 5' end of a reverse read
            rec(99, "20M3S", "ACGT" * 5 + "GGG"),                   # too short to judge
            rec(99, "35M", "ACGT" * 5 + "GGGGGGGGGGGGGGG"),         # aligned G run
        ]
        out, err = run(records)
        self.assertEqual(out, records)
        self.assertIn("hard-clipped: 0", err)


if __name__ == "__main__":
    unittest.main()

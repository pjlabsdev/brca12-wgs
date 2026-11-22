"""Read-QC notes on the sequencing platform (read-name format) and on missing duplicate marking."""
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import brca12_read_qc as qc  # noqa: E402

STATS = ("SN\traw total sequences:\t1000\nSN\treads mapped:\t990\nSN\treads properly paired:\t980\n"
         "SN\treads duplicated:\t{dup}\nSN\treads MQ0:\t3\nSN\terror rate:\t0.002\nSN\taverage length:\t150\n"
         "SN\taverage quality:\t35.0\nSN\tinsert size average:\t400.0\nIS\t390\t200\t0\t0\n")
MGI = ["E250160051L1C025R0180440{:04d}".format(i) for i in range(100)]
ILLUMINA = ["A00132:44:HCLL7DSXX:2:2417:{}:29872".format(16000 + i) for i in range(100)]
HEADER_RAW = "@HD\tVN:1.6\n@RG\tID:x\tSM:x\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\n@PG\tID:samtools\tCL:samtools sort\n"
HEADER_MD = HEADER_RAW + "@PG\tID:samtools.1\tCL:samtools markdup -s in.bam out.bam\n"


class ReadQcPlatform(unittest.TestCase):
    def run_qc(self, names, header, dup):
        d = tempfile.mkdtemp()
        paths = {}
        for key, text in (("stats", STATS.format(dup=dup)), ("names", "\n".join(names) + "\n"), ("header", header)):
            paths[key] = os.path.join(d, key)
            with open(paths[key], "w") as out:
                out.write(text)
        out = os.path.join(d, "read_qc.json")
        qc.main(["--stats", paths["stats"], "--read-names", paths["names"], "--header", paths["header"], "--out", out])
        return json.load(open(out))

    def test_mgi_reads_labelled_illumina_without_duplicate_marking(self):
        r = self.run_qc(MGI, HEADER_RAW, 0)
        self.assertEqual(r["readNamePlatform"], "MGI/DNBSEQ")
        self.assertEqual(r["readGroupPlatforms"], ["ILLUMINA"])
        self.assertFalse(r["duplicateMarkingInHeader"])
        text = " ".join(r["flags"])
        self.assertIn("MGI/DNBSEQ", text)
        self.assertIn("Read-Group nennt ILLUMINA", text)
        self.assertIn("keine Duplikatmarkierung", text)

    def test_illumina_reads_with_duplicate_marking_give_no_note(self):
        r = self.run_qc(ILLUMINA, HEADER_MD, 50)
        self.assertEqual(r["readNamePlatform"], "Illumina")
        self.assertTrue(r["duplicateMarkingInHeader"])
        self.assertEqual(r["flags"], [])

    def test_markdup_only_in_an_output_file_name_is_not_duplicate_marking(self):
        header = HEADER_RAW + "@PG\tID:samtools.2\tPN:samtools\tCL:samtools view -b -o /x/sample.markdup.bam -\n"
        r = self.run_qc(MGI, header, 0)
        self.assertFalse(r["duplicateMarkingInHeader"])
        self.assertIn("keine Duplikatmarkierung", " ".join(r["flags"]))

    def test_picard_and_sambamba_count_as_duplicate_marking(self):
        for pg in ("@PG\tID:MarkDuplicates\tPN:MarkDuplicates\tCL:MarkDuplicates INPUT=a.bam OUTPUT=b.bam\n",
                   "@PG\tID:sambamba\tPN:sambamba\tCL:/usr/bin/sambamba markdup -t 4 a.bam b.bam\n"):
            self.assertTrue(self.run_qc(ILLUMINA, HEADER_RAW + pg, 0)["duplicateMarkingInHeader"])

    def test_marked_header_but_no_duplicates_is_not_flagged(self):
        r = self.run_qc(ILLUMINA, HEADER_MD, 0)
        self.assertEqual(r["flags"], [])

    def test_too_few_or_mixed_names_give_no_platform(self):
        self.assertIsNone(self.run_qc(MGI[:10], HEADER_MD, 50)["readNamePlatform"])
        self.assertIsNone(self.run_qc(MGI[:50] + ILLUMINA[:50], HEADER_MD, 50)["readNamePlatform"])


if __name__ == "__main__":
    unittest.main()

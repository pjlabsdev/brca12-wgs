"""Cohort summary, contamination hint, XLSX export and sample sheet (bin/brca12_cohort.py etc.).

    python3 -m unittest discover -s tests
"""
import copy
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
import zipfile
from xml.etree import ElementTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)

from brca12lib import xlsx  # noqa: E402
import brca12_batch  # noqa: E402
import brca12_cohort  # noqa: E402

FIXTURE = os.path.join(ROOT, "tests", "fixtures", "brca12_results.sample.json")
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def sample(name, fraction, genotype="heterozygot"):
    with open(FIXTURE) as handle:
        result = json.load(handle)
    result = copy.deepcopy(result)
    result["sample"] = name
    variant = result["genes"]["BRCA1"]["smallVariants"][0]
    variant["alleleFraction"] = fraction
    variant["genotype"] = genotype
    return result


def write(directory, result):
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, "brca12_results.json")
    with open(path, "w") as handle:
        json.dump(result, handle)
    return path


def sheet_rows(path, index):
    with zipfile.ZipFile(path) as z:
        root = ElementTree.fromstring(z.read(f"xl/worksheets/sheet{index}.xml"))
    rows = []
    for row in root.find("m:sheetData", NS):
        cells = []
        for c in row:
            text = c.find("m:is/m:t", NS)
            value = c.find("m:v", NS)
            cells.append(text.text if text is not None else value.text if value is not None else None)
        rows.append(cells)
    return rows


class XlsxTest(unittest.TestCase):
    def test_workbook_is_valid_and_escapes_text(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "t.xlsx")
            xlsx.write(path, [("A/B:test", [["h1", "h2"], ["<x&y>", 1.5], [None, 3]]), ("A/B:test", [["x"]])])
            with zipfile.ZipFile(path) as z:
                self.assertIsNone(z.testzip())
                workbook = ElementTree.fromstring(z.read("xl/workbook.xml"))
                names = [s.get("name") for s in workbook.iter("{%s}sheet" % NS["m"])]
            self.assertEqual(names, ["A_B_test", "A_B_test (2)"])
            self.assertEqual(sheet_rows(path, 1)[1], ["<x&y>", "1.5"])

    def test_sample_export(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "s.xlsx")
            r = subprocess.run([sys.executable, os.path.join(BIN, "brca12_xlsx.py"), FIXTURE, out],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            rows = sheet_rows(out, 2)
            self.assertEqual(rows[0][:4], ["sample", "gene", "chrom", "pos_GRCh38"])
            self.assertEqual(len(rows), 1 + sum(len(g["smallVariants"]) for g in json.load(open(FIXTURE))["genes"].values()))


class CohortTest(unittest.TestCase):
    def test_contamination_hint_names_the_source(self):
        with tempfile.TemporaryDirectory() as d:
            write(os.path.join(d, "A"), sample("A", 0.08))
            write(os.path.join(d, "B"), sample("B", 0.48))
            write(os.path.join(d, "C"), sample("C", 0.51))
            out = os.path.join(d, "cohort")
            self.assertEqual(brca12_cohort.main(["--results", os.path.join(d, "*"), "--run-name", "R1", "--outdir", out]), 0)
            with open(os.path.join(out, "cohort_contamination.tsv")) as handle:
                hints = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([h["sample"] for h in hints], ["A"])
            self.assertEqual(hints[0]["possible_source"], "B, C")
            with open(os.path.join(out, "cohort_samples.tsv")) as handle:
                samples = {r["sample"]: r for r in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(samples["A"]["contamination_hints"], "1")
            self.assertEqual(samples["B"]["contamination_hints"], "0")
            with open(os.path.join(out, "cohort_variants.tsv")) as handle:
                variants = list(csv.DictReader(handle, delimiter="\t"))
            first = [v for v in variants if v["hgvsc"] == variants[0]["hgvsc"]]
            self.assertTrue(all(v["cohort_carriers"] == "3" for v in first))
            self.assertTrue(zipfile.is_zipfile(os.path.join(out, "cohort.xlsx")))

    def test_no_hint_without_a_clear_carrier(self):
        with tempfile.TemporaryDirectory() as d:
            write(os.path.join(d, "A"), sample("A", 0.08))
            write(os.path.join(d, "B"), sample("B", 0.10))
            out = os.path.join(d, "cohort")
            brca12_cohort.main(["--results", os.path.join(d, "*"), "--outdir", out])
            with open(os.path.join(out, "cohort_contamination.tsv")) as handle:
                self.assertEqual(len(handle.read().strip().splitlines()), 1)

    def test_duplicate_sample_names_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            write(os.path.join(d, "A"), sample("X", 0.5))
            write(os.path.join(d, "B"), sample("X", 0.5))
            with self.assertRaises(SystemExit):
                brca12_cohort.main(["--results", os.path.join(d, "*"), "--outdir", os.path.join(d, "c")])


class SamplesheetTest(unittest.TestCase):
    def test_parse_and_commands(self):
        with tempfile.TemporaryDirectory() as d:
            sheet = os.path.join(d, "s.csv")
            with open(sheet, "w") as handle:
                handle.write("sample,reads_r1,reads_r2,bam\nS01,a_R1.fq.gz,a_R2.fq.gz,\nS02,,,b.cram\n")
            samples = brca12_batch.read_samplesheet(sheet)
            self.assertEqual([s["sample"] for s in samples], ["S01", "S02"])
            cmd = brca12_batch.command(samples[1], d, "nextflow", ["-profile", "conda"])
            self.assertIn("--bam", cmd)
            self.assertEqual(cmd[-2:], ["-profile", "conda"])
            self.assertIn(os.path.join(d, "S02"), cmd)

    def test_invalid_rows(self):
        with tempfile.TemporaryDirectory() as d:
            for content in ("sample,reads_r1,reads_r2,bam\nS01,a.fq.gz,,\n",
                            "sample,reads_r1,reads_r2,bam\nS01,a,b,c.bam\n",
                            "sample,reads_r1,reads_r2,bam\nS 01,a,b,\n",
                            "sample,reads_r1,reads_r2,bam\nS01,a,b,\nS01,c,d,\n"):
                sheet = os.path.join(d, "s.csv")
                with open(sheet, "w") as handle:
                    handle.write(content)
                with self.assertRaises(SystemExit):
                    brca12_batch.read_samplesheet(sheet)


if __name__ == "__main__":
    unittest.main()

"""PREFIX.tsv, PREFIX.findings.tsv and PREFIX.exons.tsv (bin/brca12lib/tsv.py), in the form of hla-wgs."""
import json
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))

from brca12lib import page, tsv  # noqa: E402

with open(os.path.join(ROOT, "tests", "fixtures", "brca12_results.sample.json")) as handle:
    RESULT = json.load(handle)


def copy(result=RESULT):
    return json.loads(json.dumps(result))


def parse(text):
    lines = text.rstrip("\n").split("\n")
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    for row in rows:
        assert len(row) == len(header)
    return header, rows


class GeneTable(unittest.TestCase):
    def test_one_line_per_gene_with_the_result_of_the_page(self):
        header, rows = parse(tsv.genes_tsv(RESULT))
        self.assertEqual(tuple(header), tsv.GENE_COLUMNS)
        self.assertEqual([r["gene"] for r in rows], ["BRCA1", "BRCA2"])
        levels = page.gene_levels(RESULT)
        for row in rows:
            self.assertEqual(row["result"], levels[row["gene"]])
        brca1 = rows[0]
        self.assertEqual(brca1["result"], "attention")
        self.assertEqual(brca1["structural"], "1")
        self.assertEqual(brca1["coverage"], "limited")
        self.assertEqual(brca1["exons_uncertain"], "1")

    def test_clean_gene_and_missing_values(self):
        result = copy()
        for gene in result["genes"].values():
            gene["smallVariants"] = []
        result.update(cnvs=[], svs=[])
        for exon in result["exonCoverage"]:
            exon.update(meanDepth=35.0, status="ok")
        result["lowCoverageIntervals"] = []
        _, rows = parse(tsv.genes_tsv(result))
        for row in rows:
            self.assertEqual((row["result"], row["attention"], row["unclear"], row["coverage"]), ("none", "0", "0", "ok"))


class FindingTable(unittest.TestCase):
    def test_every_finding_with_level_and_certainty(self):
        result = copy()
        exon = result["exonCoverage"][0]
        result["exonDosageHints"] = [{"gene": exon["gene"], "exon": exon["label"], "ratio": 0.27,
                                      "type": "möglicher Verlust", "flags": []}]
        header, rows = parse(tsv.findings_tsv(result))
        self.assertEqual(tuple(header), tsv.FINDING_COLUMNS)
        small = [r for r in rows if r["kind"] == "small"]
        self.assertEqual(len(small), sum(len(g["smallVariants"]) for g in RESULT["genes"].values()))
        plp = next(r for r in small if r["clinvar"] == "Pathogenic")
        self.assertEqual((plp["level"], plp["reason"], plp["certainty"]), ("attention", "clinvar_pathogenic", "ok"))
        self.assertEqual(plp["reads"], "18/38")
        self.assertEqual(plp["genotype"], "heterozygous")
        sv = next(r for r in rows if r["kind"] == "sv")
        self.assertEqual((sv["level"], sv["certainty"]), ("attention", "ok"))
        self.assertIn("orthogonal", sv["flags"])
        hint = next(r for r in rows if r["kind"] == "exon_dosage")
        self.assertEqual((hint["level"], hint["certainty"], hint["start"]), ("hint", "hint", str(exon["start"])))
        self.assertEqual(hint["change"], "possible loss ratio 0.27")
        # same count of attention and unclear entries as the page
        entries = page.notable(RESULT)
        for level in ("attention", "unclear"):
            self.assertEqual(sum(1 for r in rows if r["level"] == level), sum(1 for e in entries if e[0] == level))

    def test_common_allele_ba1_flags_and_english_codes(self):
        result = copy()
        result["insertionSignals"] = [{"gene": "BRCA1", "chrom": "chr17", "pos": 43128247, "end": 43128372,
                                       "type": "Insertion", "zygosity": "heterozygot", "realignedReads": 12,
                                       "label": "zusätzliche 126-bp-Kopie vor BRCA1 Exon 1 (häufiges Allel)",
                                       "flags": []}]
        variant = result["genes"]["BRCA1"]["smallVariants"][0]
        variant["clinvar"]["significance"] = "Benign"
        variant["acmgHints"] = [{"code": "BA1 (Näherung)", "textDE": "", "textEN": ""},
                                {"code": "PP3/BP4 offen", "textDE": "", "textEN": ""}]
        variant["qualityFlags"] = ["bekannte Schwachstelle"]
        _, rows = parse(tsv.findings_tsv(result))
        ins = next(r for r in rows if r["kind"] == "insertion_signal")
        self.assertEqual((ins["level"], ins["reason"], ins["certainty"]), ("common", "common_allele", "hint"))
        self.assertIn("common allele", ins["change"])
        self.assertEqual(ins["reads"], "realignedReads 12")
        small = next(r for r in rows if r["kind"] == "small" and r["gene"] == "BRCA1")
        self.assertEqual((small["level"], small["reason"], small["certainty"]), ("none", "ba1", "limited"))
        self.assertEqual(small["acmg"], "BA1 (approximation); PP3/BP4 open")

    def test_tabs_and_newlines_cannot_break_a_line(self):
        result = copy()
        result["svs"][0]["flags"] = ["a\tb\nc"]
        _, rows = parse(tsv.findings_tsv(result))
        self.assertIn("a b c", next(r for r in rows if r["kind"] == "sv")["flags"])


class ExonTable(unittest.TestCase):
    def test_one_line_per_exon_with_status_and_reasons(self):
        header, rows = parse(tsv.exons_tsv(RESULT))
        self.assertEqual(tuple(header), tsv.EXON_COLUMNS)
        self.assertEqual(len(rows), len(RESULT["exonCoverage"]))
        uncertain = [r for r in rows if r["status"] == "uncertain"]
        self.assertTrue(uncertain)
        self.assertIn("x <", uncertain[0]["reasons"])


if __name__ == "__main__":
    unittest.main()

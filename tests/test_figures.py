"""SVG figures in the HTML report (bin/brca12lib/figures.py)."""
import json
import os
import sys
import unittest
from xml.etree import ElementTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))

from brca12lib import figures, render  # noqa: E402

with open(os.path.join(ROOT, "tests", "fixtures", "brca12_results.sample.json")) as handle:
    RESULT = json.load(handle)


def svg_part(text):
    return text[:text.index("</svg>") + 6]


class FiguresTest(unittest.TestCase):
    def test_coverage_svg_is_valid_and_marks_uncertain_exons(self):
        genes = {e["gene"] for e in RESULT["exonCoverage"]}
        for gene in genes:
            svg = figures.coverage_svg(RESULT, gene, "de")
            root = ElementTree.fromstring(svg_part(svg))
            self.assertTrue(root.tag.endswith("svg"))
            bars = [el for el in root.iter() if el.tag.endswith("rect")]
            exons = [e for e in RESULT["exonCoverage"] if e["gene"] == gene]
            self.assertEqual(len(bars), len(exons))
            uncertain = [e for e in exons if e["status"] != "ok"]
            self.assertEqual(sum(1 for b in bars if b.get("fill") == figures.WARN), len(uncertain))

    def test_dosage_flags_outside_band(self):
        result = json.loads(json.dumps(RESULT))
        gene = result["exonCoverage"][0]["gene"]
        result["exonCoverage"][0]["dosageRatio"] = 0.48
        svg = figures.dosage_svg(result, gene, "en")
        root = ElementTree.fromstring(svg_part(svg))
        flagged = [c for c in root.iter() if c.tag.endswith("circle") and c.get("fill") == figures.WARN]
        self.assertEqual(len(flagged), 1)

    def test_html_report_contains_figures(self):
        html = render.html_report(RESULT, "en")
        self.assertIn("<h2>Figures</h2>", html)
        self.assertIn("<svg", html)
        self.assertNotIn("<script", html)
        self.assertIn("<h2>Abbildungen</h2>", render.html_report(RESULT, "de"))

    def test_no_figures_without_exons(self):
        self.assertEqual(figures.figures_html({"exonCoverage": []}, "de"), "")


if __name__ == "__main__":
    unittest.main()

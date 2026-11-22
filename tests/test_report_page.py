"""HTML report page (bin/brca12lib/page.py): both languages, two tabs, no script."""
import json
import os
import re
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))

from brca12lib import page, render  # noqa: E402

with open(os.path.join(ROOT, "tests", "fixtures", "brca12_results.sample.json")) as handle:
    RESULT = json.load(handle)


def main_part(text, lang):
    part = text[text.index(f"<main lang='{lang}'"):]
    return part[:part.index("</main>")]


class ReportPageTest(unittest.TestCase):
    def test_both_languages_in_one_file_without_script(self):
        for first, other in (("de", "en"), ("en", "de")):
            text = render.html_report(RESULT, first)
            self.assertNotIn("<script", text)
            self.assertLess(text.index(f"<main lang='{first}'"), text.index(f"<main lang='{other}'"))
            self.assertIn(f"id='lang-{first}' aria-label", text)
            self.assertRegex(text, rf"id='lang-{first}'[^>]* checked>")
            self.assertNotRegex(text, rf"id='lang-{other}'[^>]* checked>")
        text = render.html_report(RESULT, "en")
        self.assertIn("Research use only (RUO)", main_part(text, "en"))
        self.assertIn("Nur für Forschung (RUO)", main_part(text, "de"))

    def test_tabs_result_and_details(self):
        text = main_part(render.html_report(RESULT, "de"), "de")
        result_tab = text[text.index("<section class='tab-r'>"):text.index("<section class='tab-d'>")]
        details_tab = text[text.index("<section class='tab-d'>"):]
        self.assertIn("Zu prüfen", result_tab)
        self.assertIn("ClinVar pathogen / wahrscheinlich pathogen", result_tab)
        self.assertNotIn("Tool-Versionen", result_tab)
        for title in ("BRCA1 SNVs/Indels", "Coverage je Exon", "Tool-Versionen", "Grenzen", "Abbildungen"):
            self.assertIn(title, details_tab)

    def test_print_shows_both_tabs_and_page_numbers(self):
        css = page.CSS
        printed = css[css.index("@media print"):]
        self.assertIn("#tab-r:checked~main section.tab-d,#tab-d:checked~main section.tab-r{display:block}", printed)
        self.assertIn("counter(page) ' / ' counter(pages)", css)

    def test_header_names_program_doi_and_research_page(self):
        result = dict(RESULT, pipeline={"name": "brca12-wgs", "version": "9.9.9", "doi": page.DOI})
        text = render.html_report(result, "en")
        self.assertIn("brca12-wgs 9.9.9", text)
        self.assertIn(f"https://doi.org/{page.DOI}", text)
        self.assertIn(page.RESEARCH_URL, text)
        self.assertIn("font/woff2;base64,", text)

    def test_low_depth_notice_and_quality(self):
        result = json.loads(json.dumps(RESULT))
        for exon in result["exonCoverage"]:
            exon.update(meanDepth=35.0, status="ok")
        result["readQC"] = None
        self.assertEqual(page.quality(result), (True, []))
        text = main_part(render.html_report(result, "de"), "de")
        self.assertIn("bestanden", text)
        self.assertNotIn("Geringe Tiefe", text)
        for exon in result["exonCoverage"]:
            if exon["gene"] == "BRCA2":
                exon["meanDepth"] = 11.0
        passed, reasons = page.quality(result)
        self.assertFalse(passed)
        self.assertEqual(reasons, [("q_depth", {"gene": "BRCA2", "depth": 11.0})])
        text = main_part(render.html_report(result, "de"), "de")
        self.assertIn("Geringe Tiefe: BRCA2 hat im Mittel 11x", text)
        self.assertIn("eingeschränkt", text)
        self.assertIn("Low depth: BRCA2", main_part(render.html_report(result, "de"), "en"))

    def test_nothing_to_review_keeps_the_caveat(self):
        result = json.loads(json.dumps(RESULT))
        for gene in result["genes"].values():
            gene["smallVariants"] = []
        result.update(cnvs=[], svs=[])
        self.assertEqual(page.notable_rows(result, "de"), [])
        text = main_part(render.html_report(result, "en"), "en")
        self.assertIn("does not rule out a hereditary predisposition", text)

    def test_sample_name_cannot_break_the_style_block(self):
        result = dict(RESULT, sample="a'}</style><b>x")
        text = render.html_report(result, "en")
        self.assertEqual(len(re.findall(r"</style>", text)), 1)
        self.assertNotIn("<b>x", text)


if __name__ == "__main__":
    unittest.main()

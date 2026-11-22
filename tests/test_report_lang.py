"""Report languages (German, English); the PDF step is gone (the HTML page is the report)."""
import json
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)

from brca12lib import i18n, render  # noqa: E402

with open(os.path.join(ROOT, "tests", "fixtures", "brca12_results.sample.json")) as handle:
    RESULT = json.load(handle)


class LanguageTest(unittest.TestCase):
    def test_only_german_and_english(self):
        self.assertEqual(i18n.LANGUAGES, ("de", "en"))
        self.assertFalse(hasattr(i18n, "UI_RU"))
        self.assertFalse(hasattr(render, "DISCLAIMER_RU"))

    def test_german_and_english_unchanged(self):
        self.assertEqual(i18n.ui("title", "de"), "BRCA1/BRCA2-Bericht")
        self.assertEqual(i18n.tr("heterozygot", "en"), "heterozygous")
        self.assertEqual(i18n.tr("heterozygot", "de"), "heterozygot")

    def test_report_script_rejects_ru(self):
        r = subprocess.run([sys.executable, os.path.join(BIN, "brca12_report.py"), "--lang", "ru"],
                           capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("invalid choice", r.stderr)


class MarkdownBothLanguagesTest(unittest.TestCase):
    def test_markdown_has_both_languages_chosen_one_first(self):
        text = render.markdown_both(RESULT, "en")
        en, de = text.index(i18n.ui("title", "en")), text.index(i18n.ui("title", "de"))
        self.assertLess(en, de)
        text = render.markdown_both(RESULT, "de")
        self.assertLess(text.index(i18n.ui("title", "de")), text.index(i18n.ui("title", "en")))


class PageTitleTest(unittest.TestCase):
    def test_html_title_carries_sample_and_ruo(self):
        self.assertIn(f"<title>brca12-wgs · {RESULT['sample']}</title>", render.html_report(RESULT, "en"))

    def test_no_pdf_step_left(self):
        self.assertFalse(os.path.exists(os.path.join(BIN, "brca12_pdf.sh")))
        main = open(os.path.join(ROOT, "main.nf")).read()
        self.assertNotIn("REPORT_PDF", main)
        self.assertNotIn("report_pdf", open(os.path.join(ROOT, "nextflow.config")).read())


if __name__ == "__main__":
    unittest.main()

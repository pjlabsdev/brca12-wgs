"""HTML result page (bin/brca12lib/page.py), built after the hla-wgs result page."""
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


def copy(result=RESULT):
    return json.loads(json.dumps(result))


def main_part(text, lang):
    part = text[text.index(f"<main lang='{lang}'"):]
    return part[:part.index("</main>")]


def tab(text, lang, name):
    part = main_part(text, lang)
    start = part.index(f"<section class='tab-{name}'>")
    return part[start:part.index("</section>", start)]


def clean():
    """The sample result without any finding and with good coverage."""
    result = copy()
    for gene in result["genes"].values():
        gene["smallVariants"] = []
    result.update(cnvs=[], svs=[], readQC=None)
    for exon in result["exonCoverage"]:
        exon.update(meanDepth=35.0, status="ok")
    return result


class ReportPageTest(unittest.TestCase):
    def test_both_languages_english_first_by_default(self):
        text = render.html_report(RESULT, "en")
        self.assertIn("<body data-lang='en' data-tab='r'>", text)
        self.assertLess(text.index("<main lang='en'"), text.index("<main lang='de'"))
        self.assertIn("Research use only (RUO)", main_part(text, "en"))
        self.assertIn("Nur für Forschung (RUO)", main_part(text, "de"))
        self.assertIn("<body data-lang='de' data-tab='r'>", render.html_report(RESULT, "de"))

    def test_interactive_like_hla_wgs(self):
        text = render.html_report(RESULT, "en")
        for marker in ("nav class='lang seg'", "nav class='tabs seg'", "data-tab='d'", "class=\"tip\"", "<script>",
                       "localStorage.setItem('brca12-wgs-lang'", "beforeprint", "class='copy'", "<details>"):
            self.assertIn(marker, text)
        self.assertIn("body[data-lang=de] main.l-en", page.CSS)
        self.assertIn("counter(page) ' / ' counter(pages)", page.CSS)

    def test_result_tab_and_details_per_gene(self):
        text = render.html_report(RESULT, "de")
        result_tab = tab(text, "de", "r")
        details = tab(text, "de", "d")
        for title in ("Gesamtaussage", "Zu prüfen", "ClinVar pathogen / wahrscheinlich pathogen", "Abdeckung",
                      "Methode und Grenzen"):
            self.assertIn(title, result_tab)
        self.assertNotIn("Tool-Versionen", result_tab)
        self.assertEqual(details.count("<details>"), 3)  # BRCA1, BRCA2, run
        brca1 = details[details.index("<summary><b>BRCA1</b>"):details.index("<summary><b>BRCA2</b>")]
        self.assertIn("BRCA1 SNVs/Indels", brca1)
        self.assertNotIn("BRCA2 SNVs/Indels", brca1)
        self.assertIn("43093000", brca1)
        self.assertNotIn("32338000", brca1)
        self.assertIn("Tool-Versionen", details)

    def test_report_number_is_stable_and_changes_with_the_result(self):
        a = page.report_number(RESULT)
        self.assertRegex(a, r"^BRCA-[0-9A-F]{8}$")
        later = dict(copy(), generatedAt="2030-01-01T00:00:00+00:00", nextflowSession="x")
        self.assertEqual(page.report_number(later), a)
        changed = copy()
        changed["genes"]["BRCA1"]["smallVariants"][0]["pos"] += 1
        self.assertNotEqual(page.report_number(changed), a)
        self.assertIn(a, render.html_report(RESULT, "en"))

    def test_overall_statement(self):
        self.assertEqual(page.verdict(RESULT), ("attention", 2))  # ClinVar P (also LoF, counted once) and a PASS SV
        text = main_part(render.html_report(RESULT, "en"), "en")
        self.assertIn("verdict v-attention", text)
        result = clean()
        self.assertEqual(page.verdict(result), ("none", 0))
        self.assertEqual(page.quality(result), (True, []))
        text = render.html_report(result, "en")
        self.assertIn("no sign of a hereditary predisposition in the categories checked", main_part(text, "en"))
        self.assertIn("does not rule out a predisposition", main_part(text, "en"))
        self.assertIn("keinen Hinweis auf eine erbliche Veranlagung", main_part(text, "de"))
        self.assertNotIn("The result is limited", main_part(text, "en"))
        vus = copy(RESULT["genes"]["BRCA1"]["smallVariants"][0])
        vus["clinvar"]["significance"] = "Uncertain_significance"
        vus["acmgHints"] = []
        result["genes"]["BRCA1"]["smallVariants"] = [vus]
        self.assertEqual(page.verdict(result), ("unclear", 1))

    def test_cancer_types_table(self):
        names = [row[0] for row in page.cancer_rows(clean(), "en")]
        self.assertEqual(names, ["Breast cancer (women)", "Breast cancer (men)",
                                 "Ovarian cancer (with fallopian tube and peritoneum)", "Prostate cancer",
                                 "Pancreatic cancer", "Testicular cancer"])
        rows = {row[0]: row for row in page.cancer_rows(clean(), "en")}
        self.assertEqual(rows["Breast cancer (women)"][3], "no finding in BRCA1/BRCA2 that raises this risk")
        self.assertEqual(rows["Testicular cancer"][1:], ["no established link", "no established link",
                                                         "BRCA1/BRCA2 say nothing about this"])
        # a finding in one gene only marks the cancers linked to that gene and names it
        result = clean()
        variant = copy(RESULT["genes"]["BRCA1"]["smallVariants"][0])
        variant["clinvar"]["significance"] = "Pathogenic"
        variant["gene"] = "BRCA2"
        result["genes"]["BRCA2"]["smallVariants"] = [variant]
        self.assertEqual(page.gene_levels(result), {"BRCA1": "none", "BRCA2": "attention"})
        rows = {row[0]: row for row in page.cancer_rows(result, "en")}
        self.assertTrue(rows["Prostate cancer"][3].startswith("finding in BRCA2 that may raise this risk"))
        self.assertEqual(rows["Testicular cancer"][3], "BRCA1/BRCA2 say nothing about this")
        variant["clinvar"]["significance"] = "Uncertain_significance"
        variant["acmgHints"] = []
        rows = {row[0]: row for row in page.cancer_rows(result, "de")}
        self.assertTrue(rows["Eierstockkrebs (mit Eileiter und Bauchfell)"][3].startswith("nur Einträge mit unklarer Bedeutung in BRCA2"))
        text = render.html_report(result, "en")
        self.assertIn("Cancer types and BRCA1/BRCA2", tab(text, "en", "r"))
        self.assertIn("Krebsarten und BRCA1/BRCA2", tab(text, "de", "r"))
        self.assertIn("are not checked by brca12-wgs", tab(text, "en", "r"))

    def test_low_depth_notice_and_limited_quality(self):
        result = clean()
        for exon in result["exonCoverage"]:
            if exon["gene"] == "BRCA2":
                exon["meanDepth"] = 11.0
        passed, reasons = page.quality(result)
        self.assertFalse(passed)
        self.assertEqual(reasons, [("q_depth", {"gene": "BRCA2", "depth": 11.0})])
        text = render.html_report(result, "de")
        self.assertIn("Geringe Tiefe: BRCA2 hat im Mittel 11x", main_part(text, "de"))
        self.assertIn("Die Aussagekraft ist eingeschränkt", main_part(text, "de"))
        self.assertIn("Low depth: BRCA2", main_part(text, "en"))

    def test_read_check_notes_are_shown_and_translated(self):
        source = open(os.path.join(ROOT, "bin", "brca12_read_qc.py")).read()
        flags = re.findall(r'flags\.append\(((?:\s*"[^"]*")+)\)', source)
        notes = ["".join(re.findall(r'"([^"]*)"', f)) for f in flags]
        platform_notes = [n for n in notes if "MGI" in n or "Duplikatmarkierung" in n]
        self.assertEqual(len(platform_notes), 3)
        result = clean()
        result["readQC"] = {"flags": platform_notes}
        text = render.html_report(result, "de")
        for note in platform_notes:
            self.assertIn(note, tab(text, "de", "r"))
            self.assertNotIn(note, main_part(text, "en"))
        self.assertIn("validated on Illumina data", tab(text, "en", "r"))
        self.assertIn("No duplicate marking found", tab(text, "en", "r"))

    def test_header_names_program_doi_and_research_page(self):
        result = dict(copy(), pipeline={"name": "brca12-wgs", "version": "9.9.9", "doi": page.DOI})
        text = render.html_report(result, "en")
        self.assertIn("brca12-wgs 9.9.9", text)
        self.assertIn(f"https://doi.org/{page.DOI}", text)
        self.assertIn(page.RESEARCH_URL, text)
        self.assertIn(page.SOURCE_URL, text)
        self.assertIn("font/woff2;base64,", text)

    def test_sample_name_cannot_break_the_page(self):
        result = dict(copy(), sample="a'}</style><b>x")
        text = render.html_report(result, "en")
        self.assertEqual(text.count("</style>"), 1)
        self.assertNotIn("<b>x", text)


if __name__ == "__main__":
    unittest.main()

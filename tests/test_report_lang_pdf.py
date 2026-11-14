"""Russian report language and the optional PDF step (bin/brca12_pdf.sh)."""
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


class RussianTest(unittest.TestCase):
    def test_every_ui_key_and_value_has_a_russian_text(self):
        self.assertEqual(set(i18n.UI) - set(i18n.UI_RU), set())
        self.assertEqual(set(i18n.EXACT) - set(i18n.EXACT_RU), set())
        self.assertEqual(len(i18n.PATTERNS), len(i18n.PATTERNS_RU))

    def test_report_in_russian(self):
        html = render.html_report(RESULT, "ru")
        self.assertIn("<html lang='ru'>", html)
        self.assertIn("Отчёт BRCA1/BRCA2", html)
        self.assertIn("Только для исследований (RUO)", html)
        self.assertIn("гетерозигота", html)
        self.assertNotIn("heterozygot<", html)
        md = render.markdown(RESULT, "ru")
        self.assertIn("Ограничения", md)

    def test_unknown_values_fall_back_to_english(self):
        self.assertEqual(i18n.tr("Exon 3 (Legacy 4)", "ru"), "Экзон 3 (устар. нумерация 4)")
        self.assertEqual(i18n.tr("Hinweis aus Leseenden, kein Variantenaufruf; Bestätigung mit orthogonaler Methode nötig", "ru")[:8], "указание")

    def test_german_and_english_unchanged(self):
        self.assertEqual(i18n.ui("title", "de"), "BRCA1/BRCA2-Bericht")
        self.assertEqual(i18n.tr("heterozygot", "en"), "heterozygous")
        self.assertEqual(i18n.tr("heterozygot", "de"), "heterozygot")

    def test_report_script_accepts_ru(self):
        out = subprocess.run([sys.executable, os.path.join(BIN, "brca12_report.py"), "--help"],
                             capture_output=True, text=True).stdout
        self.assertIn("ru", out)


class PdfTest(unittest.TestCase):
    def run_pdf(self, chrome, path_env):
        with tempfile.TemporaryDirectory() as d:
            html = os.path.join(d, "r.html")
            with open(html, "w") as handle:
                handle.write(render.html_report(RESULT, "de"))
            pdf = os.path.join(d, "r.pdf")
            env = dict(os.environ, PATH=path_env)
            r = subprocess.run(["bash", os.path.join(BIN, "brca12_pdf.sh"), html, pdf, chrome],
                               capture_output=True, text=True, env=env)
            return r, os.path.exists(pdf)

    def test_skipped_without_chrome(self):
        r, exists = self.run_pdf("", "/usr/bin:/bin")
        self.assertEqual(r.returncode, 0)
        self.assertFalse(exists)
        self.assertIn("skipped", r.stderr)

    def test_written_with_a_chrome(self):
        with tempfile.TemporaryDirectory() as d:
            fake = os.path.join(d, "fakechrome")
            with open(fake, "w") as handle:
                handle.write("#!/bin/sh\nfor a in \"$@\"; do case \"$a\" in --print-to-pdf=*) printf '%%PDF-1.4' > \"${a#--print-to-pdf=}\";; esac; done\n")
            os.chmod(fake, 0o755)
            r, exists = self.run_pdf(fake, "/usr/bin:/bin")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(exists)

    def test_html_title_carries_sample_and_ruo(self):
        self.assertIn(f"· {RESULT['sample']} · RUO</title>", render.html_report(RESULT, "en"))


if __name__ == "__main__":
    unittest.main()

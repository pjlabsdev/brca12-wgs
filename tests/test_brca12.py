"""Synthetische Ende-zu-Ende-Tests für die Auswerteskripte (ohne Bioinformatik-Tools).

    python3 -m unittest discover -s tests
"""
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)

from brca12lib import structural, variants  # noqa: E402
from brca12lib.vcf import genotype_label  # noqa: E402

CSQ_FIELDS = ("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|"
              "VARIANT_CLASS|CANONICAL|MANE_SELECT|PICK|MAX_AF|gnomADe_AF|gnomADg_AF|"
              "ClinVar|ClinVar_CLNSIG|ClinVar_CLNREVSTAT|ClinVar_CLNDN|ClinVar_CLNSIGCONF")


def write(path, text, gz=False):
    opener = gzip.open if gz else open
    with opener(path, "wt") as handle:
        handle.write(text)
    return path


def csq(allele, consequence, gene, transcript, exon, hgvsc, hgvsp, max_af="", clnsig="", rev=""):
    values = [allele, consequence, "HIGH", gene, "ENSG", "Transcript", transcript, "protein_coding", exon, "",
              hgvsc, hgvsp, "SNV", "YES", "NM_x", "1", max_af, "", "", "", clnsig, rev, "Hereditary_cancer", ""]
    return "|".join(values)


class PipelineScriptsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        d = self.tmp
        self.fai = write(f"{d}/ref.fa.fai", "chr13\t114364328\t0\t60\t61\nchr17\t83257441\t0\t60\t61\n")
        self.regions = write(f"{d}/regions.bed", "chr17\t43044294\t43170327\tBRCA1\nchr13\t32315076\t32400268\tBRCA2\n")
        self.exons = write(f"{d}/exons.bed",
                           "chr17\t43045000\t43045100\tBRCA1_ex23\t.\t-\n"
                           "chr17\t43091435\t43094860\tBRCA1_ex10\t.\t-\n"
                           "chr13\t32336264\t32341196\tBRCA2_ex11\t.\t+\n"
                           "chr13\t32398161\t32400268\tBRCA2_ex27\t.\t+\n")

    def run_script(self, *args):
        self.last_args = args
        return subprocess.run([sys.executable, *args], capture_output=True, text=True)

    def prepare(self):
        d = self.tmp
        r = self.run_script(f"{BIN}/brca12_prepare.py", "--fai", self.fai, "--regions", self.regions,
                            "--exons", self.exons, "--out-check", f"{d}/check.json",
                            "--out-calling", f"{d}/calling.bed", "--out-sv", f"{d}/sv.bed",
                            "--out-exact", f"{d}/exact.bed", "--out-exons", f"{d}/exons.norm.bed")
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_rejects_grch37(self):
        write(self.fai, "chr13\t115169878\t0\t60\t61\nchr17\t81195210\t0\t60\t61\n")
        r = self.run_script(f"{BIN}/brca12_prepare.py", "--fai", self.fai, "--regions", self.regions,
                            "--exons", self.exons, "--out-check", "x", "--out-calling", "x", "--out-sv", "x",
                            "--out-exact", "x", "--out-exons", "x")
        self.assertEqual(r.returncode, 2)
        self.assertIn("GRCh37", r.stderr)

    def test_rejects_region_outside_brca(self):
        write(self.regions, "chr17\t1000\t2000\tBRCA1\nchr13\t32315076\t32400268\tBRCA2\n")
        r = self.run_script(f"{BIN}/brca12_prepare.py", "--fai", self.fai, "--regions", self.regions,
                            "--exons", self.exons, "--out-check", "x", "--out-calling", "x", "--out-sv", "x",
                            "--out-exact", "x", "--out-exons", "x")
        self.assertEqual(r.returncode, 2)

    def test_bundled_assets_match_mane(self):
        d = self.tmp
        r = self.run_script(f"{BIN}/brca12_prepare.py", "--fai", self.fai,
                            "--regions", os.path.join(ROOT, "assets", "brca12_regions.bed"),
                            "--exons", os.path.join(ROOT, "assets", "brca12_exons.bed"),
                            "--out-check", f"{d}/c.json", "--out-calling", f"{d}/a", "--out-sv", f"{d}/b",
                            "--out-exact", f"{d}/e", "--out-exons", f"{d}/x")
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(f"{d}/c.json") as handle:
            self.assertEqual(json.load(handle)["exonCount"], {"BRCA1": 23, "BRCA2": 27})

    def test_last_exon_lof_is_restricted(self):
        hints = variants.acmg_hints("BRCA2", {"Consequence": "stop_gained", "EXON": "27/27",
                                              "HGVSp": "ENSP:p.Lys3326Ter"}, 0.005)
        codes = [h["code"] for h in hints]
        self.assertIn("PVS1 (eingeschränkt)", codes)
        self.assertIn("Bekannte Ausnahme", codes)
        self.assertIn("BA1 (Kandidat)", codes)

    def test_end_to_end_report(self):
        self.prepare()
        d = self.tmp
        header = ("##fileformat=VCFv4.2\n"
                  f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {CSQ_FIELDS}">\n'
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        brca1 = write(f"{d}/b1.vcf", header +
                      "chr17\t43093000\t.\tAG\tA\t500\tPASS\tMQ=60;CSQ=" +
                      csq("-", "frameshift_variant", "BRCA1", "ENST00000357654", "10/23",
                          "c.3000del", "p.Ser1000fs", "", "Pathogenic", "reviewed_by_expert_panel") +
                      "\tGT:AD:DP:GQ\t0/1:20:18:38:99\n"
                      "chr17\t43045050\t.\tC\tT\t50\tPASS\tMQ=60;CSQ=" +
                      csq("T", "missense_variant", "BRCA1", "ENST00000357654", "23/23", "c.5500G>A",
                          "p.Ala1834Thr", "0.0004", "Uncertain_significance",
                          "criteria_provided,_single_submitter") +
                      "\tGT:AD:DP:GQ\t0/1:9:3:12:15\n")
        brca2 = write(f"{d}/b2.vcf", header +
                      "chr13\t32338000\t.\tA\tG\t900\tPASS\tMQ=60;CSQ=" +
                      csq("G", "synonymous_variant", "BRCA2", "ENST00000380152", "11/27", "c.4000A>G",
                          "p.=", "0.3", "Benign", "reviewed_by_expert_panel") +
                      "\tGT:AD:DP:GQ\t1/1:0:35:35:99\n")
        manta = write(f"{d}/manta.vcf",
                      "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                      "chr17\t43090000\tMantaDEL:1\tN\t<DEL>\t200\tPASS\tSVTYPE=DEL;END=43096000;SVLEN=-6000\tGT:PR:SR\t0/1:20,8:15,6\n"
                      "chr2\t1000\tMantaDEL:2\tN\t<DEL>\t200\tPASS\tSVTYPE=DEL;END=5000\tGT:PR:SR\t0/1:20,8:15,6\n"
                      "chr5\t100\tMantaBND:3\tN\tN[chr13:32350000[\t100\tMinQUAL\tSVTYPE=BND;IMPRECISE\tGT:PR:SR\t0/1:10,4:.\n")
        cns = write(f"{d}/s.call.cns", "chromosome\tstart\tend\tgene\tlog2\tcn\tdepth\tprobes\tweight\n"
                    "chr17\t43000000\t43200000\t-\t-0.98\t1\t15\t200\t100\n"
                    "chr13\t30000000\t35000000\t-\t0.01\t2\t30\t5000\t100\n")
        cnr = write(f"{d}/s.cnr", "chromosome\tstart\tend\tgene\tdepth\tlog2\tweight\n"
                    "chr17\t43090000\t43091000\t-\t15\t-1.0\t1\nchr17\t43100000\t43101000\t-\t16\t-0.95\t1\n")
        q20 = write(f"{d}/q20.bed.gz",
                    "chr13\t32336264\t32341196\tBRCA2_ex11\t32.0\nchr13\t32398161\t32400268\tBRCA2_ex27\t31.0\n"
                    "chr17\t43045000\t43045100\tBRCA1_ex23\t8.0\nchr17\t43091435\t43094860\tBRCA1_ex10\t16.0\n", gz=True)
        q0 = write(f"{d}/q0.bed.gz",
                   "chr13\t32336264\t32341196\tBRCA2_ex11\t32.5\nchr13\t32398161\t32400268\tBRCA2_ex27\t31.0\n"
                   "chr17\t43045000\t43045100\tBRCA1_ex23\t20.0\nchr17\t43091435\t43094860\tBRCA1_ex10\t16.2\n", gz=True)
        thr = write(f"{d}/thr.bed.gz", "#chrom\tstart\tend\tregion\t10X\t20X\t30X\n"
                    "chr13\t32336264\t32341196\tBRCA2_ex11\t4932\t4932\t3000\nchr13\t32398161\t32400268\tBRCA2_ex27\t2107\t2107\t1500\n"
                    "chr17\t43045000\t43045100\tBRCA1_ex23\t20\t0\t0\nchr17\t43091435\t43094860\tBRCA1_ex10\t3425\t500\t0\n", gz=True)
        perbase = write(f"{d}/pb.bed.gz", "chr17\t43045000\t43045040\t5\nchr17\t43045040\t43045041\t3\n"
                        "chr17\t43045041\t43045100\t25\nchr13\t32336264\t32341196\t32\n", gz=True)
        median = write(f"{d}/median.txt", "chr13\t31.0\nchr17\t32.0\n")
        versions = write(f"{d}/versions.tsv", "bwa\t0.7.18-r1243\nsamtools\t1.21\ngatk\t4.6.1.0\nmanta\t1.6.0\n")
        commands = write(f"{d}/cmds.sh", "bwa mem ...\n")
        out = {k: f"{d}/out.{k}" for k in ("json", "md", "html")}
        r = self.run_script(f"{BIN}/brca12_report.py", "--sample", "S", "--input-mode", "fastq",
                            "--inputs", "S_R1.fastq.gz,S_R2.fastq.gz", "--reference", "/ref/GRCh38.fa",
                            "--reference-check", f"{d}/check.json", "--exons", f"{d}/exons.norm.bed",
                            "--regions", f"{d}/exact.bed", "--small", f"BRCA1={brca1}", "--small", f"BRCA2={brca2}",
                            "--sv", manta, "--cnv", cns, "--cnr", cnr, "--cov-q20", q20, "--cov-thresholds", thr,
                            "--cov-q0", q0, "--cov-perbase", perbase, "--chrom-median", median,
                            "--versions", versions, "--commands", commands, "--lang", "de",
                            "--out-json", out["json"], "--out-md", out["md"], "--out-html", out["html"])
        self.assertEqual(r.returncode, 0, r.stderr)
        result = json.loads(open(out["json"]).read())
        self.assertEqual(result["schema"], "genome.brca12.v1")
        self.assertIsNone(result["clinicalClassification"])
        b1 = result["genes"]["BRCA1"]
        self.assertEqual(b1["strand"], "-")
        self.assertEqual(len(b1["smallVariants"]), 2)
        frameshift = b1["smallVariants"][0]
        self.assertEqual(frameshift["clinvar"]["stars"], 3)
        self.assertTrue(any(h["code"] == "PVS1 (Kandidat)" for h in frameshift["acmgHints"]))
        missense = b1["smallVariants"][1]
        self.assertTrue(any("DP < 15" in f for f in missense["qualityFlags"]))
        self.assertTrue(any("unsicher abgedecktem Bereich" in f for f in missense["qualityFlags"]))
        self.assertTrue(any(h["code"] == "BS1 (Kandidat)" for h in missense["acmgHints"]))
        b2 = result["genes"]["BRCA2"]["smallVariants"][0]
        self.assertEqual(b2["genotype"], "homozygot")
        self.assertTrue(any(h["code"] == "BA1 (Kandidat)" for h in b2["acmgHints"]))
        ex23 = next(e for e in result["exonCoverage"] if e["label"].startswith("Exon 23"))
        self.assertEqual(ex23["status"], "unsicher")
        self.assertEqual(ex23["label"], "Exon 23 (Legacy 24)")
        self.assertEqual(len(result["cnvs"]), 1)
        self.assertEqual(result["cnvs"][0]["type"], "Deletion")
        self.assertIn("Exon 10 (Legacy 11)", result["cnvs"][0]["affectedExons"])
        self.assertEqual({s["gene"] for s in result["svs"]}, {"BRCA1", "BRCA2"})
        bnd = next(s for s in result["svs"] if s["type"] == "BND")
        self.assertIn("unpräzise Bruchpunkte", bnd["flags"])
        self.assertEqual(len(result["lowCoverageIntervals"]), 1)
        self.assertEqual(result["lowCoverageIntervals"][0]["length"], 41)
        self.assertTrue(result["summary"]["needsAttention"])
        md = open(out["md"]).read()
        html = open(out["html"]).read()
        for text in (md, html):
            self.assertIn("nicht ohne klinisch validierte Bestätigung", text)
            self.assertIn("BRCA1 SNVs/Indels", text)
            self.assertIn("BRCA2 SNVs/Indels", text)
            self.assertIn("CNVs", text)
            self.assertIn("strukturelle Varianten", text)
            self.assertIn("4.6.1.0", text)
            self.assertNotIn("—", text)
        r = self.run_script(*[a if a != "de" else "en" for a in self.last_args])
        self.assertEqual(r.returncode, 0, r.stderr)
        html_en = open(out["html"]).read()
        self.assertIn("must not be interpreted medically without clinically", html_en)
        self.assertIn("BRCA1 SNVs/indels", html_en)
        self.assertIn("Exon 23 (legacy 24)", html_en)
        self.assertIn("heterozygous", html_en)
        self.assertNotIn("unsicher", html_en)


if __name__ == "__main__":
    unittest.main()


class ReviewFixesTest(unittest.TestCase):
    def test_partial_genotypes_are_undetermined(self):
        for gt in ("0/.", "./1", ".|1", "1/."):
            self.assertEqual(genotype_label(gt), "nicht bestimmt", gt)
        self.assertEqual(genotype_label("0/1"), "heterozygot")
        self.assertEqual(genotype_label("1/2"), "compound/multiallelisch")

    def test_manta_breakend_pair_reported_once_from_gene_end(self):
        d = tempfile.mkdtemp()
        header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        rows = ("chr17\t43100000\tMantaBND:1:0\tA\tA]chr2:5000]\t100\tPASS\tSVTYPE=BND;MATEID=MantaBND:1:1\tGT:PR:SR\t0/1:20,5:20,4\n"
                "chr2\t5000\tMantaBND:1:1\tC\tC]chr17:43100000]\t100\tPASS\tSVTYPE=BND;MATEID=MantaBND:1:0\tGT:PR:SR\t0/1:20,5:20,4\n"
                "chr17\t43050000\tMantaBND:2:0\tA\tA[chr17:43060000[\t100\tPASS\tSVTYPE=BND;MATEID=MantaBND:2:1\tGT:PR:SR\t0/1:20,5:20,4\n"
                "chr17\t43060000\tMantaBND:2:1\tA\t]chr17:43050000]A\t100\tPASS\tSVTYPE=BND;MATEID=MantaBND:2:0\tGT:PR:SR\t0/1:20,5:20,4\n")
        path = write(f"{d}/sv.vcf", header + rows)
        regions = [("chr17", 43044294, 43170327, "BRCA1"), ("chr13", 32315076, 32400268, "BRCA2")]
        svs = structural.structural_variants(path, regions, [])
        self.assertEqual([(s["id"], s["chrom"]) for s in svs], [("MantaBND:1:0", "chr17"), ("MantaBND:2:0", "chr17")])

    def test_cnv_targets_equal_access_on_brca_chromosomes(self):
        d = tempfile.mkdtemp()
        check = write(f"{d}/check.json", json.dumps({"chromLengths": {"chr13": 114364328, "chr17": 83257441}}))
        access = write(f"{d}/access.bed", "chr1\t0\t100\nchr13\t10\t200\nchr17\t5\t50\n")
        r = subprocess.run([sys.executable, f"{BIN}/brca12_prepare.py", "--cnv-targets", check,
                            "--cnv-access", access, "--out-cnv-targets", f"{d}/wgs.bed"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        chroms = [line.split("\t")[0] for line in open(f"{d}/wgs.bed")]
        self.assertEqual(chroms, ["chr13", "chr17"])

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

from brca12lib import enigma, evidence, structural, variants  # noqa: E402
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
        self.assertIn("BA1 (Näherung)", codes)

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
        self.assertTrue(any(h["code"] == "PVS1" for h in frameshift["acmgHints"]))
        self.assertTrue(any(h["code"] == "PM5_PTC_Strong" for h in frameshift["acmgHints"]))
        missense = b1["smallVariants"][1]
        self.assertTrue(any("DP < 15" in f for f in missense["qualityFlags"]))
        self.assertTrue(any("unsicher abgedecktem Bereich" in f for f in missense["qualityFlags"]))
        self.assertTrue(any(h["code"] == "BS1 (Näherung)" for h in missense["acmgHints"]))
        b2 = result["genes"]["BRCA2"]["smallVariants"][0]
        self.assertEqual(b2["genotype"], "homozygot")
        self.assertTrue(any(h["code"] == "BA1 (Näherung)" for h in b2["acmgHints"]))
        ex23 = next(e for e in result["exonCoverage"] if e["label"].startswith("Exon 23"))
        self.assertEqual(ex23["status"], "unsicher")
        self.assertEqual(ex23["label"], "Exon 23 (Legacy 24)")
        self.assertEqual(len(result["cnvs"]), 1)
        self.assertEqual(result["cnvs"][0]["type"], "Deletion")
        self.assertIn("Exon 10 (Legacy 11)", result["cnvs"][0]["affectedExons"])
        self.assertTrue(any(c["code"].startswith("PVS1") for c in result["cnvs"][0]["enigmaCodes"]))
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


class EnigmaRulesTest(unittest.TestCase):
    """ClinGen ENIGMA BRCA1/BRCA2 VCEP specifications v1.2 (CSpec GN092, GN097)."""

    def codes(self, *args, **kwargs):
        return [h["code"] for h in enigma.acmg_hints(*args, **kwargs)]

    def gnomad(self, faf, present=True, non_cancer=True):
        best = None if faf is None else {"dataset": "gnomAD exomes", "population": "nfe", "faf": faf,
                                         "nonCancer": non_cancer}
        return {"gnomad": {"datasets": ["gnomAD exomes", "gnomAD genomes"], "present": present, "best": best}}

    def test_frequency_thresholds_use_filtering_allele_frequency(self):
        csq = {"Consequence": "intron_variant"}
        self.assertIn("BA1", self.codes("BRCA1", csq, None, "A", "G", self.gnomad(0.0011)))
        self.assertIn("BS1", self.codes("BRCA1", csq, None, "A", "G", self.gnomad(0.00011)))
        self.assertIn("BS1_Supporting", self.codes("BRCA1", csq, None, "A", "G", self.gnomad(0.00003)))
        # exactly at a threshold the lower code applies ("above", not "at least")
        self.assertIn("BS1", self.codes("BRCA1", csq, None, "A", "G", self.gnomad(0.001)))
        self.assertIn("BS1_Supporting", self.codes("BRCA1", csq, None, "A", "G", self.gnomad(0.0001)))
        low = self.codes("BRCA1", csq, None, "A", "G", self.gnomad(0.00002))
        self.assertFalse([c for c in low if c.startswith(("BA1", "BS1", "PM2"))])

    def test_pm2_supporting_only_for_substitutions_absent_from_gnomad(self):
        csq = {"Consequence": "missense_variant"}
        absent = self.gnomad(None, present=False)
        self.assertIn("PM2_Supporting", self.codes("BRCA1", csq, None, "A", "G", absent))
        indel = self.codes("BRCA1", {"Consequence": "inframe_deletion"}, None, "ATGC", "A", absent)
        self.assertNotIn("PM2_Supporting", indel)
        self.assertIn("PM2_Supporting nicht anwendbar", indel)
        one_file = {"gnomad": {"datasets": ["gnomAD exomes"], "present": False, "best": None}}
        self.assertIn("PM2_Supporting offen", self.codes("BRCA1", csq, None, "A", "G", one_file))
        seen_once = self.gnomad(None, present=True)
        self.assertNotIn("PM2_Supporting", self.codes("BRCA1", csq, None, "A", "G", seen_once))

    def test_without_gnomad_file_hints_are_marked_as_approximation(self):
        codes = self.codes("BRCA1", {"Consequence": "intron_variant"}, 0.00005, "A", "G", {})
        self.assertIn("BS1_Supporting (Näherung)", codes)
        self.assertNotIn("PM2_Supporting (Näherung)",
                         self.codes("BRCA1", {"Consequence": "frameshift_variant"}, None, "AG", "A", {}))

    def test_bayesdel_thresholds_inside_domains(self):
        brct = {"Consequence": "missense_variant", "Protein_position": "1700"}
        self.assertIn("PP3", self.codes("BRCA1", brct, None, "A", "G", {"spliceai": 0.01, "bayesdel": 0.28}))
        self.assertIn("BP4", self.codes("BRCA1", brct, None, "A", "G", {"spliceai": 0.1, "bayesdel": 0.15}))
        between = self.codes("BRCA1", brct, None, "A", "G", {"spliceai": 0.01, "bayesdel": 0.2})
        self.assertFalse({"PP3", "BP4", "BP1_Strong"} & set(between))
        dbd = {"Consequence": "missense_variant", "Protein_position": "2500"}
        self.assertNotIn("PP3", self.codes("BRCA2", dbd, None, "A", "G", {"spliceai": 0.01, "bayesdel": 0.29}))
        self.assertIn("PP3", self.codes("BRCA2", dbd, None, "A", "G", {"spliceai": 0.01, "bayesdel": 0.30}))
        self.assertIn("BP4", self.codes("BRCA2", dbd, None, "A", "G", {"spliceai": 0.01, "bayesdel": 0.18}))
        # benign BayesDel does not give BP4 while splicing is not excluded
        self.assertNotIn("BP4", self.codes("BRCA2", dbd, None, "A", "G", {"spliceai": 0.15, "bayesdel": 0.1}))

    def test_missense_outside_domains_is_bp1_not_bp4(self):
        outside = {"Consequence": "missense_variant", "Protein_position": "1000"}
        codes = self.codes("BRCA1", outside, None, "A", "G", {"spliceai": 0.05, "bayesdel": 0.5})
        self.assertIn("BP1_Strong", codes)
        self.assertFalse({"PP3", "BP4"} & set(codes))

    def test_spliceai_rules(self):
        silent_in = {"Consequence": "synonymous_variant", "Protein_position": "50"}
        self.assertEqual({"BP4", "BP7"} & set(self.codes("BRCA1", silent_in, None, "A", "G", {"spliceai": 0.1})),
                         {"BP4", "BP7"})
        self.assertIn("PP3 (Spleißen)", self.codes("BRCA1", silent_in, None, "A", "G", {"spliceai": 0.2}))
        deep = {"Consequence": "intron_variant", "HGVSc": "ENST:c.212+7A>G"}
        self.assertEqual({"BP4", "BP7"} & set(self.codes("BRCA1", deep, None, "A", "G", {"spliceai": 0.0})),
                         {"BP4", "BP7"})
        near = {"Consequence": "intron_variant&splice_region_variant", "HGVSc": "ENST:c.213-20T>C"}
        codes = self.codes("BRCA1", near, None, "A", "G", {"spliceai": 0.0})
        self.assertIn("BP4", codes)
        self.assertNotIn("BP7", codes)
        self.assertIn("PP3/BP4 offen", self.codes("BRCA1", deep, None, "A", "G", {}))
        canonical = {"Consequence": "splice_donor_variant", "HGVSc": "ENST:c.212+1G>A"}
        self.assertFalse({"PP3 (Spleißen)", "BP4"} & set(self.codes("BRCA1", canonical, None, "G", "A", {"spliceai": 0.9})))

    def test_gnomad_lookup_prefers_non_cancer_and_skips_founder_populations(self):
        d = tempfile.mkdtemp()
        info = ("AC=5;non_cancer_AC=4;faf95_nfe=0.5;non_cancer_faf95_nfe=0.00008;non_cancer_faf95_afr=0.00003;"
                "non_cancer_faf95_asj=0.9;non_cancer_faf95_fin=0.8")
        ex = write(f"{d}/ex.vcf", f"#h\nchr17\t43045050\t.\tC\tT\t.\tPASS\t{info}\n")
        ge = write(f"{d}/ge.vcf", "#h\n17\t43045050\t.\tC\tT\t.\tPASS\tAC_non_cancer=9;faf95_non_cancer_sas=0.0002\n")
        got = evidence.gnomad_lookup([("gnomAD exomes", ex), ("gnomAD genomes", ge)], "chr17", 43045050, "C", "T")
        self.assertTrue(got["present"])
        self.assertEqual((got["best"]["dataset"], got["best"]["population"], got["best"]["faf"]),
                         ("gnomAD genomes", "sas", 0.0002))
        only_exomes = evidence.gnomad_lookup([("gnomAD exomes", ex)], "chr17", 43045050, "C", "T")
        self.assertEqual((only_exomes["best"]["population"], only_exomes["best"]["faf"]), ("nfe", 0.00008))
        absent = evidence.gnomad_lookup([("gnomAD exomes", ex)], "chr17", 43045050, "C", "G")
        self.assertFalse(absent["present"])
        self.assertIsNone(evidence.gnomad_lookup([], "chr17", 1, "A", "C"))

    def test_spliceai_and_bayesdel_lookup(self):
        d = tempfile.mkdtemp()
        sp = write(f"{d}/sp.vcf", "#h\n17\t43045050\t.\tC\tT\t.\t.\tSpliceAI=T|BRCA1|0.01|0.35|0.00|0.02|1|2|3|4,T|NBR2|0.9|0|0|0|1|2|3|4\n")
        bd = write(f"{d}/bd.tsv", "#chr\tpos\tref\talt\tBayesDel_noAF\n17\t43045050\tC\tT\t0.31\n")
        src = evidence.Evidence(spliceai=sp, bayesdel=bd)
        found = src.for_variant("chr17", 43045050, "C", "T", "BRCA1")
        self.assertEqual((found["spliceai"], found["bayesdel"], found["gnomad"]), (0.35, 0.31, None))
        self.assertIsNone(src.for_variant("chr17", 43045050, "C", "A", "BRCA1")["spliceai"])


class InsertionSignalTest(unittest.TestCase):
    """Soft-clip signature of an insertion: clips from both sides, overlapping by the target-site duplication."""

    ALU = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGG"

    @staticmethod
    def sam(name, pos, cigar, seq):
        return "\t".join([name, "0", "chr13", str(pos), "60", cigar, "*", "0", "0", seq, "I" * len(seq)])

    def lines(self, clip):
        tail = "T" * 20 + clip[:20]
        out = [self.sam(f"r{i}", 1000 - i, f"{60 + i}M40S", "A" * (60 + i) + clip[:40]) for i in range(4)]
        out += [self.sam(f"l{i}", 1045, f"40S{60}M", (tail)[-40:] + "C" * 60) for i in range(4)]
        out += [self.sam(f"s{i}", 950 + i, "150M", "G" * 150) for i in range(10)]
        return out

    def run_scan(self, clip, **kw):
        import brca12_insertion_signals as mod
        regions = [("chr13", 0, 5000, "BRCA2")]
        return mod.candidates(*mod.scan(self.lines(clip)), regions, **kw)

    def test_alu_clips_give_a_candidate_with_target_site_duplication(self):
        found = self.run_scan(self.ALU)
        self.assertEqual(len(found), 1)
        self.assertEqual((found[0]["gene"], found[0]["pos"], found[0]["end"]), ("BRCA2", 1045, 1059))
        self.assertEqual(found[0]["targetSiteDuplication"], 15)
        self.assertIn("Alu", found[0]["signature"])
        self.assertEqual(found[0]["spanningReads"], 10)

    def test_clips_without_alu_sequence_are_dropped_by_default(self):
        clip = "ACGTTGCA" * 6
        self.assertEqual(self.run_scan(clip), [])
        self.assertEqual(len(self.run_scan(clip, keep_all=True)), 1)

    def test_report_rows_carry_exon_and_flags(self):
        import json, tempfile
        found = self.run_scan(self.ALU)
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as handle:
            json.dump(found, handle)
        exons = [{"gene": "BRCA2", "chrom": "chr13", "start": 1040, "end": 1100, "label": "Exon 3"}]
        rows = structural.insertion_signals(handle.name, exons)
        self.assertEqual(rows[0]["affectedExons"], ["Exon 3"])
        self.assertEqual(rows[0]["type"], "Insertionssignatur")
        self.assertTrue(rows[0]["flags"][0].startswith("abgeschnittene Leseenden mit Alu-Sequenz"))
        self.assertEqual(structural.insertion_signals(None, exons), [])


class SpikeInFindingsTest(unittest.TestCase):
    """Fixes that came out of the spike-in validation."""

    def test_dosage_hint_for_exon_flagged_uncertain(self):
        exons = [{"gene": "BRCA2", "label": "Exon 15", "dosageRatio": 0.49, "status": "unsicher"},
                 {"gene": "BRCA2", "label": "Exon 11", "dosageRatio": 1.02, "status": "ok"}]
        hints = structural.exon_dosage_hints(exons)
        self.assertEqual([(h["exon"], h["type"]) for h in hints], [("Exon 15", "möglicher Verlust")])
        self.assertIn("unsicher markiert", hints[0]["note"])

    def test_breakend_pair_is_flagged_as_possible_insertion(self):
        header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        rows = ("chr13\t32319152\tb1\tA\t]chr17:43152885]A\t999\tPASS\tSVTYPE=BND;MATEID=m1\tGT:PR:SR\t0/1:53,11:33,12\n"
                "chr13\t32319164\tb2\tA\tA[chr17:43152572[\t999\tPASS\tSVTYPE=BND;MATEID=m2\tGT:PR:SR\t0/1:51,20:32,13\n"
                "chr13\t32350000\tb3\tA\tA[chr5:100[\t999\tPASS\tSVTYPE=BND;MATEID=m3\tGT:PR:SR\t0/1:9,4:3,2\n")
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/sv.vcf"
            with open(path, "w") as handle:
                handle.write(header + rows)
            regions = [("chr13", 32315076, 32400268, "BRCA2")]
            found = structural.structural_variants(path, regions, [])
        flagged = [sv["pos"] for sv in found if any("mögliche Insertion" in f for f in sv["flags"])]
        self.assertEqual(flagged, [32319152, 32319164])

    def test_inserted_alu_sequence_is_flagged(self):
        header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        rows = ("chr13\t32348302\ti1\tT\tTTTGGCCGGGCGCGGTGGCTCACGCCTGTAATCC\t500\tPASS\tEND=32348302;SVTYPE=INS;SVLEN=33\tGT:PR:SR\t1/1:0,28:4,24\n"
                "chr13\t32350000\ti2\tT\tTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t500\tPASS\tEND=32350000;SVTYPE=INS;SVLEN=52\tGT:PR:SR\t0/1:9,4:3,2\n")
        with tempfile.TemporaryDirectory() as d:
            path = f"{d}/sv.vcf"
            with open(path, "w") as handle:
                handle.write(header + rows)
            found = structural.structural_variants(path, [("chr13", 32315076, 32400268, "BRCA2")], [])
        flagged = [sv["pos"] for sv in found if any("Alu-Signatur" in f for f in sv["flags"])]
        self.assertEqual(flagged, [32348302])

"""Targeted assays (exome, panel, amplicon), read QC and mosaic candidates."""
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

import brca12_read_qc  # noqa: E402
import brca12_report  # noqa: E402

STATS = ("SN\traw total sequences:\t1000\nSN\treads mapped:\t980\nSN\treads properly paired:\t960\n"
         "SN\treads duplicated:\t400\nSN\treads MQ0:\t10\nSN\terror rate:\t0.02\nSN\taverage length:\t150\n"
         "IS\t300\t10\t0\t0\nIS\t350\t30\t0\t0\nIS\t400\t20\t0\t0\n")
VCF_HEADER = ("##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: "
              "Allele|Consequence|IMPACT|SYMBOL|Feature|HGVSc|HGVSp|PICK\">\n"
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")


class ReadQcTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def path(self, name, text, gz=False):
        p = os.path.join(self.tmp, name)
        with (gzip.open(p, "wt") if gz else open(p, "w")) as handle:
            handle.write(text)
        return p

    def test_region_stats_and_flags(self):
        stats = brca12_read_qc.read_samtools_stats(self.path("s.txt", STATS))
        self.assertEqual(stats["medianInsertSize"], 350)
        self.assertEqual(stats["duplicatePct"], 40.82)
        out = os.path.join(self.tmp, "qc.json")
        brca12_read_qc.main(["--assay", "exome", "--stats", self.path("s2.txt", STATS), "--out", out])
        flags = json.load(open(out))["flags"]
        self.assertTrue(any(f.startswith("hoher Duplikatanteil") for f in flags))
        self.assertTrue(any(f.startswith("hohe Fehlerrate") for f in flags))
        # amplicon data: duplicates are expected and not flagged
        brca12_read_qc.main(["--assay", "amplicon", "--stats", self.path("s3.txt", STATS), "--out", out])
        self.assertFalse(any(f.startswith("hoher Duplikatanteil") for f in json.load(open(out))["flags"]))

    def test_target_uniformity_fold80_and_weak_targets(self):
        targets = self.path("t.bed", "chr17\t100\t110\tA\nchr17\t200\t210\tB\n")
        depth = "".join(f"chr17\t{p}\t100\n" for p in range(101, 111))
        depth += "".join(f"chr17\t{p}\t{10 if p < 206 else 40}\n" for p in range(201, 211))
        rows, summary = brca12_read_qc.target_coverage(brca12_read_qc.read_targets(targets),
                                                       self.path("d.tsv.gz", depth, gz=True), 30)
        self.assertEqual([r["status"] for r in rows], ["ok", "schwach abgedeckt"])
        self.assertEqual(rows[1]["meanDepth"], 25.0)
        self.assertEqual(summary["meanDepth"], 62.5)
        self.assertEqual(summary["uniformityPct"], 75.0)   # threshold 12.5x: the five 10x bases fail
        self.assertEqual(summary["fold80BasePenalty"], 6.25)
        self.assertEqual(summary["weakTargetCount"], 1)

    def test_primer_clip_stats(self):
        clip = brca12_read_qc.read_clip_stats(self.path("c.txt", "COMMAND: x\nTOTAL READS: 110\nTOTAL CLIPPED: 70\n"
                                                                  "EXCLUDED: 10\nFILTERED: 3\n"))
        self.assertEqual(clip, {"readsConsidered": 100, "readsWithPrimerClipped": 70, "clippedPct": 70.0,
                                "readsFilteredTooShort": 3})


class PrepareTargetsTest(unittest.TestCase):
    def test_targets_inside_genes_and_untargeted_exons(self):
        tmp = tempfile.mkdtemp()
        fai = os.path.join(tmp, "ref.fa.fai")
        open(fai, "w").write("chr13\t114364328\t0\t60\t61\nchr17\t83257441\t0\t60\t61\n")
        exons = os.path.join(ROOT, "assets", "brca12_exons.bed")
        targets = os.path.join(tmp, "targets.bed")
        # numeric contigs, one target outside the genes, one exon of BRCA1 (exon 23) left out
        with open(exons) as handle:
            rows = [l.split("\t") for l in handle if l.strip() and not l.startswith("#")]
        with open(targets, "w") as out:
            out.write("1\t1000\t2000\tother_gene\n")
            for r in rows:
                if r[3] != "BRCA1_ex23":
                    out.write(f"{r[0][3:]}\t{r[1]}\t{r[2]}\t{r[3]}\n")
        args = ["--fai", fai, "--regions", os.path.join(ROOT, "assets", "brca12_regions.bed"), "--exons", exons]
        outs = {k: os.path.join(tmp, k) for k in ("check.json", "calling.bed", "sv.bed", "exact.bed", "exons.bed",
                                                   "targets.brca.bed")}
        subprocess.run([sys.executable, os.path.join(BIN, "brca12_prepare.py"), *args,
                        "--targets", targets, "--out-targets", outs["targets.brca.bed"],
                        "--out-check", outs["check.json"], "--out-calling", outs["calling.bed"],
                        "--out-sv", outs["sv.bed"], "--out-exact", outs["exact.bed"], "--out-exons", outs["exons.bed"]],
                       check=True)
        check = json.load(open(outs["check.json"]))
        self.assertEqual(check["targets"]["intervalsInGenes"], len(rows) - 1)
        self.assertEqual(check["targets"]["exonsWithoutTarget"], ["BRCA1_ex23"])
        first = open(outs["targets.brca.bed"]).readline().split("\t")
        self.assertTrue(first[0].startswith("chr"))

    def test_exon_baseline(self):
        tmp = tempfile.mkdtemp()
        q20 = os.path.join(tmp, "q20.bed.gz")
        with gzip.open(q20, "wt") as h:
            h.write("chr17\t0\t10\tBRCA1_ex01\t100\nchr17\t20\t30\tBRCA1_ex02\t50\nchr13\t0\t10\tBRCA2_ex01\t80\n"
                    "chr13\t20\t30\tBRCA2_ex02\t0\n")
        out = os.path.join(tmp, "median.txt")
        subprocess.run([sys.executable, os.path.join(BIN, "brca12_prepare.py"), "--exon-baseline", q20,
                        "--out-chrom-median", out], check=True)
        self.assertEqual(open(out).read(), "chr13\t80.000\nchr17\t80.000\n")


class MosaicTest(unittest.TestCase):
    def test_candidates_exclude_germline_and_high_fraction(self):
        tmp = tempfile.mkdtemp()
        path = os.path.join(tmp, "m.vcf")
        body = ("chr17\t100\t.\tC\tT\t.\tPASS\t.\tGT:AD:DP:AF\t0/1:90,10:100:0.1\tCSQ=T|missense_variant|MODERATE|BRCA1|T1|c.1C>T|p.X|1\n"
                "chr17\t200\t.\tG\tA\t.\tPASS\t.\tGT:AD:DP:AF\t0/1:50,50:100:0.5\n"
                "chr17\t300\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:AF\t0/1:97,3:100:0.03\n"
                "chr17\t400\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:AF\t0/1:80,20:100:0.2\n")
        # move CSQ into INFO for the first record (column 8)
        lines = []
        for line in body.splitlines():
            f = line.split("\t")
            if len(f) == 11:
                f[7] = f.pop()
            lines.append("\t".join(f))
        open(path, "w").write(VCF_HEADER + "\n".join(lines) + "\n")
        germline = [{"chrom": "chr17", "pos": 400, "ref": "A", "alt": "G"}]
        rows = brca12_report.mosaic_candidates({"BRCA1": path}, germline, [], None, 0.02, 0.30)
        self.assertEqual([r["pos"] for r in rows], [100, 300])
        self.assertEqual(rows[0]["genotype"], "Mosaik-Kandidat")
        # 3 alt reads are the minimum; 4 needed -> only position 100 stays
        rows = brca12_report.mosaic_candidates({"BRCA1": path}, germline, [], None, 0.02, 0.30, min_alt=4)
        self.assertEqual([r["pos"] for r in rows], [100])
        self.assertFalse(any(f.startswith("Allelverhältnis") for f in rows[0]["qualityFlags"]))


if __name__ == "__main__":
    unittest.main()

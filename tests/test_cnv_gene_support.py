"""CNVkit segments reaching beyond a gene count as candidates only if the bins inside the gene agree."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin"))
from brca12lib import structural  # noqa: E402

REGIONS = [("chr17", 43044294, 43170327, "BRCA1")]
EXONS = [{"gene": "BRCA1", "chrom": "chr17", "start": 43045678, "end": 43045802, "label": "Exon 23"},
         {"gene": "BRCA1", "chrom": "chr17", "start": 43124017, "end": 43124115, "label": "Exon 2"}]


def run(segments, bins):
    with tempfile.TemporaryDirectory() as d:
        cns, cnr = os.path.join(d, "s.call.cns"), os.path.join(d, "s.cnr")
        with open(cns, "w") as h:
            h.write("chromosome\tstart\tend\tgene\tlog2\tcn\tprobes\n")
            for s in segments:
                h.write("chr17\t%d\t%d\t-\t%s\t%d\t100\n" % s)
        with open(cnr, "w") as h:
            h.write("chromosome\tstart\tend\tgene\tlog2\n")
            for start, log2 in bins:
                h.write("chr17\t%d\t%d\t-\t%s\n" % (start, start + 1000, log2))
        return structural.cnv_segments(cns, cnr, REGIONS, EXONS)


def gene_bins(log2):
    return [(p, log2) for p in range(43044000, 43171000, 1000)]


class GeneSupport(unittest.TestCase):
    def test_arm_level_segment_with_neutral_gene_bins_is_no_candidate(self):
        # 24 Mb segment at log2 -0.59, bins inside the gene near zero
        (seg,) = run([(22405868, 46087826, -0.593, 1)], gene_bins(-0.10) + [(30000000, -0.9)])
        self.assertEqual(seg["type"], "kein Hinweis")
        self.assertEqual(seg["segmentType"], "Deletion")
        self.assertFalse(seg["supportedInGene"])
        self.assertEqual(seg["flags"], [])
        self.assertAlmostEqual(seg["overlapBinMeanLog2"], -0.10)

    def test_whole_gene_deletion_stays_a_candidate(self):
        # larger segment, bins inside the gene show the loss too
        (seg,) = run([(42000000, 44000000, -0.95, 1)], gene_bins(-0.97))
        self.assertEqual(seg["type"], "Deletion")
        self.assertTrue(seg["supportedInGene"])
        self.assertEqual(seg["affectedExons"], ["Exon 23", "Exon 2"])

    def test_segment_inside_the_gene_is_not_checked(self):
        # a focal segment inside the gene is taken as CNVkit calls it
        (seg,) = run([(43045000, 43047000, -0.26, 1)], [(43045000, -0.2), (43046000, -0.2)])
        self.assertEqual(seg["type"], "Deletion")
        self.assertTrue(seg["supportedInGene"])

    def test_duplication_reaching_beyond_the_gene(self):
        (seg,) = run([(42000000, 44000000, 0.5, 3)], gene_bins(0.55))
        self.assertEqual(seg["type"], "Duplikation")
        (seg,) = run([(42000000, 44000000, 0.3, 3)], gene_bins(0.05))
        self.assertEqual(seg["type"], "kein Hinweis")
        self.assertEqual(seg["segmentType"], "Duplikation")


if __name__ == "__main__":
    unittest.main()

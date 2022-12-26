import unittest
from pathlib import Path

from snprimer.primer import Primer
from snprimer.primer_pair import PrimerPair


class TestPrimerPair(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fasta_file = Path(__file__).parent.resolve() / "ref" / "hg19.fasta"
        f_p = Primer("GGAGATGTACAGCGTGCCATAC", "hg19")
        r_p = Primer("TACATCTTGCTGAGGGGAAGGC", "hg19")
        cls.one_match_primer_pair = PrimerPair(f_p, r_p)
        f_p = Primer("tgacacagtctgcgtttgtaagt", "hg19")
        r_p = Primer("tccaagccatatggtctctctct", "hg19")
        cls.multiple_match_primer_pair = PrimerPair(f_p, r_p)

    def test_ref(self):
        a = Primer("ACGTTGCGATGCTACGATGCTA", "hg38")
        b = Primer("AGTCGCGATGCTACGATGCTA", "hg19")
        with self.assertRaises(Exception):
            PrimerPair(a, b)

    def test_group_position_ranges(self):
        self.assertEqual(1, len(self.one_match_primer_pair._make_mate_pair()))
        self.assertEqual(12, len(self.multiple_match_primer_pair._make_mate_pair()))

    def test_pcr_no_hit(self):
        self.one_match_primer_pair.compatibility = {"max_length": 0}
        sequences = self.one_match_primer_pair.make_pcr(self.fasta_file)
        self.assertEqual(0, len(sequences))
        self.one_match_primer_pair.compatibility = {}

    def test_pcr_1_hit(self):
        hits = self.one_match_primer_pair.make_pcr(self.fasta_file)
        self.assertEqual(
            1,
            len(hits),
        )
        self.assertEqual(
            "GGAGATGTACAGCGTGCCATACAGCACTGGGCAGGGGCAGCCTCAGCAGCAGCAGTTGCCCCCAGCCCAGCCCCAGCCTGCCAGCCAGCAACAAGCTGCCCAGCCTTCCCCTCAGCAAGATGTA",
            hits[0].seq,
        )

    def test_pcr_multiple_hits(self):
        # value come from http://genome.cse.ucsc.edu/
        expected_hits = [
            ("chr3", 197902633, 197902763, "-"),
            ("chr22", 51226056, 51226186, "-"),
            ("chr19", 59099312, 59099442, "-"),
            ("chr16", 90287511, 90287641, "-"),
            ("chr10", 135504547, 135504678, "-"),
            ("chr8", 15328, 15458, "+"),
            ("chr8", 153879, 154009, "+"),
            ("chr2", 114380730, 114380860, "+"),
        ]
        hits = self.multiple_match_primer_pair.make_pcr(self.fasta_file)
        self.assertEqual(12, len(hits))
        for expected_hit in expected_hits:
            self.assertIn(expected_hit[0], [hit.chr for hit in hits])
            self.assertIn(expected_hit[1], [hit.start for hit in hits])
            self.assertIn(expected_hit[2], [hit.end for hit in hits])
            self.assertEqual(
                expected_hit[3],
                [hit.strand for hit in hits if hit.chr == expected_hit[0] and (hit.start == expected_hit[1])][0],
            )

        # strand -
        self.assertEqual(
            "tgacacagtctgcgtttgtaagtaaagttgtaatgggacacagccaatacatgtgttacataatgtctctggctactttcatggtataatggaagagctgagtcattgagagagagaccatatggcttgga",
            hits[0].seq,
        )
        # strand +
        self.assertEqual(
            "tgacacagtctgcgtttgtaagtaaagttgtaataggacacagccaatacatgtgttacataatgtctctggctactttcatggtataatggaagagctgagtcattgagagagagaccatatggcttgga",
            hits[6].seq,
        )

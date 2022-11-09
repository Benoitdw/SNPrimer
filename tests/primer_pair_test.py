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
        f_p = Primer("TGTTCAAACTGATGGGACCCAC", "hg19")
        r_p = Primer("ATTTCTTCATGAAGACCTCACAGT", "hg19")
        cls.multiple_match_primer_pair = PrimerPair(f_p, r_p)

    def test_ref(self):
        a = Primer("ACGTTGCGATGCTACGATGCTA", "hg38")
        b = Primer("AGTCGCGATGCTACGATGCTA", "hg19")
        with self.assertRaises(Exception):
            PrimerPair(a, b)

    def test_group_position_ranges(self):
        self.assertEqual(1, len(self.one_match_primer_pair._make_mate_pair()))
        self.assertEqual(2, len(self.multiple_match_primer_pair._make_mate_pair()))

    def test_pcr_no_hit(self):
        self.one_match_primer_pair.compatibility = {"max_length": 0}
        sequences = self.one_match_primer_pair.make_pcr(self.fasta_file)
        self.assertEqual(0, len(sequences))
        self.one_match_primer_pair.compatibility = {}

    def test_pcr_1_hit(self):
        sequences = self.one_match_primer_pair.make_pcr(self.fasta_file)
        self.assertEqual(
            "GGAGATGTACAGCGTGCCATACAGCACTGGGCAGGGGCAGCCTCAGCAGCAGCAGTTGCCCC"
            "CAGCCCAGCCCCAGCCTGCCAGCCAGCAACAAGCTGCCCAGCCTTCCCCTCAGCAAGATGTA",
            sequences[0].seq,
        )

    def test_pcr_multiple_hits(self):
        sequences = self.multiple_match_primer_pair.make_pcr(self.fasta_file)
        self.assertIn("chr7", [seq.name for seq in sequences])
        self.assertIn("chrX", [seq.name for seq in sequences])

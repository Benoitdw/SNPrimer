import unittest
from pathlib import Path

from snprimer.position_range import PositionRange
from snprimer.primer_maker import PrimerMaker


class TestPrimerMaker(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fasta_file = Path(__file__).parent.resolve() / "ref" / "hg19.fasta"
        cls.seq_args = {}
        cls.global_args = {
            "PRIMER_EXPLAIN_FLAG": 1,
            "PRIMER_MAX_END_STABILITY": 9.0,
            "PRIMER_MAX_LIBRARY_MISPRIMING": 18.00,
            "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING": 50.00,
            "PRIMER_MIN_SIZE": 17,
            "PRIMER_OPT_SIZE": 22,
            "PRIMER_MAX_SIZE": 27,
            "PRIMER_MIN_TM": 58,
            "PRIMER_OPT_TM": 61,
            "PRIMER_MAX_TM": 63,
            "PRIMER_MAX_DIFF_TM": 100,
            "PRIMER_MIN_GC": 20,
            "PRIMER_MAX_GC": 80,
            "PRIMER_SELF_END": 9.00,
            "PRIMER_MAX_POLY_X": 5,
            "PRIMER_GC_CLAMP": 0,
            "PRIMER_PRODUCT_SIZE_RANGE": [[70, 160]],
        }

    def test_make_primers(self):
        target = PositionRange("chr8", 15421, 15422)
        primer_pairs = PrimerMaker(self.fasta_file, self.seq_args, self.global_args, "hg19").design_primers(
            target, n_primers=2
        )
        self.assertIsInstance(primer_pairs, list)
        self.assertIsInstance(primer_pairs[0], dict)
        self.assertIsInstance(primer_pairs[0]["target"], PositionRange)
        self.assertEqual(12, len(primer_pairs[0]["hits"]))
        self.assertIsInstance(primer_pairs[0]["hits"][0], PositionRange)
        self.assertEqual(
            primer_pairs[0]["hits"][0].seq,
            "tgacacagtctgcgtttgtaagtaaagttgtaatgggacacagccaatacatgtgttacataatgtctctggctactttcatggtataatggaagagctgagtcattgagagagagaccatatggcttgga",
        )

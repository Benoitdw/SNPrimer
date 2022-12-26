from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, List

import primer3
from pyfaidx import Fasta

from snprimer.position_range import PositionRange
from snprimer.primer import Primer
from snprimer.primer_pair import PrimerPair

RELATIVE_POSITION = (95, 90)

if TYPE_CHECKING:
    from utils.custom_types import DesignedPrimer


class PrimerMaker:
    def __init__(self, ref_fasta: Path, seq_args: dict, global_args: dict, ref_version: str = "hg19"):
        self.ref_fasta = ref_fasta
        self.seq_args = seq_args
        self.global_args = global_args
        self.ref_version = ref_version

    def design_primers(self, target_range: PositionRange, n_primers: int = 1) -> List[DesignedPrimer]:
        primers_designed = []
        for primer_pair in self.get_primers(target_range, n_primers):
            primers_designed.append(
                {"primers_pair": primer_pair, "target": target_range, "hits": primer_pair.make_pcr(self.ref_fasta)}
            )
        return primers_designed

    def get_primers(self, target_range: PositionRange, n_primers: int = 1) -> list[PrimerPair]:
        sequence = Fasta(self.ref_fasta)[target_range.chr][
            target_range.start - RELATIVE_POSITION[0] : target_range.end + RELATIVE_POSITION[1]
        ]

        seq_args = self.seq_args
        seq_args["SEQUENCE_TEMPLATE"] = sequence.seq
        seq_args["SEQUENCE_TARGET"] = (RELATIVE_POSITION[0], len(target_range))

        global_args = self.global_args
        global_args["PRIMER_NUM_RETURN"] = n_primers

        primers_pair = []
        primers = primer3.bindings.designPrimers(seq_args, global_args)
        for i in range(n_primers):
            right_primer = Primer(primers.get(f"PRIMER_LEFT_{i}_SEQUENCE"), self.ref_version)
            left_primer = Primer(primers.get(f"PRIMER_RIGHT_{i}_SEQUENCE"), self.ref_version)
            primers_pair.append(PrimerPair(right_primer, left_primer))
        return primers_pair


if __name__ == "__main__":
    seq_args = {}
    global_args = {
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
    target_range = PositionRange("chr1", 27100919, 27100919)

    # .make_pcr("/home/benoit/Projects/SNPrimer/tests/ref/hg19.fasta"))

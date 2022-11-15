from pathlib import Path

import primer3
from pyfaidx import Fasta

from snprimer.position_range import PositionRange
from snprimer.primer import Primer
from snprimer.primer_pair import PrimerPair

RELATIVE_POSITION = (95, 90)


class PrimerMaker:
    def __init__(self, ref_fasta: Path, seq_args: dict, global_args: dict, ref_version: str = "hg19"):
        self.ref_fasta = ref_fasta
        self.seq_args = seq_args
        self.global_args = global_args
        self.ref_version = ref_version

    def get_primers(self, position_range: PositionRange, n_primers: int = 1) -> list[PrimerPair]:
        sequence = Fasta(self.ref_fasta)[position_range.chr][
            position_range.start - RELATIVE_POSITION[0] : position_range.end + RELATIVE_POSITION[1]
        ]

        seq_args = self.seq_args
        seq_args["SEQUENCE_TEMPLATE"] = sequence.seq
        seq_args["SEQUENCE_TARGET"] = (RELATIVE_POSITION[0], len(position_range))

        global_args = self.global_args
        global_args["PRIMER_NUM_RETURN"] = n_primers

        primers_pair = []
        primers = primer3.bindings.designPrimers(seq_args, global_args)
        for i in range(n_primers):
            right_primer = Primer(primers.get(f"PRIMER_LEFT_{i}_SEQUENCE"), self.ref_version)
            left_primer = Primer(primers.get(f"PRIMER_RIGHT_{i}_SEQUENCE"), self.ref_version)
            primers_pair.append(PrimerPair(right_primer, left_primer))
        return primers_pair

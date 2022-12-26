from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, List

from dataclasses_json import dataclass_json

from snprimer import PositionRange, Primer

if TYPE_CHECKING:
    pass


@dataclass
@dataclass_json
class PrimerPair:
    forward_primer: Primer
    reverse_primer: Primer
    name: str
    compatibility = dict

    def __init__(self, forward_primer: Primer, reverse_primer: Primer, name: str = None, compatibility: dict = None):
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.name = name
        self.compatibility = compatibility or {}
        self._check_compatibility()

    def _check_compatibility(self):
        if self.forward_primer.ref_version != self.reverse_primer.ref_version:
            raise Exception("Primers in a primer pair should have the same genome of reference.")

    def _make_mate_pair(self) -> list[tuple[PositionRange, PositionRange, str]]:
        mate_pairs = []
        for f_position_range in self.forward_primer.hybridization_place:
            for r_position_range in self.reverse_primer.hybridization_place:
                if self._assert_condition_for_mate_pair(f_position_range, r_position_range):
                    if f_position_range.strand == "-" and r_position_range.strand == "+":
                        mate_pairs.append((r_position_range, f_position_range, "-"))
                    else:
                        mate_pairs.append((f_position_range, r_position_range, "+"))
        return mate_pairs

    def _assert_condition_for_mate_pair(self, f_mate: PositionRange, r_mate: PositionRange) -> bool:
        #  TODO a lot of verification but need to read papers (max len, order -> reverse complement?)
        if f_mate.chr == r_mate.chr:
            return abs(r_mate.start - f_mate.start) < self.compatibility.get("max_length", 150)
        return False

    def make_pcr(self, ref_fasta_file: Path) -> List[PositionRange]:
        hits = []
        for mate_pair in self._make_mate_pair():
            hit = PositionRange(mate_pair[0].chr, mate_pair[0].start, mate_pair[1].end, mate_pair[2])
            hit.set_seq(ref_fasta_file)
            hits.append(hit)
        return hits


if __name__ == "__main__":
    a = Primer("GGAGATGTACAGCGTGCCATAC", "hg19")
    b = Primer("TACATCTTGCTGAGGGGAAGGC", "hg19")
    pp = PrimerPair(a, b)

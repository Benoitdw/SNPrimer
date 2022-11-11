from pathlib import Path

from pyfaidx import Fasta
from pyfaidx import Sequence

from snprimer import PositionRange
from snprimer import Primer


class PrimerPair:
    def __init__(self, forward_primer: Primer, reverse_primer: Primer, name: str = None, compatibility: dict = None):
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.name = name
        self.compatibility = compatibility or {}
        self._check_compatibility()

    def _check_compatibility(self):
        if self.forward_primer.ref_version != self.reverse_primer.ref_version:
            raise Exception("Primers in a primer pair should have the same genome of reference.")

    def _make_mate_pair(self) -> list[tuple[PositionRange, PositionRange]]:
        mate_pairs = []
        for f_position_range in self.forward_primer.position_ranges:
            for r_position_range in self.reverse_primer.position_ranges:
                if self._assert_condition_for_mate_pair(f_position_range, r_position_range):
                    mate_pairs.append((f_position_range, r_position_range))
        return mate_pairs

    def _assert_condition_for_mate_pair(self, f_mate: PositionRange, r_mate: PositionRange) -> bool:
        #  TODO a lot of verification but need to read papers (max len, order -> reverse complement?)
        return (f_mate.chr == r_mate.chr) and (
            r_mate.start - f_mate.end < self.compatibility.get("max_length", 150)
        )  # TODO take from primer design

    def make_pcr(self, ref_fasta_file: Path) -> list[Sequence]:
        hits = []
        fasta_handler = Fasta(ref_fasta_file)
        for mate_pair in self._make_mate_pair():
            hit = fasta_handler[mate_pair[0].chr][mate_pair[0].start - 1 : mate_pair[1].end]
            hit.mate_pair = mate_pair
            hits.append(hit)  # -1 because of difference of index
        return hits


if __name__ == "__main__":
    a = Primer("GGAGATGTACAGCGTGCCATAC", "hg19")
    b = Primer("TACATCTTGCTGAGGGGAAGGC", "hg19")
    pp = PrimerPair(a, b)
    print(pp.make_pcr(Path("/home/benoit/Projects/SNPrimer/tests/ref/hg19.fasta")))
    print(a)
    print(b)

import logging
from dataclasses import dataclass
from dataclasses import field

import gget

from snprimer.position_range import PositionRange


@dataclass
class Primer:
    seq: str
    ref_version: str = "hg38"
    position_ranges: list[PositionRange] = field(default_factory=list, init=False)

    def __post_init__(self):
        if len(self.seq) <= 20:
            logging.warning(f"Primer {self} is too small to apply blast and therefore to get more informations.")
        else:
            self._make_blast(self.ref_version)

    def _make_blast(self, ref_version: str):
        blast = gget.blat(self.seq, seqtype="DNA", assembly=ref_version, json=True)
        if blast:
            for hit in gget.blat(self.seq, seqtype="DNA", assembly=ref_version, json=True):
                self.position_ranges.append(
                    PositionRange(hit["chromosome"], hit["start"], hit["end"], hit["strand"], init_snp=False)
                )


if __name__ == "__main__":
    a = Primer("CACACAGATCAGAGGGCCAAC")
    a.infos(0)

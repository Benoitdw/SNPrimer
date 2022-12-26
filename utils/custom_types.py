from typing import List, TypedDict

from snprimer.position_range import PositionRange
from snprimer.primer_pair import PrimerPair


class DesignedPrimer(TypedDict):
    primers_pair: PrimerPair
    target: PositionRange
    hits: List[PositionRange]

# import logging
from dataclasses import dataclass
from dataclasses import field
from dataclasses import InitVar

import myvariant

from snprimer.snp import SNP


@dataclass
class PositionRange:
    chr: str
    start: int
    end: int
    strand: str = None
    init_snp: InitVar[bool] = True
    snp: list[SNP] = field(default_factory=list)

    def __post_init__(self, init_snp):
        if init_snp:
            self.search_snp()

    def search_snp(self):
        for hit in myvariant.MyVariantInfo().query(f"{self.chr}:{self.start}-{self.end}", fields="dbsnp")["hits"]:
            if (
                "dbsnp" in hit.keys() and "alleles" in hit["dbsnp"]
            ):  # SNP with no alleles informations are not processed
                self.snp.append(SNP(hit["_id"], hit["dbsnp"]["rsid"], hit["dbsnp"]["alleles"], hit["dbsnp"]["ref"]))

    def get_snp(self, max_vaf: float = 0.1) -> list[SNP]:
        if not self.snp:
            self.search_snp()
        return [snp for snp in self.snp if snp.vaf >= max_vaf]

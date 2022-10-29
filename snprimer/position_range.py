import logging
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
        else:
            logging.info(f"{self} is not initialised")

    def search_snp(self) -> list[SNP]:
        for hit in myvariant.MyVariantInfo().query(f"{self.chr}:{self.start}-{self.end}", fields="dbsnp")["hits"]:
            if (
                "dbsnp" in hit.keys() and "alleles" in hit["dbsnp"]
            ):  # SNP with no alleles informations are not processed
                self.snp.append(SNP(hit["_id"], hit["dbsnp"]["rsid"], hit["dbsnp"]["alleles"], hit["dbsnp"]["ref"]))
        return self.snp

    def get_snp(self, max_vaf: float = 0.1) -> list[SNP]:
        return [snp for snp in self.snp if snp.vaf >= max_vaf]


if __name__ == "__main__":
    a = PositionRange("chr8", 19818430, 19818440)
    for snp in a.get_snp(max_vaf=0.05):
        print(snp)

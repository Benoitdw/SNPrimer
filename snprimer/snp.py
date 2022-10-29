from dataclasses import dataclass
from dataclasses import InitVar


@dataclass
class SNP:
    id: str
    rsid: str
    alleles: InitVar[list]
    ref: InitVar[str]
    vaf: float = 0

    def __post_init__(self, alleles, ref):
        for allele in alleles:
            if allele["allele"] == ref:
                self.vaf = self._get_vaf(allele["freq"])

    @staticmethod
    def _get_vaf(allele_freq):
        for key in ["1000g", "exac", "gnomad", "gnomad_exomes", "topmed"]:  # top db
            if allele_freq.get(key):
                return 1 - allele_freq.get(key)
        return 0  # If no match we can estimate the vaf is nul.

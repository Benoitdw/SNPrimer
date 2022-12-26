# import logging
from dataclasses import InitVar, dataclass, field
from pathlib import Path
from typing import Optional

import myvariant
from pyfaidx import Fasta

from snprimer.snp import SNP


@dataclass
class PositionRange:
    chr: str
    start: int
    end: int
    strand: Optional[str] = None
    init_snp: InitVar[bool] = False
    snp: list[SNP] = field(default_factory=list)
    seq: Optional[str] = None

    def __post_init__(self, init_snp):
        if init_snp:
            self.search_snp()

    def __len__(self):
        return self.end - self.start

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

    def set_seq(self, ref_fasta_file: Path):
        fasta_handler = Fasta(ref_fasta_file)
        if self.strand == "-":
            self.seq = fasta_handler[self.chr][self.start - 1 : self.end].reverse.complement.seq
        else:
            self.seq = fasta_handler[self.chr][self.start - 1 : self.end].seq

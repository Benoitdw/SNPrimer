import logging

import click

from snprimer import Primer
from snprimer import PrimerPair


@click.group()
def cli():
    pass


@click.command()
def check():
    print("ok")


@click.command()
@click.option("--snp/--no-snp", default=False, help="Check for the presence of SNP in primer.")
@click.option("--ref_fasta_file", type=click.Path(exists=True))  # TODO make mandatory
@click.option("--ref", type=click.Choice(["hg19", "hg38"], case_sensitive=False), default="hg38")
@click.argument("primer_file", type=click.File("r"))
def pcr(snp, ref_fasta_file, ref, primer_file):
    logging.info("Start In silico PCR...")
    for line in primer_file.readlines():
        l_p, r_p = [Primer(a.strip(), "hg19") for a in line.split("\t")]
        primer_pair = PrimerPair(l_p, r_p)
        pcr_hits = primer_pair.make_pcr(ref_fasta_file)
        if len(pcr_hits) > 1:
            logging.warning("Multiple pcr hit for primer pair")
        for pcr_hit in pcr_hits:
            logging.info(pcr_hit.seq)


cli.add_command(check)
cli.add_command(pcr)

if __name__ == "__main__":
    cli()


# SNPrimer check primer.tsv --snp

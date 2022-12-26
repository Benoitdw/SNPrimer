import click

from snprimer import PositionRange


@click.group()
def cli():
    pass


@click.command()
@click.option("--vaf", default=0.1, type=float, help="minimum snp vaf to report")
@click.argument("position_file", type=click.File("r"))
@click.argument("output", type=click.File("w"), default="-")
def check(vaf, position_file, output):
    header = ["chr", "start", "end", "found_snp"]
    header_stringify = "\t".join(header)
    output.write(f"{header_stringify}\n")
    for row in position_file.readlines():
        row = [c.strip() for c in row.split("\t")]
        position_range = PositionRange(row[0], row[1], row[2])
        output_line = [position_range.chr, position_range.start, position_range.end]
        output_line.append(str(";".join([snp.rsid for snp in position_range.get_snp(vaf)])))
        output_line += "\n"
        output.write("\t".join(output_line))


cli.add_command(check)

if __name__ == "__main__":
    cli()


# SNPrimer check primer.tsv --snp

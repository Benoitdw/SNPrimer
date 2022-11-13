# SNPrimer

Small python library to search snp in primer by position or by sequence.

## Installation

`pip install snprimer`

## Usage

### 🎯 By position 🎯

```python
from snprimer import PositionRange

position = PositionRange("chr8",19818430, 19818440)
print(position)
# PositionRange(chr='chr8', start=19818430, end=19818440, strand=None, snp=[SNP(id='chr8:g.19818436C>G', rsid='rs316', vaf=0.15300000000000002), SNP(id='chr8:g.19818436C>T', rsid='rs316', vaf=0.15300000000000002)])
for snp in position.get_snp(max_vaf=0.05):
    print(snp)
#SNP(id='chr8:g.19818436C>G', rsid='rs316', vaf=0.15300000000000002)
#SNP(id='chr8:g.19818436C>T', rsid='rs316', vaf=0.15300000000000002)
```

### 🔤 By sequence 🔤
```python
from snprimer import Primer

primer = Primer("CACACAGATCAGAGGGCCAAC")
print(primer)
#Primer(seq='CACACAGATCAGAGGGCCAAC', position_ranges=[PositionRange(chr='chr1', start=26774827, end=26774847, strand='+', snp=[SNP(id='chr1:g.26774827G>A', rsid='rs2075289787', vaf=0), SNP(id='chr1:g.26774830A>G', rsid='rs986550282', vaf=0.0), SNP(id='chr1:g.26774842T>C', rsid='rs1440652363', vaf=0.0)])])
```

### 🖥️In silico PCR🖥️

```python
from pathlib import Path
from snprimer import Primer, PrimerPair

a = Primer("GGAGATGTACAGCGTGCCATAC", "hg19")
b = Primer("TACATCTTGCTGAGGGGAAGGC", "hg19")
pp = PrimerPair(a, b)
pcr = pp.make_pcr(Path('/path/to/hg19.fa'))
print(pcr)
```



# CLI

You can also use SNPrimer in cli with two options :

* check will check the presence of snp in a position range
* pcr will make the insilico pcr

## Installation

```bash
git clone git@github.com:Benoitdw/SNPrimer.git
cd SNPrimer
python -m venv venv # recommended work in a virtual environment
source venv/bin/activate # Recommended activate the virtual environment


poetry install
# OR pip install .
```



## Check

```bash
Usage: main.py check [OPTIONS] POSITION_FILE [OUTPUT]

Options:
  --vaf FLOAT  minimum snp vaf to report
```

`python main.py check --vaf 0.00 tests/files/position.tsv tests/files/output_position.tsv`

will produce the `output_position.tsv` file you can see in `tests/files` directory



## PCR

```bash
Usage: main.py pcr [OPTIONS] PRIMER_FILE [OUTPUT]

Options:
  --snp / --no-snp       Check for the presence of SNP in primer.
  --vaf FLOAT            minimum snp vaf to report
  --ref_fasta_file PATH  Path to the reference fasta file.  [required]
  --ref [hg19|hg38]
  --help                 Show this message and exit.

```

`main.py pcr --snp --vaf 0 --ref_fasta_file tests/ref/hg19.fasta --ref hg19 tests/files/test_primer.tsv tests/files/output.tsv `

will produce the `output.tsv` you can see in `tests/files` directory

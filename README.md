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

### Design
```python
from snprimer import PrimerMaker, PositionRange
from pathlib import Path
seq_args = {
        }
global_args = {
        "PRIMER_EXPLAIN_FLAG": 1,
        "PRIMER_MAX_END_STABILITY": 9.0,
        "PRIMER_MAX_LIBRARY_MISPRIMING": 18.00,
        "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING": 50.00,
        "PRIMER_MIN_SIZE": 17,
        "PRIMER_OPT_SIZE": 22,
        "PRIMER_MAX_SIZE": 27,
        "PRIMER_MIN_TM": 58,
        "PRIMER_OPT_TM": 61,
        "PRIMER_MAX_TM": 63,
        "PRIMER_MAX_DIFF_TM": 100,
        "PRIMER_MIN_GC": 20,
        "PRIMER_MAX_GC": 80,
        "PRIMER_SELF_END": 9.00,
        "PRIMER_MAX_POLY_X": 5,
        "PRIMER_GC_CLAMP": 0,
        "PRIMER_PRODUCT_SIZE_RANGE": [[70, 160]],
    }
position_range = PositionRange("chr1",27100919,27100919)
ref_fasta = Path("/path/to/fasta")
print(PrimerMaker(ref_fasta,
                  seq_args, global_args).get_primers(position_range)[0]
      .make_pcr(ref_fasta))
#[>chr1:27100876-27100999
#GGAGATGTACAGCGTGCCATACAGCACTGGGCAGGGGCAGCCTCAGCAGCAGCAGTTGCCCCCAGCCCAGCCCCAGCCTGCCAGCCAGCAACAAGCTGCCCAGCCTTCCCCTCAGCAAGATGTA]
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

## 🎯 Check 🎯

```bash
Usage: main.py check [OPTIONS] POSITION_FILE [OUTPUT]

Options:
  --vaf FLOAT  minimum snp vaf to report
```

`python main.py check --vaf 0.00 tests/files/position.tsv tests/files/output_position.tsv`

will produce the `output_position.tsv` file you can see in `tests/files` directory

## 🖥️ PCR 🖥️

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

## 🔎 DESIGN 🔎
This utilities aims to produce a list of primer for a range of position.

```bash
Usage: main.py design [OPTIONS] BED_FILE [OUTPUT]

Options:
  --snp / --no-snp       Check for the presence of SNP in primer.
  --vaf FLOAT            minimum snp vaf to report
  --n_primers INTEGER    number of primers to generate by position
  --ref_fasta_file PATH  Path to the reference fasta file.  [required]
  --ref [hg19|hg38]
  --config FILENAME
  --help                 Show this message and exit.
```
An example of config file can be find in tests/files/config.json. This file stores the configuration needed by primer3 to design primer.
If no config file is given. default config is the following:

```json
{
                "PRIMER_EXPLAIN_FLAG": 1,
                "PRIMER_MAX_END_STABILITY": 9.0,
                "PRIMER_MAX_LIBRARY_MISPRIMING": 18.00,
                "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING": 50.00,
                "PRIMER_MIN_SIZE": 17,
                "PRIMER_OPT_SIZE": 22,
                "PRIMER_MAX_SIZE": 27,
                "PRIMER_MIN_TM": 58,
                "PRIMER_OPT_TM": 61,
                "PRIMER_MAX_TM": 63,
                "PRIMER_MAX_DIFF_TM": 100,
                "PRIMER_MIN_GC": 20,
                "PRIMER_MAX_GC": 80,
                "PRIMER_SELF_END": 9.00,
                "PRIMER_MAX_POLY_X": 5,
                "PRIMER_GC_CLAMP": 0,
                "PRIMER_PRODUCT_SIZE_RANGE": [[70, 160]]
}
```

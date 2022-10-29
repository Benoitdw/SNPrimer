# SNPrimer

Small python library to search snp in primer by position or by sequence.

## Installation

`pip install snprimer``

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
primer.infos(max_vaf=0)
#CACACAGATCAGAGGGCCAAC has snp with vaf > 0 : [SNP(id='chr1:g.26774827G>A', rsid='rs2075289787', vaf=0), SNP(id='chr1:g.26774830A>G', rsid='rs986550282', vaf=0.0), SNP(id='chr1:g.26774842T>C', rsid='rs1440652363', vaf=0.0)]
```

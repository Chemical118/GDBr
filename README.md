GDBr : Genome identification tool for Double-strand Break Repair
================
<img src="logo/gdbr.svg" alt="GDBr logo" align="right" height="160" style="display: inline-block;"> GDBr (pronounced _GDB; Genome Debugger_) is tool that identify Double-strand Break Repair (DSBR) using genome and variant. GDBr goes through three processes to identify DSBR. First step is preprocess the genome using [`RagTag`](https://github.com/malonge/RagTag) and [`svim-asm`](https://github.com/eldariont/svim-asm) and make sure they have same chromosome name with reference. Second step is correct the variant using [`BLAST`](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and save a csv file. Last step is to segregate the corrected variants into the appropriate DSBRs. 

[![CI](https://github.com/Chemical118/GDBr/workflows/CI/badge.svg)](https://github.com/Chemical118/GDBr/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/Chemical118/GDBr/branch/master/graph/badge.svg?token=NA5V5H52M6)](https://codecov.io/gh/Chemical118/GDBr)
[![anaconda](https://anaconda.org/chemical118/gdbr/badges/version.svg)](https://anaconda.org/Chemical118/gdbr)

You need only reference sequence and query sequences file to use `GDBr`.

### Install
We strongly recommend using `conda` package manager to install `GDBr`
```sh
conda install -c bioconda -c chemical118 -c conda-forge gdbr
```

### Preprocess
By using `RagTag` and `svim-asm`, `GDBr` preprocess data and return properly scaffolded query `.fa` sequence file and variant `.vcf` file.
```sh
gdbr preprocess -r <reference.fa> -q <query1.fa query2.fa ...> -o prepro -t <threads>
```
> Preprocess step needs lots of memory, turn on `--low_memory` if you run out of memory
### Correct
By using `BLAST`, `GDBr` correct the variant file to analysis DSBR accurately.
```sh
gdbr correct -r <reference.fa> -q prepro/*.PRE.fa -v prepro/*.PRE.vcf -o sv -t <threads>
```

### Analysis
`GDBr` analysis the variant and identify DSBR mechanism.
```sh
gdbr analysis -r <reference.fa> -q prepro/*.PRE.fa -v sv/*.COR.csv -o dsbr -t <threads>
```

#### Final output
`GDBr`'s final ouput is `<query basename>.ANL.tsv`. This is simple description of the final output.

| Field             | Description                                          |
|-------------------|------------------------------------------------------|
| ID                | GDBr.\<query order\>.\<variant order\>               |
| SV_TYPE           | INS, DEL, SUB, etc (fail to find variant)            |
| CHR               | variant chromosome                                   |
| REF_START         | variant reference start location                     |
| REF_END           | variant reference end location                       |
| QRY_START         | variant query start location                         |
| QRY_END           | variant query end location                           |
| REPAIR_TYPE       | DSBR mechanism type                                  |
| HOM_LEN/HOM_START | INDEL : homology length / SUB : left homology length |
| HOM_END           | right homology length                                |
| DSBR_CHR          | different locus DSBR chromosome                      |
| DSBR_START        | different locus DSBR start                           |
| DSBR_END          | different locus DSBR end                             |
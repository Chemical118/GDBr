GDBr : Genome identification tool for Double-strand Break Repair
================================================================

<img src="logo/gdbr.svg" alt="GDBr logo" align="right" height="160" style="display: inline-block;"> GDBr (pronounced _Genome Debugger_) is tool that identify Double-strand Break Repair (DSBR) using genome and variant. GDBr goes through three processes to identify DSBR. First step is preprocess the genome using [`RagTag`](https://github.com/malonge/RagTag) and [`svim-asm`](https://github.com/eldariont/svim-asm) and make sure they have same chromosome name with reference. Second step is correcting the variant using [`BLAST`](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and filtering the variant which have repeat bt [`TRF`](https://github.com/Benson-Genomics-Lab/TRF) and [`RepeatMasker`](https://github.com/rmhubley/RepeatMasker), then save a csv file. Last step is to segregate the corrected variants into the appropriate DSBRs.

[![CI](https://github.com/Chemical118/GDBr/workflows/CI/badge.svg)](https://github.com/Chemical118/GDBr/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/Chemical118/GDBr/branch/master/graph/badge.svg?token=NA5V5H52M6)](https://codecov.io/gh/Chemical118/GDBr)
[![anaconda](https://anaconda.org/chemical118/gdbr/badges/version.svg)](https://anaconda.org/Chemical118/gdbr)

You need only reference sequence and query sequences file to use `GDBr`.

### Install

We strongly recommend using `conda` package manager to install `GDBr`.

```sh
conda create -n GDBr -c conda-forge -c bioconda -c chemical118 gdbr
conda activate GDBr
gdbr --version
```
Also, you can use [`mamba`](https://github.com/conda-forge/miniforge) package mamager to install `GDBr` quickly.
```sh
mamba create -n GDBr -c conda-forge -c bioconda -c chemical118 gdbr
mamba activate GDBr
gdbr --version
```
### Quick Start
We require references and pangenomes made from long-read sequencing for accurate results. It is recommended that reference is assembled at the chromosome-level, and query should be assembled at the scaffold-level.

```sh
gdbr analysis -r <reference.fa> -q <query1.fa query2.fa ...> -s <species of data> -t <number of threads>
```

### Example

```sh
mkdir gdbr_test
cd gdbr_test

# download reference
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

# download pangenome
wget s3://human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_f1_assembly_v2_genbank/HG002.paternal.f1_assembly_v2_genbank.fa.gz
wget s3://human-pangenomics/working/HPRC_PLUS/HG005/assemblies/year1_f1_assembly_v2_genbank/HG005.paternal.f1_assembly_v2_genbank.fa.gz

# decompress genome
gzip -d chm13v2.0.fa.gz
gzip -d HG002.paternal.f1_assembly_v2_genbank.fa.gz
gzip -d HG005.paternal.f1_assembly_v2_genbank.fa.gz

# install GDBr
conda create -n GDBr -c conda-forge -c bioconda -c chemical118 gdbr
conda activate GDBr

# run GDBr
gdbr analysis -r chm13v2.0.fa -q HG002.paternal.f1_assembly_v2_genbank.fa HG005.paternal.f1_assembly_v2_genbank.fa -s human -o gdbr_output -t 10
```

### Steps of GDBr

The above command executes the following three processes simultaneously. If you want to redo some of the processes, you can manually run the command below.

#### Preprocess

By using `RagTag` and `svim-asm`, `GDBr` preprocess data and return properly scaffolded query `.fa` sequence file and variant `.vcf` file.

```sh
gdbr preprocess -r <reference.fa> -q <query1.fa query2.fa ...> -o prepro -t <number of threads>
```

The preprocess step requires the use of a sorting program to do the scaffolding and variant calling, and even though they are allocated a lot of threads, they still don't use them all. An optimization approach is to distribute multiple queries to a small number of threads. However, this approach requires very high memory usage, so GDBr was developed to allow the user to freely choose this optimization by providing the `--low_memory` option.

#### Correct

By using `BLAST`, `GDBr` correct the variant file to analysis DSBR accurately. And, filter the repeat by using `TRF`, `RepeatMasker`.

```sh
gdbr correct -r <reference.fa> -q prepro/query/*.GDBr.preprocess.fa -v prepro/vcf/*.GDBr.preprocess.vcf -s <species of data> -o sv -t <number of threads>
```

#### Analysis

`GDBr` analysis the variant and identify DSBR mechanism.

```sh
gdbr analysis -r <reference.fa> -q prepro/query/*.GDBr.preprocess.fa -v sv/*.GDBr.correct.csv -o dsbr -t <number of threads>
```

You can turn on different locus DSBR analysis by `--diff_locus_dsbr_analysis`, however analysis can give false positives due to partial homology on the sex chromosomes.

### Final output

`GDBr`'s final ouput is `<query basename>.GDBr.result.tsv`. This is simple description of the final output.

| Field                 | Description                                              |
| --------------------- | -------------------------------------------------------- |
| ID                    | GDBr.\<query order\>.\<variant order\>                   |
| CALL_TYPE             | Variant type : INS, DEL, etc                             |
| SV_TYPE               | Corrected variant type : INS, DEL, SUB, etc              |
| CHR                   | variant chromosome                                       |
| REF_START             | variant reference start location                         |
| REF_END               | variant reference end location                           |
| QRY_START             | variant query start location                             |
| QRY_END               | variant query end location                               |
| GDBR_TYPE             | GDBr variant type                                        |
| HOM_LEN/HOM_START_LEN | INDEL : homology length / SUB : left homology length     |
| HOM_END_LEN           | SUB : right homology length                              |
| TEMP_INS_SEQ_LOC      | templated insertion sequence location (REF or QRY)       |
| DSBR_CHR              | different locus DSBR chromosome                          |
| DSBR_START            | different locus DSBR start                               |
| DSBR_END              | different locus DSBR end                                 |
| HOM_SEQ/HOM_START_SEQ | INDEL : homology sequence / SUB : left homology sequence |
| HOM_END_SEQ           | SUB : right homology sequence                            |
| PUTATIVE_MECHANISM    | GDBr DSB repair putative mechanism                       |

### Benckmarking

You can benchmark any command in GDBr with the `--benchmark` option by `GNU time` and `psutil`. It provides user time, system time, average CPU usage, multiprocessing efficiency, maximum RAM usage and wall clock time.

```
...
[2023-08-18 13:44:16] GDBr benchmark complete
User time (seconds) : 8007.44
System time (seconds) : 19901.00
Percent of CPU this job got : 7267%
Multiprocessing efficiency : 0.5118
Wall clock time (h:mm:ss or m:ss) : 6:23.99
Max memory usage (GB) : 23.5805
```

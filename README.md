GDBr : Genome identification tool for Double-strand Break Repair
================
<img src="logo/gdbr.svg" alt="GDBr logo" align="right" height="160" style="display: inline-block;"> GDBr (pronounced _GDB; Genome Debugger_) is tool that identify Double-strand Break Repair (DSBR) using genome and variant. GDBr goes through two processes to identify DSBR. The first step is to blastn the variants and save the corrected variants to a csv file. The second process is to segregate the corrected variants into the appropriate DSBRs. 

[![CI](https://github.com/Chemical118/GDBr/workflows/CI/badge.svg)](https://github.com/Chemical118/GDBr/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/Chemical118/GDBr/branch/master/graph/badge.svg?token=NA5V5H52M6)](https://codecov.io/gh/Chemical118/GDBr)

You need reference and query sequence file and variant calling by [SVIM-asm](https://github.com/eldariont/svim-asm). Also, GDBr human pangenome preprocess pipeline is provided in `example/` folder.
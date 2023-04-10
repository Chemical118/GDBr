GDBr : Genome identification tool for Double strand Break Repair
================
<img src="logo/gdbr.svg" alt="GDBr logo" align="right" height="160" style="display: inline-block;"> GDBr (pronounced _GDB; Genome Debugger_) is tool that identify Double Strand Break Repair (DSBR) using genome and variant. GDBr goes through two processes to identify DSBR. The first step is to blastn the variants and save the corrected variants to a csv file. The second process is to segregate the corrected variants into the appropriate DSBRs. 

You need reference and query sequence file and variant calling by [SVIM-asm](https://github.com/eldariont/svim-asm). Also, GDBr human pangenome preprocess pipeline is provided in `example/` folder.
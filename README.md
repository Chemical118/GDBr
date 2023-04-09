GDBr : Genome identification tool for Double strand Break Repair
================
GDBr (pronounced _GDB_) is tool that identify DSBs using genome and variant. GDBr goes through two processes to identify Double Strand Break Repair. The first step is to blastn the variants and save the corrected variants to a csv file. The second process is to segregate the corrected variants into the appropriate DSBRs.

You need reference and query sequence file and variant calling by [SVIM-asm](https://github.com/eldariont/svim-asm). Also, GDBr human pangenome preprocess pipeline is provided in `example/` folder.
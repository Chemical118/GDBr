ref=chm13v2.0.fa
qry_save=qryseq
var_save=vcfs
var_min_size=50
core=16

mkdir $qry_save
mkdir $var_save

mkdir svim_asm
mkdir ragtag

for pan in `ls *.f1_assembly_v2_genbank.fa.gz`; do
name=$(basename -s .f1_assembly_v2_genbank.fa.gz $pan)

echo "<$name genome preprocess>"

gunzip -dk $pan
panfa=$(basename -s .gz $pan)

ragtag.py scaffold -t $core -u -C -o "ragtag/$name" $ref $panfa

pre_panfa="$qry_save/$name.fa"
faidx -e "lambda t: t.replace('_RagTag', '')" -g "[^Chr0]" "ragtag/$name/ragtag.scaffold.fasta" -o $pre_panfa

svim_asm_save="svim_asm/$name"
mkdir $svim_asm_save

minimap2 -a -x asm5 --cs -r2k -t $core $ref $pre_panfa > "$svim_asm_save/$name.sam"
samtools sort -m4G -@ 8 -o "$svim_asm_save/$name.sorted.bam" "$svim_asm_save/$name.sam"
samtools index "$svim_asm_save/$name.sorted.bam"
svim-asm haploid $svim_asm_save "$svim_asm_save/$name.sorted.bam" $ref --min_sv_size $var_min_size --tandem_duplications_as_insertions --interspersed_duplications_as_insertions

cp "$svim_asm_save/variants.vcf" "$var_save/$name.vcf"
done
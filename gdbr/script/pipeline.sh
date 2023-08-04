pan=$1
ref=$2
qry_save=$3
var_save=$4
ragtag_save=$5
svim_asm_save=$6
var_min_size=$7
core=$8
trust=$9

name=$(basename $pan)

echo "<$name genome preprocess>"

if [ "$trust" -eq 0 ]; then
    ragtag.py scaffold -t $core -u -C -o $ragtag_save $ref $pan
fi

pre_panfa="$qry_save/$name.GDBr.preprocess.fa"
faidx -e "lambda t: t.replace('_RagTag', '')" -g "[^Chr0]" "$ragtag_save/ragtag.scaffold.fasta" -o $pre_panfa


if [ "$trust" -eq 0 ]; then
    minimap2 -a -x asm5 --cs -r2k -t $core $ref $pre_panfa > "$svim_asm_save/$name.sam"
    samtools sort -m4G -@ $core -o "$svim_asm_save/$name.sorted.bam" "$svim_asm_save/$name.sam"
    rm "$svim_asm_save/$name.sam"
    samtools index "$svim_asm_save/$name.sorted.bam"
fi

svim-asm haploid $svim_asm_save "$svim_asm_save/$name.sorted.bam" $ref --min_sv_size $var_min_size --tandem_duplications_as_insertions --interspersed_duplications_as_insertions

cp "$svim_asm_save/variants.vcf" "$var_save/$name.GDBr.preprocess.vcf"
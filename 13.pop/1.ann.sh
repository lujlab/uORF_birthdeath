#uORF annotation with UTRannotator
vep -i mel1356.vcf --gtf Drosophila_melanogaster.BDGP6.32.105.sort.gtf.gz --fasta ~/lcl/dataset/dm6/dm6.24/dmel-all-chromosome-r6.24.fasta --plugin UTRannotator -o mel1356_vep.out

cat mel1356_vep.out |awk '$7~/5_prime_UTR_variant/' |grep -v "five_prime_UTR_variant_annotation=-;five_prime_UTR_variant_consequence=-" > 1356_uORF_var.txt
sed 's/;/\t/g' 1356_uORF_var.txt|cut -f 1-5,8,15,21|sed 's/five_prime_UTR_variant_consequence=//g' >1356_uORF_var.txt2
nohup python uORF_var_type_filter.py 1356_uORF_var.txt2 1356_uORF_var_type.bed &

bedtools intersect -a 1356_uORF_var_type.bed -b freq/all_pop_freq.bed -wa -wb >1356_uORF_var_type_freq.bed

cut -f4 1356_uORF_var_type_freq.bed|sed 's/_/\t/g'|sed 's$/$\t$g' |paste - 1356_uORF_var_type_freq.bed|awk -vOFS='\t' '{if ($3!="-"&&$4!="-" &&$4==$22) print $0;if ($3=="-"&&length($22)==length($4)+1) print $0;if ($4=="-"&&length($21)==length($3)+1) print $0 }'>1356_uORF_var_type_freq2.bed

cut -f 5-8,11,15-46 1356_uORF_var_type_freq2.bed |sort -k1,1V -k2,2n -k3,3n |uniq >1356_uORF_var_type_freq3.bed
awk '$11/(2*1356)>=0.9' 1356_uORF_var_type_freq3.bed >1356_uORF_var_type_freq_missing0.9.bed
awk -F'\t' '{split($4,a,"/");split(a[1],b,"_");if (b[3]=="-"&&substr($13,2)==a[2]) print $0; else if (a[2]=="-") print $0;else if (a[2]!="-" && b[3]!="-" &&a[2]==$13) print $0}' 1356_uORF_var_type_freq_missing0.9.bed >1356_uORF_var_type_freq_missing0.9.bed2

nohup bedtools intersect -a 1356_uORF_var_type_freq_missing0.9.bed -b SNP.bed >1356_uORF_var_type_SNPonly.bed & #high quality filtered SNPs (SNP.bed) were derived from the previous study (Chen et al, 2024)
awk 'length($12)>1||length($13)>1' 1356_uORF_var_type_freq_missing0.9.bed2 > 1356_uORF_var_type_INDELonly.bed

#remove regions overlapped with CDS
bedtools intersect -a 1356_uORF_var_type_SNPonly.bed -b /gpfs2/liucl/lcl/dataset/annotation_file/gtf/dmel-all-r6.24.gtf.cds.bed -v >1356_uORF_var_type_SNPonly_rmCDS.bed
bedtools intersect -a 1356_uORF_var_type_INDELonly.bed -b /gpfs2/liucl/lcl/dataset/annotation_file/gtf/dmel-all-r6.24.gtf.cds.bed -v >1356_uORF_var_type_INDELonly_rmCDS.bed




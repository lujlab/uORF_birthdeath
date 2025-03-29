bcftools view -i 'ID=@id_retained_compound' -o compound_uORF.vcf.gz -Oz merged_uORF_rmdup_fixploidy.vcf.gz
bcftools view -i 'ID=@id_retained_uAUG_gained' -o uAUG_gained_uORF.vcf.gz -Oz merged_uORF_rmdup_fixploidy.vcf.gz
bcftools view -i 'ID=@id_retained_uAUG_lost' -o uAUG_lost_uORF.vcf.gz -Oz merged_uORF_rmdup_fixploidy.vcf.gz
bcftools view -i 'ID=@id_retained_uFrameShift' -o uFrameShift_uORF.vcf.gz -Oz merged_uORF_rmdup_fixploidy.vcf.gz
bcftools view -i 'ID=@id_retained_uSTOP_gained' -o uSTOP_gained_uORF.vcf.gz -Oz merged_uORF_rmdup_fixploidy.vcf.gz
bcftools view -i 'ID=@id_retained_uSTOP_lost' -o uSTOP_lost_uORF.vcf.gz -Oz merged_uORF_rmdup_fixploidy.vcf.gz

#tajimasD
nohup vcftools --gzvcf compound_uORF.vcf.gz --TajimaD 50000000 --out compound_uORF &
nohup vcftools --gzvcf uAUG_gained_uORF.vcf.gz --TajimaD 50000000 --out uAUG_gained_uORF &
nohup vcftools --gzvcf uAUG_lost_uORF.vcf.gz --TajimaD 50000000 --out uAUG_lost_uORF &
nohup vcftools --gzvcf uFrameShift_uORF.vcf.gz --TajimaD 50000000 --out uFrameShift_uORF &
nohup vcftools --gzvcf uSTOP_gained_uORF.vcf.gz --TajimaD 50000000 --out uSTOP_gained_uORF &
nohup vcftools --gzvcf uSTOP_lost_uORF.vcf.gz --TajimaD 50000000 --out uSTOP_lost_uORF &

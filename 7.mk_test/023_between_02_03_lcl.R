#### Extract polymorphic differences in different regions (1356samples) ####
# only include snps
setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/9.compensation/MK_test/pop_AF")
#for i in `ls /gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/1.vcf_snp_and_indel_new/*_snp_and_indel_split_decompose.vcf`;do j=`basename $i _snp_and_indel_split_decompose.vcf`; echo "vcftools --vcf $i --remove-indels --min-alleles 2 --max-alleles 2 --freq --max-missing 0.9 --out ${j} ">>run1.sh;done
#ParaFly -c run1.sh -CPU 30
#cat *frq|grep -v "CHROM" |sed 's/:/\t/g' |awk -vOFS='\t' '{print "chr"$1,$2-1,$2,$1"_"$2"_SNP",".","+"}' |sort -k1,1V -k2,2n >1356_snp_dm6.bed
#for i in `ls /gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/1.vcf_snp_and_indel_new/*_snp_and_indel_split_decompose.vcf`;do j=`basename $i _snp_and_indel_split_decompose.vcf`; echo "vcftools --vcf $i --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.9 --counts --out ${j} ">>run2.sh;done
#cat *count |grep -v "CHR"|sed 's/:/\t/g' |awk -vOFS='\t' '{print "chr"$1,$2,$1"_"$2"_SNP",$5,$7,"999\tPASS","ALTCOUNT="$8";REFCOUNT="$6";AC="$8";AN="($6+$8)}' |sort -k1,1V -k2,2n >1356_all.fakevcf

#manually add head: 1356_all.fake.vcf
#cat head 1356_all.fakevcf|bcftools view - -Oz -o 1356_all.fake.vcf.gz
#tabix 1356_all.fake.vcf.gz


## Make ancestral tables based on pairwise alignment
# liftOver 1356_snp_dm6.bed /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/2.lastz/pacbioSim/dm6.PacBioSim.all.chain.gz 1356_snp_dm6_pacbiosim.bed dsim.log
# liftOver 1356_snp_dm6.bed /gpfs2/zhangh/backup/database/phylogeny/dm6/dm6.droYak3.all.chain 1356_snp_dm6_dyak.bed dyak.log

# bedtools getfasta -fi /gpfs2/zhangh/backup/database/phylogeny/dm6/dm6.fa -bed 1356_snp_dm6.bed -fo 1356_snp_dm6.fa -name -tab -s 
# bedtools getfasta -fi /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/0.genome/Drosophila_simulans_Genome.fasta -bed 1356_snp_dm6_pacbiosim.bed -fo 1356_snp_dm6_pacbiosim.fa -name -tab -s 
# bedtools getfasta -fi /gpfs2/zhangh/backup/database/phylogeny/dm6/droYak3.fa -bed 1356_snp_dm6_dyak.bed -fo 1356_snp_dm6_dyak.fa -name -tab -s

# python3 mergeBedFastab.py 1356_snp_dm6.bed 1356_snp_dm6.fa >1356_snp_dm6.reftab &
# python3 mergeBedFastab.py 1356_snp_dm6_pacbiosim.bed 1356_snp_dm6_pacbiosim.fa >1356_snp_dm6_pacbiosim.reftab &
# python3 mergeBedFastab.py 1356_snp_dm6_dyak.bed 1356_snp_dm6_dyak.fa >1356_snp_dm6_dyak.reftab &
# python3 ancestral_base_upper.py 2way 1356_snp_dm6.reftab 1356_snp_dm6_pacbiosim.reftab >1356_snp_dm6_anc.txt #
# python3 ancestral_base_upper.py 3way 1356_snp_dm6.reftab 1356_snp_dm6_pacbiosim.reftab 1356_snp_dm6_dyak.reftab >1356_snp_dm6_anc3way.txt


## Extract polymorphic sites (remove those overlapping with CDS)
# filter CDS overlapping sites;
# sed 's/^/chr/' ../dmel-all-r6.04.gtf.CDS.bed| bedtools sort >dmel.r6.04.CDS.chrm.bed

# bcftools view -G -R ../short_intron_strict_neutral.chr.bed -v snps -q 0.00001:minor -O v 1356_all.fake.vcf.gz | bedtools intersect -a stdin -b ../dmel.r6.04.CDS.chrm.bed -v -header | gzip -c  >1356_short_intron_strict_neutral.vcf.gz &
# bcftools view -G -R ../short_intron_looser_neutral.chr.bed -v snps -q 0.00001:minor -O v 1356_all.fake.vcf.gz | bedtools intersect -a stdin -b ../dmel.r6.04.CDS.chrm.bed -v -header | gzip -c  >1356_short_intron_looser_neutral.vcf.gz &
# cat ../dmel-all-r6.04.gtf.5UTR.bed | sed 's/^/chr/' | bcftools view -G -R -  -v snps -q 0.00001:minor -O v 1356_all.fake.vcf.gz | bedtools intersect -a stdin -b ../dmel.r6.04.CDS.chrm.bed -v -header | gzip -c >1356_utr5.vcf.gz &
# cat ../dmel-all-r6.04.gtf.3UTR.bed | sed 's/^/chr/' | bcftools view -G -R -  -v snps -q 0.00001:minor -O v 1356_all.fake.vcf.gz | bedtools intersect -a stdin -b ../dmel.r6.04.CDS.chrm.bed -v -header | gzip -c >1356_utr3.vcf.gz &
# cat ../dmel-all-r6.04.gtf.CDS.bed | sed 's/^/chr/' | bcftools view -G -R -  -v snps -q 0.00001:minor -O z -H 1356_all.fake.vcf.gz >1356_cds.vcf.gz

## Determine the effect of polymorphic mutations
# sed 's/^/chr/' ../dmel-all-r6.04.gtf.5UTR.bed | bedtools intersect -a 1356_utr5.vcf.gz -b stdin -wo | python3 find_polymorphic_triplet.py /lustre/user/lulab/luoshq/syq/reference/dmel-all-transcript-r6.04.fasta >poly_diff_utr5.txt
# sed 's/^/chr/' ../dmel-all-r6.04.gtf.3UTR.bed | bedtools intersect -a 1356_utr3.vcf.gz -b stdin -wo | python3 find_polymorphic_triplet.py /lustre/user/lulab/luoshq/syq/reference/dmel-all-transcript-r6.04.fasta >poly_diff_utr3.txt
# zcat 1356_short_intron_strict_neutral.vcf.gz | bedtools intersect -a stdin -b ../short_intron_strict_neutral.chr.bed  -wo | python3 find_polymorphic_triplet_intron.py  ../short_intron_strict_neutral.fa >poly_diff_short_intron_strict_neutral.txt
# zcat 1356_short_intron_looser_neutral.vcf.gz | bedtools intersect -a stdin -b ../short_intron_looser_neutral.chr.bed  -wo | python3 find_polymorphic_triplet_intron.py  ../short_intron_looser_neutral.fa >poly_diff_short_intron_looser_neutral.txt
# not done: java -jar ~/lcl/software/snpEff/snpEff.jar -nodownload -canon -onlyProtein -no-downstream -no-upstream -no-intron -no-utr -no-intergenic -v BDGP6.82 1064_cds.vcf.gz | java -jar ~/lcl/software/snpEff/SnpSift.jar extractFields -s "," -e "." -  CHROM POS ID REF ALT "ANN[*].EFFECT" "ANN[*].HGVS_P" "ANN[*].GENEID" > 1064_cds.snpeff


#sed 's/ALTCOUNT=//g' 1356_all.fakevcf|sed 's/;REFCOUNT=/\t/g' |sed 's/;AC=/\t/g' | sed 's/;AN=/\t/g' |awk -vOFS='\t' '{print $1,$2,$3,$10/$11}' >1356_snp_dm6_all_coods_nohead.txt



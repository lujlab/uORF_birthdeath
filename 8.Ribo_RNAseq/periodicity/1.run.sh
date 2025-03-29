less ../../0.ann/dmel-all-exon-r6.04.gtf|sed 's/; /\t/g'|sed 's/"//g' |sed 's/ /\t/g' |cut -f 1,4,5,7,10,14|awk -vOFS='\t'  '{print $1,$2-1,$3,$5,$6,$4}' >dm6.04_exon_FBgn_FBtr.bed

awk -vOFS='\t' '$2=="CDS" {print $1,$3-1,$4}' /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/0.ann/tr-r6.04-cds_ut_tr-pos.txt >all_cds_tpos.bed
awk -vOFS='\t' '{print $1,$2,$3}' /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/1.STAR_new2/find_nonoverlaped_uATGpos/dmel_uORF_nogap_dmonly.bed.corrected|sed 's/dm6.//g' >all_uorf_tpos.bed

for i in `ls ../2.map3/coverage/*ribo_sorted_riboseq_strand.bedgraph`; do j=`basename $i _sorted_riboseq_strand.bedgraph`; bedtools intersect -a $i -b dm6.04_exon_FBgn_FBtr.bed -wa -wb |awk -vOFS='\t' '$6==$12{print $1,$2,$3,$11,$10,$6,$5}' >${j}_ribo_Psite.bed; python /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/1.STAR_new2/find_nonoverlaped_uORF_region/gpos2tpos_bed6.py /lustre/user/lulab/luoshq/syq/reference/dmel-all-r6.04.gtf ${j}_ribo_Psite.bed > ${j}_ribo_Psite.bed.tmp; awk -vOFS='\t' '{if ($8-1>=0) print $4,$8-1,$9,$7;else print $4,$8,$9,$7}' ${j}_ribo_Psite.bed.tmp|bedtools intersect -a - -b all_cds_tpos.bed -wa -wb|awk -vOFS='\t' '{print $0,$3-$6,($3-$6)%3}' >${j}_ribo_Psite.bed.frame;awk -vOFS='\t' '{if ($8-1>=0) print $4,$8-1,$9,$7;else print $4,$8,$9,$7}' ${j}_ribo_Psite.bed.tmp|bedtools intersect -a - -b all_uorf_tpos.bed -wa -wb|awk -vOFS='\t' '{print $0,$3-$6,($3-$6)%3}' >${j}_ribo_uORF_Psite.bed.frame;done


ls *_ribo_Psite.bed.frame >framefile;awk -vOFS='\t' '{print $1,$1}' framefile >framefile2 #manually add sample col
ls *_ribo_uORF_Psite.bed.frame >framefile_uORF;awk -vOFS='\t' '{print $1,$1}' framefile_uORF >framefile_uORF2 #manually add sample col


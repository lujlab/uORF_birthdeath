setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/4.count_rpkm_shift1")
library(data.table)

#dmel potential startcodon
count <- fread('/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.map3/all_potential_start_cnt_shift1/all_potential_start_cnt.tab',header=T)
colnames(count)
count=count[,c("trid","feature_len","name","S2_harr.Psite_potential_start_cnt.txt")]
f=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/10.GLOOME/non-canonical-start/total_potential_start_sort.bed",header=F)
f=cbind(f,count)
f=f[,c(4:6,10)]
colnames(f)=c("id","pos","codon","S2_harr")
f2=dcast(f,id+codon~pos)
f2[is.na(f2)]=0
f2$TIS=0
f2[(f2$pos0>f2$pos1 & f2$pos0>(f2$`pos-1`+f2$`pos-2`)),]$TIS=1
fwrite(f2,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/4.count_rpkm_shift1/TIS.txt")

f3=f2[f2$TIS==1,]

f4=f3[f3$pos0>3,]

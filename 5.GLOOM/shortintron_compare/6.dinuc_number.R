#sed 's/-//g' dm6_PacBioSim_27species_shortintron.fa>dm6_PacBioSim_27species_shortintron_nogap.fa
#~/lcl/software/ucsc-kent/faCount -dinuc dm6_PacBioSim_27species_shortintron_nogap.fa > dm6_PacBioSim_27species_shortintron.dinuc.txt
library(data.table)
f=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/14.shortintron_ATG/dm6_PacBioSim_27species_shortintron.dinuc.txt",header=T)
colnames(f)[1]="seq"
f=f[f$seq!="total",]
f[, c('sp', 'intronid') := tstrsplit(seq, "[.]", fixed = FALSE)]
sp=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/7.new_27_maf/pacbioSim/species_prefix_in_maf.txt",header=F)$V1
sp_out=data.table()

#sum(f[f$sp=="droSim1",]$len)
for (i in c(1:length(sp))){
  sp_tmp=sp[i]
  tmp=f[sp==sp_tmp,]
  col_sums <- as.data.table(t(colSums(tmp[,2:24])))
  
  colnames(col_sums)=colnames(tmp[,2:24])
  col_sums$sp=sp_tmp
  sp_out=rbind(sp_out,col_sums)
}
fwrite(sp_out,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/14.shortintron_ATG/shortintron_dinuc_stat_merged_27sp.txt",sep='\t',col.names = T)

#CG content
f$CG_content=(f$C+f$G)/f$len
f_mel=f[sp=="dm6",]

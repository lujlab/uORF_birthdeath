library(data.table)
library(dplyr)
f=fread("/data/tr_len_27sp/dm6_PacBioSim_27species_merged_5UTR_nongap.dinuc.txt",header=T)
colnames(f)[1]="seq"
f=f[f$seq!="total",]
f[, c('sp', 'FBid', 'number') := tstrsplit(seq, "[._]", fixed = FALSE)]
sp=fread("/data/maf/species_prefix_in_maf.txt",header=F)$V1

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
fwrite(sp_out,"/results/stat_merged_27sp.txt",sep='\t',col.names = T)


#for each gene in dm6
sp_tmp="dm6"
tmp=f[sp==sp_tmp,]
fwrite(tmp,"/results/dm6_block_stat.txt",sep='\t',col.names = T)
#python stat_dm6gene_5UTR_dinuc.py dm6_block_stat.txt dm6_block_stat_split.txt
#tmp=fread("/results/dm6_block_stat_split.txt",sep='\t')
col_sums <- as.data.table(t(colSums(tmp[,2:24])))
tmp2=tmp[, lapply(.SD, sum),by = FBid, .SDcols = colnames(tmp)[2:24]]
fwrite(tmp2,"/results/dm6_merge_gene_stat.txt",sep='\t',col.names = T)



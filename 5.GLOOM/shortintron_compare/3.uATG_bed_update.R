#setwd("~/syq/project/uorf_MK/gain_loss_syq/")
setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/14.shortintron_ATG")
library(data.table)
#intron bed
intron=fread("short_intron_strict_neutral.bed2",header=F,sep='\t')
intron=intron[,c(4,6)]
colnames(intron)=c("intronid","strand")
#### generate uATG bed ####
hh<-read.csv("./ATG_matrix_dm6_PacBioSim_27species_shortintron.csv",header = T,stringsAsFactors = F)
names(hh)[1]<-"id"
as.data.table(hh)->hh1

hh1[, intronid := sub("_\\d+\\(\\d+\\).*", "", id)]
hh1=merge(hh1,intron,by="intronid")

hh1[, c("seqname","start0","end","align_pos") := tstrsplit(id, "_")]
hh1[, c("position_withoutgap","start0","end","align_pos") := tstrsplit(id, "_")]
hh1[, position_withoutgap := sub(".*\\((\\d+)\\).*", "\\1", id)]
hh1$position_withoutgap=as.numeric(hh1$position_withoutgap)
hh1$start0=as.numeric(hh1$start0)
hh1$end=as.numeric(hh1$end)
hh1$gene_pos=hh1$start0+hh1$position_withoutgap
hh1[strand=="-",]$gene_pos=hh1[strand=="-",]$end-hh1[strand=="-",]$position_withoutgap
hh1[,align_pos:=NULL]

fwrite(hh1,"ATG_matrix_dm6_PacBioSim_27species_shortintron_update.csv",col.names = T,row.names = F,sep=",",quote = F)



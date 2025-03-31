library(data.table)
#for initiation site: 1) Psite peak 2) Ribo rpkm >10 for the uorf 3)expression mRNA rpkm > 1 for the NATGN region 4)not overlapped with cds
TIS=fread("/data/rpkm/TIS.txt",sep=',') #with overlapped with cds
colnames(TIS)[1]=c("uorf_id")

orfpos=fread("/data/rpkm/potential_start_id_pos_FBtr_noCDS.bed",sep='\t',header=F) #rmCDS
orfpos=orfpos[,c(4,1,5,7)]
colnames(orfpos)=c("uorf_id","seqname","gene_pos","codon")

f2=merge(orfpos,TIS,by=c("uorf_id","codon")) #rm start codon overlapped with cds

#initiation signal
is=fread("/data/rpkm/dmel_potential_uORF_start_signal_S2harr_rpkm.txt",sep='\t',header=T) #rmCDS
is=is[,c("id","codon","S2_NXXXN_Ribo","S2_NXXXN_mRNA","S2_NXXXN_mRNA_rpkm","S2_NXXXN_Ribo_rpkm","S2_initiationsignal")]
colnames(is)[1]="uorf_id"
f3=merge(f2,is,by=c("uorf_id","codon"))



# most abundant tr
dm6_most_tr=fread("/data/rpkm/tr_mRNA_merge_rank_S2.txt",header=T,sep='\t')
S2_most_tr=dm6_most_tr[S2_tr_rank==1 & S2_tr_mRNA_rpkm>1,]$trid

f3[, c("FBtr", "tmp2") := tstrsplit(uorf_id, "_", fixed = TRUE)]

tmp <- f3[f3$FBtr %in% S2_most_tr & f3$S2_NXXXN_mRNA>1,]

tmp2 <- tmp[tmp$TIS==1 & tmp$pos0>=3,]
#tmp2 = tmp

tmp2$color="red"
tmp2[codon!="ATG",]$color="grey"

library(ggplot2)
library(dplyr)
tmp_median_sorted <- tmp2 %>%
  group_by(codon) %>%
  summarize(median_value = median(S2_initiationsignal)) %>%
  arrange(median_value)
tmp2$codon <- factor(tmp2$codon, levels = rev(tmp_median_sorted$codon))

pdf("/results/figS27TIS.pdf",width=6.6,height=2.5)
ggplot(tmp2,aes(x=codon,y=S2_initiationsignal,fill=color))+
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian(ylim=c(0,200))+
  scale_fill_manual(values=c("#A9A9A9","#FF4040"))+
  theme_classic() #S2_initiation_signal
dev.off()

table(tmp2$codon)
#ATG  ATT  ATA  CTG  TTG  ACG  GTG  AAG  ATC  AGG 
#1602 2876 1996 1845 2290 2180 2431 2849 1847  966
wilcox.test(tmp[codon=="ATG",]$S2_initiationsignal,tmp[codon!="ATG",]$S2_initiationsignal) #p-value < 2.2e-16



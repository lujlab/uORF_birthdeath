library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(Rmisc)
library(RColorBrewer)
#group
uATG_group=fread("/data/rpkm/group/uATG_genopos_TEgroup_combined.txt",sep='\t')[,c("geno_pos","sum","group")]

#uORF_matrix
uORF_matrix=fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",sep=',',header=T)
colnames(uORF_matrix)[1]="uorf_align"
uORF_matrix=uORF_matrix[uORF_matrix$dm6>0 |uORF_matrix$PacBioSim>0 , ]
uORF_matrix=uORF_matrix[,c("uorf_align","dm6","PacBioSim","droYak3","transcriptID","position_withoutgap",
                           "seqname","geneID","gene_pos","genesymbol")]
uORF_matrix$geno_pos=paste0(uORF_matrix$seqname,"_",uORF_matrix$gene_pos)
uORF_matrix_group=merge(uORF_matrix,uATG_group,by="geno_pos")
uORF_matrix_group[, c("tmp1", "tmp2") := tstrsplit(uorf_align, "(", fixed = TRUE)]
uORF_matrix_group[, c("tmp1", "tmp2") := tstrsplit(tmp1, "_", fixed = TRUE)]
uORF_matrix_group$tmp2=as.numeric(uORF_matrix_group$tmp2)+1
uORF_matrix_group$uorf_id=paste(uORF_matrix_group$tmp1,uORF_matrix_group$tmp2,sep="_")
uORF_matrix_group=uORF_matrix_group[,c("uorf_align","uorf_id","dm6","transcriptID","position_withoutgap",
                                       "seqname","geneID","gene_pos","genesymbol","geno_pos","group","sum")]
# uORF TE
uORF_TE <- fread("/data/rpkm/dmel_uORF_mRNA_ribo_merge.txt",sep='\t',header=T)
uORF_matrix_TE=merge(uORF_matrix_group,uORF_TE,by="uorf_id") #rm overlap with cds


#BLS
BLS <- fread("/data/uORF_matrix_triCas2_ATG_updata_BLS.txt",header=T,sep='\t')
colnames(BLS)[1]="uorf_align"
BLS=BLS[dm6>0 |PacBioSim>0 , ]
#BLS=BLS[transcriptID %in% dm6_canon$FBtr,] #canon
BLS=BLS[,c("uorf_align","dm6","PacBioSim","BLS")]
BLS=merge(uORF_matrix_group,BLS,by=c("uorf_align","dm6"))

sample<-c("Kronja_mature_oocyte","Kronja_activated_egg",
          "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
          "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
          "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva","dmel_pupa",
          "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
          "Dunn_em02")
sample0<-c("Kronja_mature_oocyte","Kronja_activated_egg",
           "Samuels_GSC_AHS0","Samuels_GSC_AHS0.5","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
           "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
           "dmel_em_0_2h","dmel_em_2_6h","dmel_em_6_12h","dmel_em_12_24h","dmel_larva","dmel_pupa",
           "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
           "Dunn_em02")



####all BLS
BLS_tmp=BLS[,c("geno_pos","group","BLS","sum")]
BLS_tmp=BLS_tmp[!duplicated(BLS_tmp$geno_pos),]


#####cld
elements <- c("group1", "group2","group3", "group4")
combinations <- CJ(V1 = elements, V2 = elements)
combinations=combinations[V1<V2]
combinations$p=NA
for (i in c(1:nrow(combinations))){
  x=combinations$V1[i]
  y=combinations$V2[i]
  combinations$p[i]=wilcox.test(BLS_tmp[group==x,]$BLS,BLS_tmp[group==y,]$BLS)$p.value
}
median_table=data.table()
median_table <- BLS_tmp %>%
  group_by(group) %>%
  summarise(median = median(BLS))
median_table=median_table[order(-median_table$median),]

out=data.table(group=c("group1", "group2","group3", "group4"),letter="")
letterlist=c("a","b","c","d")
flag=1
for (g in median_table$group){
  if (out[out$group==g,]$letter==""){
    l1=combinations[(V1==g|V2==g)&p>0.05,]$V1
    l2=combinations[(V1==g|V2==g)&p>0.05,]$V2
    l=unique(c(l1,l2,g))
    out[group %in% l, letter := paste0(letter, letterlist[flag])]
    flag=flag+1
  }
}
pdf("/results/figS20_group_BLS.pdf",width=6.6,height=2.5)
ggplot(data=BLS_tmp,aes(x=group,y=BLS))+
  geom_violin(aes(x=group,y=BLS,fill=group),width=1)+
  geom_boxplot(fill="white",width=0.1,outlier.shape =NA)+
  geom_text(data=out,aes(x=group,label=letter,y=0.8),vjust = 2)+
  coord_cartesian(ylim=c(0,0.75)) + 
  #scale_fill_manual(values=c("#A9A9A9","#FF4040"))+
  theme_classic() #fig20

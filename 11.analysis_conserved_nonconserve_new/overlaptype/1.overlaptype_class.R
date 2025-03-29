library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(Rmisc)
#class
pos_class=fread("/data/rpkm/overlaptype/pos_uORF_class.txt",header=T)
pos_class$class="complex"
pos_class[nonoverlap==1,]$class="class1"
pos_class[inframe_overlapCDS==1 & outframe_overlapCDS==0 & nonoverlap==0 & inframe_overlapuORF==0 & outframe_overlapuORF==0,]$class="class2"
pos_class[inframe_overlapCDS==0 & outframe_overlapCDS==1 & nonoverlap==0 & inframe_overlapuORF==0 & outframe_overlapuORF==0,]$class="class3"
pos_class[inframe_overlapCDS==0 & outframe_overlapCDS==0 & nonoverlap==0 & inframe_overlapuORF==1 & outframe_overlapuORF==0,]$class="class4"
pos_class[inframe_overlapCDS==0 & outframe_overlapCDS==0 & nonoverlap==0 & inframe_overlapuORF==0 & outframe_overlapuORF==1,]$class="class5"
uATG_class=pos_class
#table(uATG_class$class)
# class1  class2  class3  class4  class5 complex 
# 12365   227     1053    4480    9952    8488 


#uORF_matrix
uORF_matrix <- fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",sep=',',header=T)
colnames(uORF_matrix)[1]="uorf_align"
uORF_matrix=uORF_matrix[uORF_matrix$dm6>0 |uORF_matrix$PacBioSim>0 , ]
uORF_matrix=uORF_matrix[,c("uorf_align","dm6","PacBioSim","droYak3","transcriptID","position_withoutgap",
                           "seqname","geneID","gene_pos","genesymbol")]
uORF_matrix$geno_pos=paste0(uORF_matrix$seqname,"_",uORF_matrix$gene_pos)
uORF_matrix_class=merge(uORF_matrix,uATG_class,by="geno_pos")
uORF_matrix_class[, c("tmp1", "tmp2") := tstrsplit(uorf_align, "(", fixed = TRUE)]
uORF_matrix_class[, c("tmp1", "tmp2") := tstrsplit(tmp1, "_", fixed = TRUE)]
uORF_matrix_class$tmp2=as.numeric(uORF_matrix_class$tmp2)+1
uORF_matrix_class$uorf_id=paste(uORF_matrix_class$tmp1,uORF_matrix_class$tmp2,sep="_")
uORF_matrix_class=uORF_matrix_class[,c("uorf_align","uorf_id","dm6","transcriptID","position_withoutgap",
                                       "seqname","geneID","gene_pos","genesymbol","geno_pos","class")]
# uORF TE
uORF_TE <- fread("/data/rpkm/dmel_uORF_mRNA_ribo_merge.txt",sep='\t',header=T)
uORF_matrix_TE=merge(uORF_matrix_class,uORF_TE,by="uorf_id") #rm overlap with cds

# CDS TE
temp_mRNA_ribo_merge <- fread("/data/rpkm/dmel_CDS_mRNA_ribo_merge.txt",sep='\t',header=T)
colnames(temp_mRNA_ribo_merge)[1]="FBtr"
# most abundant tr
dm6_most_tr=fread("/data/rpkm/tr_mRNA_merge_rank_21sample.txt",header=T,sep='\t')
mean_most=dm6_most_tr[dm6_most_tr$dmel_mean_tr_rank==1,]$trid
# kozak score
uORF_matrix_kozak <- fread("/data/uORF_kozakscore.txt2_new",header=T,sep='\t')

uORF_matrix_kozak=merge(uORF_matrix_class,uORF_matrix_kozak[,c("uorf_align","dm6_nogap_pos","kozak_seq_dm6","kozak_score_dm6","dsim_nogap_pos","kozak_seq_PacBioSim","kozak_score_PacBioSim")],by="uorf_align")#rm overlap with cds
#uORF_matrix_kozak=uORF_matrix_kozak[uORF_matrix_kozak$transcriptID %in% mean_most,] #mean most abundant

#BLS
BLS <- fread("/data/uORF_matrix_triCas2_ATG_updata_BLS.txt",header=T,sep='\t')
colnames(BLS)[1]="uorf_align"
BLS=BLS[dm6>0 |PacBioSim>0 , ]
#BLS=BLS[transcriptID %in% dm6_canon$FBtr,] #canon
BLS=BLS[,c("uorf_align","dm6","PacBioSim","BLS")]
BLS=merge(uORF_matrix_class,BLS,by=c("uorf_align","dm6"))

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

#all ATG
uORF_matrix_kozak_tmp=uORF_matrix_kozak
uORF_matrix_kozak_tmp_clean <- na.omit(uORF_matrix_kozak_tmp)
uORF_matrix_kozak_tmp_clean=uORF_matrix_kozak_tmp_clean[!duplicated(uORF_matrix_kozak_tmp_clean$geno_pos),]
#####cld
elements <- c("class1", "class2","class3", "class4","class5","complex")
combinations <- CJ(V1 = elements, V2 = elements)
combinations=combinations[V1<V2]
combinations$p=NA
for (i in c(1:nrow(combinations))){
  x=combinations$V1[i]
  y=combinations$V2[i]
  combinations$p[i]=wilcox.test(uORF_matrix_kozak_tmp_clean[class==x,]$kozak_score_dm6,uORF_matrix_kozak_tmp_clean[class==y,]$kozak_score_dm6)$p.value
}
median_table=data.table()
median_table <- uORF_matrix_kozak_tmp_clean %>%
  group_by(class) %>%
  summarise(median = median(kozak_score_dm6))
median_table=median_table[order(-median_table$median),]

out=data.table(class=c("class1", "class2","class3", "class4","class5","complex"),letter="")
letterlist=c("a","b","c","d","e","f")
flag=1
for (g in median_table$class){
  if (out[out$class==g,]$letter==""){
    l1=combinations[(V1==g|V2==g)&p>0.05,]$V1
    l2=combinations[(V1==g|V2==g)&p>0.05,]$V2
    l=unique(c(l1,l2,g))
    out[class %in% l, letter := paste0(letter, letterlist[flag])]
    flag=flag+1
  }
}

p=ggplot(data=uORF_matrix_kozak_tmp_clean,aes(x=class,y=kozak_score_dm6))+
  geom_violin(aes(x=class,y=kozak_score_dm6,fill=class),width=0.8)+
  geom_boxplot(fill="white",width=0.2,outlier.shape =NA)+
  geom_text(data=out,aes(x=class,label=letter,y=6))+
  #geom_hline(yintercept=median(uORF_matrix_kozak_tmp_clean[class=="class1",]$kozak_score_dm6),color="red",linetype="dashed")+
  #stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..))+
  #coord_cartesian(ylim=c(-6,6)) + 
  #scale_fill_manual(values=c("#A9A9A9","#FF4040"))+
  theme_classic() +
  theme(legend.position="none") #fig_kozak_overlaptype_letter
pdf("/results/fig_kozak_overlaptype_letter.pdf",width=4,height=4)
print(p)
dev.off()

####all BLS
#BLS_tmp=BLS[(transcriptID %in% mean_most),c("geno_pos","class","BLS")]
BLS_tmp=BLS[,c("geno_pos","class","BLS")]
BLS_tmp=BLS_tmp[!duplicated(BLS_tmp$geno_pos),]
#####cld
elements <- c("class1", "class2","class3", "class4","class5","complex")
combinations <- CJ(V1 = elements, V2 = elements)
combinations=combinations[V1<V2]
combinations$p=NA
for (i in c(1:nrow(combinations))){
  x=combinations$V1[i]
  y=combinations$V2[i]
  combinations$p[i]=wilcox.test(BLS_tmp[class==x,]$BLS,BLS_tmp[class==y,]$BLS)$p.value
}
median_table=data.table()
median_table <- BLS_tmp %>%
  group_by(class) %>%
  summarise(median = median(BLS))
median_table=median_table[order(-median_table$median),]

out=data.table(class=c("class1", "class2","class3", "class4","class5","complex"),letter="")
letterlist=c("a","b","c","d","e","f")
flag=1
for (g in median_table$class){
  if (out[out$class==g,]$letter==""){
    l1=combinations[(V1==g|V2==g)&p>0.05,]$V1
    l2=combinations[(V1==g|V2==g)&p>0.05,]$V2
    l=unique(c(l1,l2,g))
    out[class %in% l, letter := paste0(letter, letterlist[flag])]
    flag=flag+1
  }
}


p=ggplot(data=BLS_tmp,aes(x=class,y=BLS))+
  stat_boxplot(geom ='errorbar', width = 0.3 )+
  geom_boxplot(aes(x=class,y=BLS,fill=class),outlier.shape =NA)+
  coord_cartesian(ylim=c(0,0.75)) + 
  geom_text(data=out,aes(x=class,label=letter,y=0.8),vjust = 2)+
  #scale_fill_manual(values=c("#A9A9A9","#FF4040"))+
  theme_classic()+
  theme(legend.position="none")#fig_BLS_overlaptype_letter
pdf("/results/fig_BLS_overlaptype_letter.pdf",width=4,height=4)
print(p)
dev.off()

##TE for each sample
#TE of conserved uORF vs species specific uORF,only for mel
###cld
type_TE <- uORF_matrix_TE

for(i in 1:21){
  k <- sample[i]
  FBtr <- temp_mRNA_ribo_merge$FBtr[temp_mRNA_ribo_merge[[paste0(k,"_CDS_mRNA_rpkm")]]>1 ]
  tmp <- type_TE
  tmp <- tmp[tmp$transcriptID %in% FBtr,]
  tmp <- tmp[tmp$transcriptID %in% dm6_most_tr[get(paste0(k,"_tr_rank"))==1,]$trid,]
  selected_columns <- c("uorf_align","class",paste0(k, "_TE"),paste0( k, "_uORF_mRNA_rpkm"))
  TE <- tmp[,..selected_columns]
  TE=TE[TE[[paste0(k,"_uORF_mRNA_rpkm")]]>1]
  #TE[is.na(TE)] <- 0
  colnames(TE)=c("uorf_align","class","dmel","uORF_mRNA_rpkm")
  TE <- TE[,c("uorf_align","class","dmel")]
  TE2=melt(TE,id.vars=c("uorf_align","class"))
  elements <- c("class1", "class2","class3", "class4","class5","complex")
  combinations <- CJ(V1 = elements, V2 = elements)
  combinations=combinations[V1<V2]
  combinations$p=NA
  for (j in c(1:nrow(combinations))){
    x=combinations$V1[j]
    y=combinations$V2[j]
    combinations$p[j]=wilcox.test(TE2[class==x,]$value,TE2[class==y,]$value)$p.value
  }
  median_table=data.table()
  median_table <- TE2 %>%
    group_by(class) %>%
    summarise(median = median(value))
  median_table=median_table[order(-median_table$median),]
  
  out=data.table(class=c("class1", "class2","class3", "class4","class5","complex"),letter="")
  letterlist=c("a","b","c","d")
  flag=1
  for (g in median_table$class){
    if (out[out$class==g,]$letter==""){
      l1=combinations[(V1==g|V2==g)&p>0.05,]$V1
      l2=combinations[(V1==g|V2==g)&p>0.05,]$V2
      l=unique(c(l1,l2,g))
      out[class %in% l, letter := paste0(letter, letterlist[flag])]
      flag=flag+1
    }
  }
  nn=paste0("p",i)
  #compaired=list(c("mel-gain","sim-loss"),c("mel-loss","sim-gain"))
  p=ggplot(data=TE2,aes(x=class,y=value,fill=class))+
    #geom_boxplot(linetype = "dashed",outlier.shape = NA)+
    stat_boxplot(geom ='errorbar', width = 0.3 )+
    stat_boxplot(aes(ymin = after_stat(lower), ymax =after_stat(upper)),outlier.shape = NA)+
    coord_cartesian(ylim=c(0,8))+
    geom_text(data=out,aes(x=class,label=letter,y=8),vjust = 2)+
    labs(y="TE",title = k)+
    #geom_signif(comparisons = compaired,test = wilcox.test,y_position =1.8,tip_length = 0,map_signif_level=TRUE)+
    theme_classic()+
    theme(legend.position="none")
  assign(nn,p)
}
pdf("/results/fig_TE_typeoverlap_21sample_letter.pdf",width=20,height=12)
multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,p11, p12, p13, p14, p15, p16, p17, p18, p19, p20,p21,cols=5) #fig_TE_typeoverlap_21sample_letter,20*12
dev.off()


####kozak vs BLS
uORF_matrix_kozak_tmp=uORF_matrix_kozak
uORF_matrix_kozak_tmp_clean <- na.omit(uORF_matrix_kozak_tmp)
uORF_matrix_kozak_tmp_clean=uORF_matrix_kozak_tmp_clean[!duplicated(uORF_matrix_kozak_tmp_clean$geno_pos),c("geno_pos","class","kozak_score_dm6")]
BLS_tmp=BLS[,c("geno_pos","class","BLS")]
BLS_tmp=BLS_tmp[!duplicated(BLS_tmp$geno_pos),]
kozak_BLS=merge(uORF_matrix_kozak_tmp_clean,BLS_tmp,by=c("geno_pos","class"))
#for different type
kozak_BLS_cor <- matrix(ncol=6,nrow=2)
colnames(kozak_BLS_cor) <- c("class1","class2","class3","class4","class5","complex")
rownames(kozak_BLS_cor) <- c("Rho","p")
for(i in 1:6){
  k <- c("class1","class2","class3","class4","class5","complex")[i]
  c=cor.test(kozak_BLS[class==k,]$kozak_score_dm6,kozak_BLS[class==k,]$BLS,method="spearman",exact = FALSE)#rho=0.9497054,p-value < 2.2e-16
  kozak_BLS_cor["p",k]=c$p.value
  kozak_BLS_cor["Rho",k]=c$estimate
  
  nn=paste0("p",i)
  #compaired=list(c("mel-gain","sim-loss"),c("mel-loss","sim-gain"))
  p=ggplot(kozak_BLS[class==k,],aes(x=kozak_score_dm6,y=BLS))+
    geom_point(color=rgb(1,0.5,0,0.3))+
    #coord_cartesian(xlim=c(-6,5),ylim=c(-6,5))+
    theme_classic()+
    labs(x="kozak",y="BLS",title = k)
  assign(nn,p)
}
#kozak_BLS_cor
#class1       class2       class3       class4       class5      complex
#Rho 3.477132e-01 1.831059e-01 1.243040e-01 1.030206e-01 8.313532e-02 8.646318e-02
#p   1.157323e-06 7.755517e-09 4.837791e-42 7.310058e-12 2.099970e-16 3.192228e-15

pdf("/results/fig_kozak_BLS_overlaptype.pdf",width=6,height=4.3)
multiplot(p1, p2, p3, p4, p5, p6,cols=3) ##fig_kozak_BLS_overlaptype, 6*4.3
dev.off()


## uORF TE vs kozak, only mel
rho_kozak_TE <- matrix(ncol=21,nrow=6)
colnames(rho_kozak_TE) <- sample
rownames(rho_kozak_TE) <- c("class1","class2","class3","class4","class5","complex")
rho_kozak_TE_p <- matrix(ncol=21,nrow=6)
colnames(rho_kozak_TE_p) <- sample
rownames(rho_kozak_TE_p) <- c("class1","class2","class3","class4","class5","complex")
#
for(i in 1:21){
  k <- sample[i]
  k2 <- sample0[i]
  FBtr <- temp_mRNA_ribo_merge$FBtr[temp_mRNA_ribo_merge[[paste0(k,"_CDS_mRNA_rpkm")]]>1 ]
  tmp <- uORF_matrix_TE
  tmp <- tmp[tmp$transcriptID %in% FBtr,]
  #tmp <- tmp[tmp$transcriptID %in% dm6_canon$FBtr,]
  tmp <- tmp[tmp$transcriptID %in% dm6_most_tr[get(paste0(k,"_tr_rank"))==1,]$trid,]
  
  tmp2 <- merge(tmp,uORF_matrix_kozak,by=c("uorf_align","uorf_id","class","dm6","geno_pos"),sort=F)
  
  selected_columns <- c(paste0(k, "_uORF_mRNA_rpkm"),paste0( k, "_TE"), "kozak_score_dm6","class")
  x_dmel <- tmp2[tmp2$dm6==1,..selected_columns]
  x_dmel=x_dmel[x_dmel[[paste0(k,"_uORF_mRNA_rpkm")]]>1]
  cleaned_x_dmel <- na.omit(x_dmel)
  colnames(cleaned_x_dmel)=c("mRNA_rpkm","TE","kozak","class")
  #class
  rho_kozak_TE["class1",k] <- cor.test(cleaned_x_dmel[class=="class1",]$TE,cleaned_x_dmel[class=="class1",]$kozak,method="spearman")$estimate
  rho_kozak_TE_p["class1",k] <- cor.test(cleaned_x_dmel[class=="class1",]$TE,cleaned_x_dmel[class=="class1",]$kozak,method="spearman")$p.value
  rho_kozak_TE["class2",k] <- cor.test(cleaned_x_dmel[class=="class2",]$TE,cleaned_x_dmel[class=="class2",]$kozak,method="spearman")$estimate
  rho_kozak_TE_p["class2",k] <- cor.test(cleaned_x_dmel[class=="class2",]$TE,cleaned_x_dmel[class=="class2",]$kozak,method="spearman")$p.value
  rho_kozak_TE["class3",k] <- cor.test(cleaned_x_dmel[class=="class3",]$TE,cleaned_x_dmel[class=="class3",]$kozak,method="spearman")$estimate
  rho_kozak_TE_p["class3",k] <- cor.test(cleaned_x_dmel[class=="class3",]$TE,cleaned_x_dmel[class=="class3",]$kozak,method="spearman")$p.value
  rho_kozak_TE["class4",k] <- cor.test(cleaned_x_dmel[class=="class4",]$TE,cleaned_x_dmel[class=="class4",]$kozak,method="spearman")$estimate
  rho_kozak_TE_p["class4",k] <- cor.test(cleaned_x_dmel[class=="class4",]$TE,cleaned_x_dmel[class=="class4",]$kozak,method="spearman")$p.value
  rho_kozak_TE["class5",k] <- cor.test(cleaned_x_dmel[class=="class5",]$TE,cleaned_x_dmel[class=="class5",]$kozak,method="spearman")$estimate
  rho_kozak_TE_p["class5",k] <- cor.test(cleaned_x_dmel[class=="class5",]$TE,cleaned_x_dmel[class=="class5",]$kozak,method="spearman")$p.value
  rho_kozak_TE["complex",k] <- cor.test(cleaned_x_dmel[class=="complex",]$TE,cleaned_x_dmel[class=="complex",]$kozak,method="spearman")$estimate
  rho_kozak_TE_p["complex",k] <- cor.test(cleaned_x_dmel[class=="complex",]$TE,cleaned_x_dmel[class=="complex",]$kozak,method="spearman")$p.value
}
#plot kozak vs TE
#rho_kozak_TE2=t(as.data.table(rho_kozak_TE))
#rho_kozak_TE_p2=t(as.data.table(rho_kozak_TE_p))
#rho_kozak_TE=data.table(cbind(rho_kozak_TE2,rho_kozak_TE_p2))
#rho_kozak_TE$sample=sample

rho_kozak_TE2=as.data.table(rho_kozak_TE)
rho_kozak_TE2$class=rownames(rho_kozak_TE)
rho_kozak_TE=melt(rho_kozak_TE2,id.vars="class")
colnames(rho_kozak_TE)=c("class","stage","Rho")



rho_kozak_TE_p2=as.data.table(rho_kozak_TE_p)
rho_kozak_TE_p2$class=rownames(rho_kozak_TE_p)
rho_kozak_TE_p=melt(rho_kozak_TE_p2,id.vars="class")
colnames(rho_kozak_TE_p)=c("class","stage","p")

rho_kozak_TE=merge(rho_kozak_TE,rho_kozak_TE_p,by=c("class","stage"))

custom_order=rev(sample)
rho_kozak_TE$stage=factor(rho_kozak_TE$stage,levels=custom_order)
rho_kozak_TE$color="grey"
rho_kozak_TE[rho_kozak_TE$p<0.05,]$color="#52A4DF"

rho_kozak_TE_p$sig="ns"
rho_kozak_TE_p[rho_kozak_TE_p$p<0.001,]$sig="***"
rho_kozak_TE_p[rho_kozak_TE_p$p>0.001&rho_kozak_TE_p$p<0.01,]$sig="**"
rho_kozak_TE_p[rho_kozak_TE_p$p>0.01&rho_kozak_TE_p$p<0.05,]$sig="*"

###plot
for (i in c(1:6)){
  k=c("class1","class2","class3","class4","class5","complex")[i]
  tmp=rho_kozak_TE[rho_kozak_TE$class==k,]
  nn=paste0("p",i)
  p=ggplot(tmp, aes(x = stage, y = Rho)) +
    geom_bar(fill=tmp$color,stat = "identity", position = "dodge",width = 0.7) +
    labs(title = paste0("Kozak score vs TE, ",k),x = "Stage", y = "Rho") +
    #scale_fill_manual(values=c("#52A4DF"))+
    theme_classic()+
    coord_flip()
  assign(nn,p)
}

pdf("/results/fig_kozak_TE_overlaptype.pdf",width=16,height=8)
multiplot(p1, p2, p3, p4, p5, p6,cols=3) ##fig_kozak_TE_overlaptype, 16*8
dev.off()
######


## uORF TE vs BLS
rho_BLS_TE <- matrix(ncol=21,nrow=6)
colnames(rho_BLS_TE) <- sample
rownames(rho_BLS_TE) <- c("class1","class2","class3","class4","class5","complex")
rho_BLS_TE_p <- matrix(ncol=21,nrow=6)
colnames(rho_BLS_TE_p) <- sample
rownames(rho_BLS_TE_p) <- c("class1","class2","class3","class4","class5","complex")
#
for(i in 1:21){
  k <- sample[i]
  k2 <- sample0[i]
  FBtr <- temp_mRNA_ribo_merge$FBtr[temp_mRNA_ribo_merge[[paste0(k,"_CDS_mRNA_rpkm")]]>1 ]
  tmp <- uORF_matrix_TE
  tmp <- tmp[tmp$transcriptID %in% FBtr,]
  #tmp <- tmp[tmp$transcriptID %in% dm6_canon$FBtr,]
  tmp <- tmp[tmp$transcriptID %in% dm6_most_tr[get(paste0(k,"_tr_rank"))==1,]$trid,]
  
  tmp2 <- merge(tmp,BLS[,c("uorf_align","uorf_id","class","dm6","geno_pos","BLS")],by=c("uorf_align","uorf_id","class","dm6","geno_pos"),sort=F)
  
  selected_columns <- c(paste0(k, "_uORF_mRNA_rpkm"),paste0( k, "_TE"), "BLS","class")
  x_dmel <- tmp2[tmp2$dm6==1,..selected_columns]
  x_dmel=x_dmel[x_dmel[[paste0(k,"_uORF_mRNA_rpkm")]]>1]
  cleaned_x_dmel <- na.omit(x_dmel)
  colnames(cleaned_x_dmel)=c("mRNA_rpkm","TE","BLS","class")
  #class
  rho_BLS_TE["class1",k] <- cor.test(cleaned_x_dmel[class=="class1",]$TE,cleaned_x_dmel[class=="class1",]$BLS,method="spearman")$estimate
  rho_BLS_TE_p["class1",k] <- cor.test(cleaned_x_dmel[class=="class1",]$TE,cleaned_x_dmel[class=="class1",]$BLS,method="spearman")$p.value
  rho_BLS_TE["class2",k] <- cor.test(cleaned_x_dmel[class=="class2",]$TE,cleaned_x_dmel[class=="class2",]$BLS,method="spearman")$estimate
  rho_BLS_TE_p["class2",k] <- cor.test(cleaned_x_dmel[class=="class2",]$TE,cleaned_x_dmel[class=="class2",]$BLS,method="spearman")$p.value
  rho_BLS_TE["class3",k] <- cor.test(cleaned_x_dmel[class=="class3",]$TE,cleaned_x_dmel[class=="class3",]$BLS,method="spearman")$estimate
  rho_BLS_TE_p["class3",k] <- cor.test(cleaned_x_dmel[class=="class3",]$TE,cleaned_x_dmel[class=="class3",]$BLS,method="spearman")$p.value
  rho_BLS_TE["class4",k] <- cor.test(cleaned_x_dmel[class=="class4",]$TE,cleaned_x_dmel[class=="class4",]$BLS,method="spearman")$estimate
  rho_BLS_TE_p["class4",k] <- cor.test(cleaned_x_dmel[class=="class4",]$TE,cleaned_x_dmel[class=="class4",]$BLS,method="spearman")$p.value
  rho_BLS_TE["class5",k] <- cor.test(cleaned_x_dmel[class=="class5",]$TE,cleaned_x_dmel[class=="class5",]$BLS,method="spearman")$estimate
  rho_BLS_TE_p["class5",k] <- cor.test(cleaned_x_dmel[class=="class5",]$TE,cleaned_x_dmel[class=="class5",]$BLS,method="spearman")$p.value
  rho_BLS_TE["complex",k] <- cor.test(cleaned_x_dmel[class=="complex",]$TE,cleaned_x_dmel[class=="complex",]$BLS,method="spearman")$estimate
  rho_BLS_TE_p["complex",k] <- cor.test(cleaned_x_dmel[class=="complex",]$TE,cleaned_x_dmel[class=="complex",]$BLS,method="spearman")$p.value
}
#plot kozak vs TE
#rho_kozak_TE2=t(as.data.table(rho_kozak_TE))
#rho_kozak_TE_p2=t(as.data.table(rho_kozak_TE_p))
#rho_kozak_TE=data.table(cbind(rho_kozak_TE2,rho_kozak_TE_p2))
#rho_kozak_TE$sample=sample

rho_BLS_TE2=as.data.table(rho_BLS_TE)
rho_BLS_TE2$class=rownames(rho_BLS_TE)
rho_BLS_TE=melt(rho_BLS_TE2,id.vars="class")
colnames(rho_BLS_TE)=c("class","stage","Rho")



rho_BLS_TE_p2=as.data.table(rho_BLS_TE_p)
rho_BLS_TE_p2$class=rownames(rho_BLS_TE_p)
rho_BLS_TE_p=melt(rho_BLS_TE_p2,id.vars="class")
colnames(rho_BLS_TE_p)=c("class","stage","p")

rho_BLS_TE=merge(rho_BLS_TE,rho_BLS_TE_p,by=c("class","stage"))

custom_order=rev(sample)
rho_BLS_TE$stage=factor(rho_BLS_TE$stage,levels=custom_order)
rho_BLS_TE$color="grey"
rho_BLS_TE[rho_BLS_TE$p<0.05,]$color="#52A4DF"

rho_BLS_TE_p$sig="ns"
rho_BLS_TE_p[rho_BLS_TE_p$p<0.001,]$sig="***"
rho_BLS_TE_p[rho_BLS_TE_p$p>0.001&rho_BLS_TE_p$p<0.01,]$sig="**"
rho_BLS_TE_p[rho_BLS_TE_p$p>0.01&rho_BLS_TE_p$p<0.05,]$sig="*"
###plot
for (i in c(1:6)){
  k=c("class1","class2","class3","class4","class5","complex")[i]
  tmp=rho_BLS_TE[rho_BLS_TE$class==k,]
  nn=paste0("p",i)
  p=ggplot(tmp, aes(x = stage, y = Rho)) +
    geom_bar(fill=tmp$color,stat = "identity", position = "dodge",width = 0.7) +
    labs(title = paste0("BLS vs TE, ",k),x = "Stage", y = "Rho") +
    #scale_fill_manual(values=c("#52A4DF"))+
    theme_classic()+
    coord_flip()
  assign(nn,p)
}
pdf("/results/fig_BLS_TE_overlaptype.pdf",width=16,height=8)
multiplot(p1, p2, p3, p4, p5, p6,cols=3) ##fig_BLS_TE_overlaptype, 16*8
dev.off()

library(data.table)
library(ggplot2)
#library(ggsignif)
library(plyr)
library(Rmisc)
#uORF_matrix
uORF_matrix <- fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",sep=',',header=T)
colnames(uORF_matrix)[1]="uorf_align"
uORF_matrix=uORF_matrix[uORF_matrix$dm6>0 |uORF_matrix$PacBioSim>0 , ]
uORF_matrix=uORF_matrix[,c("uorf_align","dm6","PacBioSim","droYak3","transcriptID","position_withoutgap",
                           "seqname","geneID","gene_pos","genesymbol")]
uORF_matrix[, c("tmp1", "tmp2") := tstrsplit(uorf_align, "(", fixed = TRUE)]
uORF_matrix[, c("tmp1", "tmp2") := tstrsplit(tmp1, "_", fixed = TRUE)]
uORF_matrix$tmp2=as.numeric(uORF_matrix$tmp2)+1
uORF_matrix$uorf_id=paste(uORF_matrix$tmp1,uORF_matrix$tmp2,sep="_")
uORF_matrix=uORF_matrix[,c("uorf_align","uorf_id","dm6","PacBioSim","droYak3","transcriptID","position_withoutgap",
                           "seqname","geneID","gene_pos","genesymbol")]
# uORF TE
uORF_TE <- fread("/data/rpkm/dmel_uORF_mRNA_ribo_merge.txt",sep='\t',header=T)
uORF_matrix_TE=merge(uORF_matrix,uORF_TE,by="uorf_id") #rm overlap with cds

# CDS TE
temp_mRNA_ribo_merge <- fread("/data/rpkm/dmel_CDS_mRNA_ribo_merge.txt",sep='\t',header=T)
colnames(temp_mRNA_ribo_merge)[1]="FBtr"
# most abundant tr
dm6_most_tr=fread("/data/rpkm/tr_mRNA_merge_rank_21sample.txt",header=T,sep='\t')
mean_most=dm6_most_tr[dm6_most_tr$dmel_mean_tr_rank==1,]$trid
# kozak score
uORF_matrix_kozak <- fread("/data/uORF_kozakscore.txt2_new",header=T,sep='\t')

uORF_matrix_kozak=merge(uORF_matrix,uORF_matrix_kozak[,c("uorf_align","dm6_nogap_pos","kozak_seq_dm6","kozak_score_dm6","dsim_nogap_pos","kozak_seq_PacBioSim","kozak_score_PacBioSim")])#rm overlap with cds
#uORF_matrix_kozak=uORF_matrix_kozak[uORF_matrix_kozak$transcriptID %in% mean_most,] #mean most abundant


#BLS
BLS <- fread("/data/uORF_matrix_triCas2_ATG_updata_BLS.txt",header=T,sep='\t')
colnames(BLS)[1]="uorf_align"
BLS=BLS[dm6>0 |PacBioSim>0 , ]
#BLS=BLS[transcriptID %in% dm6_canon$FBtr,] #canon
BLS=BLS[,c("uorf_align","dm6","PacBioSim","BLS")]


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

## uORF TE vs kozak, only mel
rho_kozak_TE <- matrix(ncol=21,nrow=1)
colnames(rho_kozak_TE) <- sample
rownames(rho_kozak_TE) <- c("dmel")
rho_kozak_TE_p <- matrix(ncol=21,nrow=1)
colnames(rho_kozak_TE_p) <- sample
rownames(rho_kozak_TE_p) <- c("dmel")
#
for(i in 1:21){
  k <- sample[i]
  k2 <- sample0[i]
  FBtr <- temp_mRNA_ribo_merge$FBtr[temp_mRNA_ribo_merge[[paste0(k,"_CDS_mRNA_rpkm")]]>1 ]
  tmp <- uORF_matrix_TE
  tmp <- tmp[tmp$transcriptID %in% FBtr,]
  #tmp <- tmp[tmp$transcriptID %in% dm6_canon$FBtr,]
  tmp <- tmp[tmp$transcriptID %in% dm6_most_tr[get(paste0(k,"_tr_rank"))==1,]$trid,]
  tmp2 <- merge(tmp,uORF_matrix_kozak,by=c("uorf_align","uorf_id","dm6","PacBioSim"),sort=F)
  selected_columns <- c(paste0(k, "_uORF_mRNA_rpkm"),paste0( k, "_TE"), "kozak_score_dm6")
  x_dmel <- tmp2[tmp2$dm6==1,..selected_columns]
  x_dmel=x_dmel[x_dmel[[paste0(k,"_uORF_mRNA_rpkm")]]>1]
  cleaned_x_dmel <- na.omit(x_dmel)
  #nrow(cleaned_x_dmel[cleaned_x_dmel$Samuels_GSC_AHS05_TE>0,])

  rho_kozak_TE["dmel",k] <- cor.test(cleaned_x_dmel[[2]],cleaned_x_dmel[[3]],method="spearman")$estimate
  rho_kozak_TE_p["dmel",k] <- cor.test(cleaned_x_dmel[[2]],cleaned_x_dmel[[3]],method="spearman")$p.value
}
#plot kozak vs TE
#rho_kozak_TE2=t(as.data.table(rho_kozak_TE))
#rho_kozak_TE_p2=t(as.data.table(rho_kozak_TE_p))
#rho_kozak_TE=data.table(cbind(rho_kozak_TE2,rho_kozak_TE_p2))
#rho_kozak_TE$sample=sample

rho_kozak_TE2=as.data.table(rho_kozak_TE)
rho_kozak_TE2$species=rownames(rho_kozak_TE)
rho_kozak_TE=melt(rho_kozak_TE2,id.vars="species")
colnames(rho_kozak_TE)=c("species","stage","Rho")

#max(rho_kozak_TE_p) #0.05051736

rho_kozak_TE_p2=as.data.table(rho_kozak_TE_p)
rho_kozak_TE_p2$species=rownames(rho_kozak_TE_p)
rho_kozak_TE_p=melt(rho_kozak_TE_p2,id.vars="species")
colnames(rho_kozak_TE_p)=c("species","stage","p")

rho_kozak_TE=merge(rho_kozak_TE,rho_kozak_TE_p,by=c("species","stage"))

custom_order=rev(sample)
rho_kozak_TE$stage=factor(rho_kozak_TE$stage,levels=custom_order)
rho_kozak_TE_mel=rho_kozak_TE[rho_kozak_TE$species=="dmel",]
rho_kozak_TE_mel$color="grey"
rho_kozak_TE_mel[rho_kozak_TE_mel$p<0.05,]$color="#52A4DF"
p=ggplot(rho_kozak_TE_mel, aes(x = stage, y = Rho)) +
  geom_bar(fill=rho_kozak_TE_mel$color,stat = "identity", position = "dodge",width = 0.7) +
  labs(title = "Kozak score vs TE",x = "Stage", y = "Rho") +
  #scale_fill_manual(values=c("#52A4DF"))+
  theme_classic()+
  #scale_x_discrete(labels = rev(c("embryo 0-2h","embryo 2-6h","embryo 6-12h","embryo 12-24h","larva","pupa","female body","male body","female head","male head")))+
  coord_flip() #fig2A_kozakTE_21sample

pdf("/results/fig2e_kozakTE_21sample.pdf",width=4,height=6)
print(p)
dev.off()


# conserved uATG, Dmel vs Dsim kozak
uORF_matrix_kozak_conserved=uORF_matrix_kozak[uORF_matrix_kozak$dm6==1 &uORF_matrix_kozak$PacBioSim==1,]
uORF_matrix_kozak_conserved <- uORF_matrix_kozak_conserved[!is.na(uORF_matrix_kozak_conserved$kozak_score_dm6) & !is.na(uORF_matrix_kozak_conserved$kozak_score_PacBioSim),]
uORF_matrix_kozak_conserved=uORF_matrix_kozak_conserved[uORF_matrix_kozak_conserved$transcriptID %in% mean_most,]
ggplot(uORF_matrix_kozak_conserved,aes(x=uORF_matrix_kozak_conserved$kozak_score_dm6,y=uORF_matrix_kozak_conserved$kozak_score_PacBioSim))+
  geom_point(color=rgb(1,0.5,0,0.3))+
  coord_cartesian(xlim=c(-6,5),ylim=c(-6,5))+
  theme_classic()+
  labs(x="D. mel",y="D. sim",title = "Kozak score of conserved uORF") #fig2B_kozak_2sp
#3.5*3.5
cor.test(uORF_matrix_kozak_conserved$kozak_score_dm6,uORF_matrix_kozak_conserved$kozak_score_PacBioSim,method="spearman")#rho=0.9502957,p-value < 2.2e-16
cor.test(uORF_matrix_kozak_conserved$kozak_score_dm6,uORF_matrix_kozak_conserved$kozak_score_PacBioSim,method="pearson") #0.9518118


#compare species-specific and conserved uORF kozak
uORF_matrix_kozak_tmp=uORF_matrix_kozak[uORF_matrix_kozak$transcriptID %in% mean_most,c("uorf_align","dm6","PacBioSim","droYak3","dm6_nogap_pos","kozak_seq_dm6","kozak_score_dm6","dsim_nogap_pos","kozak_seq_PacBioSim","kozak_score_PacBioSim")]
uORF_matrix_kozak_tmp$type="unknown"
uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$dm6==1 & uORF_matrix_kozak_tmp$PacBioSim==1,]$type="conserved"
uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$dm6==1 & uORF_matrix_kozak_tmp$PacBioSim==0 & uORF_matrix_kozak_tmp$droYak3==0,]$type="mel_gain"
uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$dm6==1 & uORF_matrix_kozak_tmp$PacBioSim==0 & uORF_matrix_kozak_tmp$droYak3==1,]$type="sim_loss"

uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$dm6==0 & uORF_matrix_kozak_tmp$PacBioSim==1 & uORF_matrix_kozak_tmp$droYak3==0,]$type="sim_gain"
uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$dm6==0 & uORF_matrix_kozak_tmp$PacBioSim==1 & uORF_matrix_kozak_tmp$droYak3==1,]$type="mel_loss"


uORF_matrix_kozak_tmp2=uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="conserved"|uORF_matrix_kozak_tmp$type=="mel_gain"|uORF_matrix_kozak_tmp$type=="sim_loss",c("kozak_score_dm6","type")]
uORF_matrix_kozak_tmp3=uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="conserved"|uORF_matrix_kozak_tmp$type=="sim_gain"|uORF_matrix_kozak_tmp$type=="mel_loss",c("kozak_score_PacBioSim","type")]
colnames(uORF_matrix_kozak_tmp2)=c("kozak_score","type")
colnames(uORF_matrix_kozak_tmp3)=c("kozak_score","type")
uORF_matrix_kozak_tmp2[uORF_matrix_kozak_tmp2$type=="conserved",]$type="mel_cons"
uORF_matrix_kozak_tmp3[uORF_matrix_kozak_tmp3$type=="conserved",]$type="sim_cons"
uORF_matrix_kozak_tmp=rbind(uORF_matrix_kozak_tmp2,uORF_matrix_kozak_tmp3)
uORF_matrix_kozak_tmp$type=factor(uORF_matrix_kozak_tmp$type,levels=c("mel_gain","sim_loss","mel_cons","sim_cons","mel_loss","sim_gain"))
uORF_matrix_kozak_tmp$c="grey"
uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_cons" | uORF_matrix_kozak_tmp$type=="mel_cons" ,]$c="red"
uORF_matrix_kozak_tmp=uORF_matrix_kozak_tmp[!is.na(uORF_matrix_kozak_tmp$kozak_score),]

table(uORF_matrix_kozak_tmp$type)

p=ggplot(data=uORF_matrix_kozak_tmp,aes(x=type,y=kozak_score))+
  geom_violin(aes(x=type,y=kozak_score,fill=c),width=0.8)+
  geom_boxplot(fill="white",width=0.2,outlier.shape =NA)+
  #stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..))+
  #coord_cartesian(ylim=c(-6,6)) +
  scale_fill_manual(values=c("#A9A9A9","#FF4040"),guide="none")+
  theme_classic() #fig2d
pdf("/results/fig2d.pdf",width=4,height=4)
print(p)
dev.off()

wilcox.test(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_cons", ]$kozak_score,uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_loss", ]$kozak_score) #0.003784 **
wilcox.test(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_cons", ]$kozak_score,uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_gain", ]$kozak_score) #p-value < 2.2e-16 ***
wilcox.test(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_loss", ]$kozak_score,uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_gain", ]$kozak_score) #0.137
wilcox.test(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_cons", ]$kozak_score,uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_loss"|uORF_matrix_kozak_tmp$type=="mel_gain", ]$kozak_score) #p-value < 2.2e-16 ***

#mean(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_cons", ]$kozak_score)
#mean(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_loss", ]$kozak_score)
#mean(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_gain", ]$kozak_score)

#mean(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_cons", ]$kozak_score)
#mean(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_loss", ]$kozak_score)
#mean(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_gain", ]$kozak_score)

#median(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_cons", ]$kozak_score)
#median(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_loss", ]$kozak_score)
#median(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_gain", ]$kozak_score)

#median(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_cons", ]$kozak_score)
#median(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_loss", ]$kozak_score)
#median(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_gain", ]$kozak_score)

wilcox.test(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_cons", ]$kozak_score,uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_loss", ]$kozak_score) #0.3811
wilcox.test(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_cons", ]$kozak_score,uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_gain", ]$kozak_score) #7.486e-07 ***
wilcox.test(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_loss", ]$kozak_score,uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_gain", ]$kozak_score) #0.06343 
wilcox.test(uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="sim_cons", ]$kozak_score,uORF_matrix_kozak_tmp[uORF_matrix_kozak_tmp$type=="mel_loss"|uORF_matrix_kozak_tmp$type=="sim_gain", ]$kozak_score) #4.209e-06 ***



#TE of conserved uORF vs species specific uORF,only for mel
type_TE <- uORF_matrix_TE
type_TE$type="unknown"
type_TE<- type_TE[type_TE$dm6==1 | type_TE$PacBioSim==1,]
type_TE[type_TE$dm6==1 & type_TE$PacBioSim==1,]$type="conserved"
type_TE[type_TE$dm6==1 & type_TE$PacBioSim==0 & type_TE$droYak3==0,]$type="mel-gain"
type_TE[type_TE$dm6==1 & type_TE$PacBioSim==0 & type_TE$droYak3==1,]$type="sim-loss"
type_TE=type_TE[type_TE$type=="conserved"|type_TE$type=="mel-gain"|type_TE$type=="sim-loss",]
table(type_TE$type)

TE_p <- matrix(ncol=21,nrow=2)
colnames(TE_p) <- sample
rownames(TE_p) <- c("dmel_cons","dmel_gain-sim_loss")

for(i in 1:21){
  k <- sample[i]
  FBtr <- temp_mRNA_ribo_merge$FBtr[temp_mRNA_ribo_merge[[paste0(k,"_CDS_mRNA_rpkm")]]>1 ]
  tmp <- type_TE
  tmp <- tmp[tmp$transcriptID %in% FBtr,]
  tmp <- tmp[tmp$transcriptID %in% dm6_most_tr[get(paste0(k,"_tr_rank"))==1,]$trid,]
  selected_columns <- c("uorf_align","type",paste0(k, "_TE"),paste0( k, "_uORF_mRNA_rpkm"))
  TE <- tmp[,..selected_columns]
  TE=TE[TE[[paste0(k,"_uORF_mRNA_rpkm")]]>1]
  #TE[is.na(TE)] <- 0
  colnames(TE)=c("uorf_align","type","dmel","uORF_mRNA_rpkm")
  TE <- TE[,c("uorf_align","type","dmel")]
  TE2=melt(TE,id.vars=c("uorf_align","type"))
  TE2$c="grey"
  TE2[(TE2$type=="conserved" ),]$c="red"
  TE2$type=factor(TE2$type,levels=c("mel-gain","sim-loss","conserved"))
  TE_p["dmel_cons",k] <- wilcox.test(TE2[TE2$type=="mel-gain" |TE2$type=="sim-loss", ]$value,TE2[TE2$type=="conserved", ]$value)$p.value
  TE_p["dmel_gain-sim_loss",k] <- wilcox.test(TE2[TE2$type=="mel-gain" , ]$value,TE2[TE2$type=="sim-loss", ]$value)$p.value
  nn=paste0("p",i)
  #compaired=list(c("mel-gain","sim-loss"),c("mel-loss","sim-gain"))
  p=ggplot(data=TE2,aes(x=type,y=value,fill=c))+
    #geom_boxplot(linetype = "dashed",outlier.shape = NA)+
    stat_boxplot(geom ='errorbar', width = 0.3 )+
    stat_boxplot(aes(ymin = after_stat(lower), ymax =after_stat(upper)),outlier.shape = NA)+
    #coord_cartesian(ylim=c(0,20))+
    coord_cartesian(ylim=c(0,8))+
    scale_fill_manual(values=c("#A9A9A9","#FF4040"))+
    labs(y="TE",title = k)+
    #geom_signif(comparisons = compaired,test = wilcox.test,y_position =1.8,tip_length = 0,map_signif_level=TRUE)+
    theme_classic()
  assign(nn,p)
}

k <- sample[1]
FBtr <- temp_mRNA_ribo_merge$FBtr[temp_mRNA_ribo_merge[[paste0(k,"_CDS_mRNA_rpkm")]]>1 ]
tmp <- type_TE
tmp <- tmp[tmp$transcriptID %in% FBtr,]
tmp <- tmp[tmp$transcriptID %in% dm6_most_tr[get(paste0(k,"_tr_rank"))==1,]$trid,]
selected_columns <- c("uorf_align","type",paste0(k, "_TE"),paste0( k, "_uORF_mRNA_rpkm"))
TE <- tmp[,..selected_columns]
TE=TE[TE[[paste0(k,"_uORF_mRNA_rpkm")]]>1]
#TE[is.na(TE)] <- 0
colnames(TE)=c("uorf_align","type","dmel","uORF_mRNA_rpkm")
TE <- TE[,c("uorf_align","type","dmel")]
TE2=melt(TE,id.vars=c("uorf_align","type"))
TE2$c="grey"
TE2[(TE2$type=="conserved" ),]$c="red"
TE2$type=factor(TE2$type,levels=c("mel-gain","sim-loss","conserved"))
TE_p["dmel_cons",k] <- wilcox.test(TE2[TE2$type=="mel-gain" |TE2$type=="sim-loss", ]$value,TE2[TE2$type=="conserved", ]$value)$p.value
TE_p["dmel_gain-sim_loss",k] <- wilcox.test(TE2[TE2$type=="mel-gain" , ]$value,TE2[TE2$type=="sim-loss", ]$value)$p.value
#compaired=list(c("mel-gain","sim-loss"),c("mel-loss","sim-gain"))
p1=ggplot(data=TE2,aes(x=type,y=value,fill=c))+
    #geom_boxplot(linetype = "dashed",outlier.shape = NA)+
    stat_boxplot(geom ='errorbar', width = 0.3 )+
    stat_boxplot(aes(ymin = after_stat(lower), ymax =after_stat(upper)),outlier.shape = NA)+
    #coord_cartesian(ylim=c(0,20))+
    coord_cartesian(ylim=c(0,15))+
    scale_fill_manual(values=c("#A9A9A9","#FF4040"))+
    labs(y="TE",title = k)+
    #geom_signif(comparisons = compaired,test = wilcox.test,y_position =1.8,tip_length = 0,map_signif_level=TRUE)+
    theme_classic()


pdf("/results/fig2f_S22.pdf",width=17,height=8)
multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,p11, p12, p13, p14, p15, p16, p17, p18, p19, p20,p21,cols=6) #17*8
dev.off()

#TE_p

#uORF TE and BLS correlation
uORF_TE_BLS <- matrix(ncol=21,nrow=1)
colnames(uORF_TE_BLS) <- sample
rownames(uORF_TE_BLS) <- c("dmel")
uORF_TE_BLS_p=uORF_TE_BLS
for(i in 1:21){
  k <- sample[i]
  FBtr <- temp_mRNA_ribo_merge$FBtr[temp_mRNA_ribo_merge[[paste0(k,"_CDS_mRNA_rpkm")]]>1]
  tmp <- uORF_matrix_TE
  tmp <- tmp[tmp$transcriptID %in% FBtr,]
  tmp <- tmp[tmp$transcriptID %in% dm6_most_tr[get(paste0(k,"_tr_rank"))==1,]$trid,]
  tmp2 <- merge(tmp,BLS,by=c("uorf_align","dm6","PacBioSim"),sort=F)
  
  selected_columns <- c(paste0( k, "_TE"), "BLS",paste0(k,"_uORF_mRNA_rpkm"))
  x_dmel <- tmp2[dm6==1,..selected_columns]
  x_dmel=x_dmel[x_dmel[[paste0(k,"_uORF_mRNA_rpkm")]]>1]
  
  # Dmel TE vs kozak
  uORF_TE_BLS["dmel",k] <- cor.test(x_dmel[[2]],x_dmel[[1]],method="spearman")$estimate
  uORF_TE_BLS_p["dmel",k] <- cor.test(x_dmel[[2]],x_dmel[[1]],method="spearman")$p.value
}


#plot TE vs BLS
uORF_TE_BLS2=as.data.table(uORF_TE_BLS)
uORF_TE_BLS2$species=rownames(uORF_TE_BLS)
uORF_TE_BLS=melt(uORF_TE_BLS2,id.vars="species")
colnames(uORF_TE_BLS)=c("species","stage","Rho")
#max(uORF_TE_BLS_p) # 0.2653976
#     Kronja_mature_oocyte Kronja_activated_egg Samuels_GSC_AHS0 Samuels_GSC_AHS05 Samuels_GSC_AHS5 Samuels_GSC_AHS9
#dmel         4.647605e-07         6.710012e-12     2.699548e-12      5.880322e-08     0.0002864002     0.0007893596
#Samuels_GSC_AHS18 Samuels_GSC_AHS28 Samuels_GSC_AHS38 Douka_S2_UTD2    dmel_em02 dmel_em26   dmel_em612
#dmel      4.078332e-11      1.314942e-09      4.274251e-09  4.146083e-13 2.291684e-05 0.2653976 8.121934e-31
#dmel_em1224   dmel_larva    dmel_pupa dmel_female_body dmel_male_body dmel_female_head dmel_male_head
#dmel 1.43787e-77 2.329219e-10 3.165374e-28     1.876207e-26    1.30888e-18     3.869656e-30   4.967681e-29
#Dunn_em02
#dmel 7.367453e-26
uORF_TE_BLS_p2=as.data.table(uORF_TE_BLS_p)
uORF_TE_BLS_p2$species=rownames(uORF_TE_BLS_p)
uORF_TE_BLS_p=melt(uORF_TE_BLS_p2,id.vars="species")
colnames(uORF_TE_BLS_p)=c("species","stage","p")

uORF_TE_BLS=merge(uORF_TE_BLS,uORF_TE_BLS_p,by=c("species","stage"))

custom_order=rev(sample)
uORF_TE_BLS$stage=factor(uORF_TE_BLS$stage,levels=custom_order)
uORF_TE_BLS$color="grey"
uORF_TE_BLS[uORF_TE_BLS$p<0.05,]$color="#52A4DF"

p=ggplot(uORF_TE_BLS, aes(x = stage, y = Rho)) +
  geom_bar(stat = "identity", position = "dodge",fill=uORF_TE_BLS$color,width=0.7) +
  labs(title = "TE vs BLS",x = "Stage", y = "Rho") +
  #scale_fill_manual(values=c("#478AFA","#FF4040"))+
  theme_classic()+
  #scale_x_discrete(labels = rev(c("embryo 0-2h","embryo 2-6h","embryo 6-12h","embryo 12-24h","larva","pupa","female body","male body","female head","male head")))+
  coord_flip() #fig2E_TE_BLS_21sample
pdf("/results/fig2a_TE_BLS_21sample.pdf",width=4,height=6)
print(p)
dev.off()

#4*5
#eg: male head
k = "dmel_male_head"
FBtr <- temp_mRNA_ribo_merge$FBtr[temp_mRNA_ribo_merge[[paste0(k,"_CDS_mRNA_rpkm")]]>1]
tmp <- uORF_matrix_TE
tmp <- tmp[tmp$transcriptID %in% FBtr,]
tmp <- tmp[tmp$transcriptID %in% dm6_most_tr[get(paste0(k,"_tr_rank"))==1,]$trid,]
tmp2 <- merge(tmp,BLS,by=c("uorf_align","dm6","PacBioSim"),sort=F)
selected_columns <- c(paste0(k, "_TE"), "BLS",paste0(k,"_uORF_mRNA_rpkm"))
x_dmel <- tmp2[dm6==1,..selected_columns]
x_dmel=x_dmel[x_dmel[[paste0(k,"_uORF_mRNA_rpkm")]]>1]
#x_dmel[is.na(x_dmel)] <- 0
cor.test(x_dmel[[2]],x_dmel[[1]],method="spearman",exact = FALSE) #0.1131372 ,p-value < 2.2e-16

x_dmel$log2_TE=log2(x_dmel$dmel_male_head_TE)
x_dmel$log2_BLS=log2(x_dmel$BLS)
p=ggplot(x_dmel,aes(x=log2(BLS),y=log2(dmel_male_head_TE)))+
  geom_point(color="#FF4040",alpha=0.3)+
  theme_classic()+
  #coord_cartesian(xlim=c(-5,5),ylim=c(-5,0))+
  labs(x="log2(BLS)",y="log2(TE)",title = "male head")
#4*4 #fig2F_TE_BLS_MH
pdf("/results/fig2b_TE_BLS_MH.pdf",width=4,height=4)
print(p)
dev.off()








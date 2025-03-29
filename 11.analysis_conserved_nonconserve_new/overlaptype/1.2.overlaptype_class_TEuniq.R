library(data.table)
setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/5.result_shift1/overlaptype_compare")
#class
pos_class=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/13.uORF_overlap/overlap_with_cds/pos_uORF_class.txt",header=T)
pos_class$class="complex"
pos_class[nonoverlap==1,]$class="class1"
pos_class[inframe_overlapCDS==1 & outframe_overlapCDS==0 & nonoverlap==0 & inframe_overlapuORF==0 & outframe_overlapuORF==0,]$class="class2"
pos_class[inframe_overlapCDS==0 & outframe_overlapCDS==1 & nonoverlap==0 & inframe_overlapuORF==0 & outframe_overlapuORF==0,]$class="class3"
pos_class[inframe_overlapCDS==0 & outframe_overlapCDS==0 & nonoverlap==0 & inframe_overlapuORF==1 & outframe_overlapuORF==0,]$class="class4"
pos_class[inframe_overlapCDS==0 & outframe_overlapCDS==0 & nonoverlap==0 & inframe_overlapuORF==0 & outframe_overlapuORF==1,]$class="class5"
uATG_class=pos_class
#table(uATG_class$class)
# class1  class2  class3  class4  class5 complex 
#12365     227    1053    4480    9952    8488


#uORF_matrix
uORF_matrix <- fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/8.new_27_mfa/PacBioSim/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",sep=',',header=T)
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
# uORF TE,这里面没去掉和cds有overlap的
uORF_TE <- fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/4.count_rpkm_shift1/dmel_uORF_uniq_mRNA_ribo_merge.txt",sep='\t',header=T)
colnames(uORF_TE)[2]="uorf_align"
uORF_matrix_TE=merge(uORF_matrix_class,uORF_TE,by="uorf_align") #rm overlap with cds

# CDS TE
temp_mRNA_ribo_merge <- fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/4.count_rpkm_shift1/dmel_CDS_mRNA_ribo_merge.txt",sep='\t',header=T)
colnames(temp_mRNA_ribo_merge)[1]="FBtr"
# most abundant tr
dm6_most_tr=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/4.count_rpkm_shift1/tr_mRNA_merge_rank_21sample.txt",header=T,sep='\t')
mean_most=dm6_most_tr[dm6_most_tr$dmel_mean_tr_rank==1,]$trid
# kozak score
uORF_matrix_kozak <- fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/3.kozak/uORF_kozakscore.txt2_new",header=T,sep='\t')

uORF_matrix_kozak=merge(uORF_matrix_class,uORF_matrix_kozak[,c("uorf_align","dm6_nogap_pos","kozak_seq_dm6","kozak_score_dm6","dsim_nogap_pos","kozak_seq_PacBioSim","kozak_score_PacBioSim")],by="uorf_align")#rm overlap with cds
#uORF_matrix_kozak=uORF_matrix_kozak[uORF_matrix_kozak$transcriptID %in% mean_most,] #mean most abundant

#canonical
#canonical <- fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/3.kozak/uORF_kozakscore.txt2_new",header=T,sep='\t')

#BLS
BLS <- fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/5.BLS/uORF_matrix_triCas2_ATG_updata_BLS.txt",header=T,sep='\t')
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



####compare kozak
#most tr
#uORF_matrix_kozak_tmp=uORF_matrix_kozak[uORF_matrix_kozak$transcriptID %in% mean_most,c("uorf_align","geno_pos","class","kozak_score_dm6")]
#uORF_matrix_kozak_tmp_clean <- na.omit(uORF_matrix_kozak_tmp)

#ggplot(data=uORF_matrix_kozak_tmp_clean,aes(x=class,y=kozak_score_dm6))+
#  geom_violin(aes(x=class,y=kozak_score_dm6,fill=class),width=0.8)+
#  geom_boxplot(fill="white",width=0.2,outlier.shape =NA)+
#  #stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..))+
#  #coord_cartesian(ylim=c(-6,6)) + #不会删除点
#  #scale_fill_manual(values=c("#A9A9A9","#FF4040"))+
#  theme_classic() #fig_kozak_class

##TE for each sample
#TE of conserved uORF vs species specific uORF,only for mel
type_TE <- uORF_matrix_TE

TE_p <- matrix(ncol=21,nrow=15)
colnames(TE_p) <- sample
rownames(TE_p) <- c("class1_2","class1_3","class1_4","class1_5","class1_c","class2_3","class2_4","class2_5","class2_c",
                    "class3_4","class3_5","class3_c","class4_5","class4_c","class5_c")
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
  TE_p["class1_2",k] <- wilcox.test(TE2[TE2$class=="class1", ]$value,TE2[TE2$class=="class2", ]$value)$p.value
  TE_p["class1_3",k] <- wilcox.test(TE2[TE2$class=="class1", ]$value,TE2[TE2$class=="class3", ]$value)$p.value
  TE_p["class1_4",k] <- wilcox.test(TE2[TE2$class=="class1", ]$value,TE2[TE2$class=="class4", ]$value)$p.value
  TE_p["class1_5",k] <- wilcox.test(TE2[TE2$class=="class1", ]$value,TE2[TE2$class=="class5", ]$value)$p.value
  TE_p["class1_c",k] <- wilcox.test(TE2[TE2$class=="class1", ]$value,TE2[TE2$class=="complex", ]$value)$p.value
  TE_p["class2_3",k] <- wilcox.test(TE2[TE2$class=="class2", ]$value,TE2[TE2$class=="class3", ]$value)$p.value
  TE_p["class2_4",k] <- wilcox.test(TE2[TE2$class=="class2", ]$value,TE2[TE2$class=="class4", ]$value)$p.value
  TE_p["class2_5",k] <- wilcox.test(TE2[TE2$class=="class2", ]$value,TE2[TE2$class=="class5", ]$value)$p.value
  TE_p["class2_c",k] <- wilcox.test(TE2[TE2$class=="class2", ]$value,TE2[TE2$class=="complex", ]$value)$p.value
  TE_p["class3_4",k] <- wilcox.test(TE2[TE2$class=="class3", ]$value,TE2[TE2$class=="class4", ]$value)$p.value
  TE_p["class3_5",k] <- wilcox.test(TE2[TE2$class=="class3", ]$value,TE2[TE2$class=="class5", ]$value)$p.value
  TE_p["class3_c",k] <- wilcox.test(TE2[TE2$class=="class3", ]$value,TE2[TE2$class=="complex", ]$value)$p.value
  TE_p["class4_5",k] <- wilcox.test(TE2[TE2$class=="class4", ]$value,TE2[TE2$class=="class5", ]$value)$p.value
  TE_p["class4_c",k] <- wilcox.test(TE2[TE2$class=="class4", ]$value,TE2[TE2$class=="complex", ]$value)$p.value
  TE_p["class5_c",k] <- wilcox.test(TE2[TE2$class=="class5", ]$value,TE2[TE2$class=="complex", ]$value)$p.value
  nn=paste0("p",i)
  #compaired=list(c("mel-gain","sim-loss"),c("mel-loss","sim-gain"))
  p=ggplot(data=TE2,aes(x=class,y=value,fill=class))+
    #geom_boxplot(linetype = "dashed",outlier.shape = NA)+
    stat_boxplot(geom ='errorbar', width = 0.3 )+
    stat_boxplot(aes(ymin = after_stat(lower), ymax =after_stat(upper)),outlier.shape = NA)+
    coord_cartesian(ylim=c(0,8))+
    labs(y="TE",title = k)+
    #geom_signif(comparisons = compaired,test = wilcox.test,y_position =1.8,tip_length = 0,map_signif_level=TRUE)+
    theme_classic()
  assign(nn,p)
}
multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,p11, p12, p13, p14, p15, p16, p17, p18, p19, p20,p21,cols=5) #fig_TE_typeoverlap_21sample_TEuniq,20*12

TE_p2=as.data.frame(TE_p)
TE_p2$compare=row.names(TE_p)
fwrite(TE_p2,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/5.result_shift1/overlaptype_compare/TE_class_p_TEuniq.txt",sep='\t')
TE_p3=TE_p2
TE_p3[TE_p3<0.001]="***"
TE_p3[TE_p3<0.01 & TE_p3!="***"]="**"
TE_p3[TE_p3<0.05 & TE_p3!="***" & TE_p3!="**"]="*"




###TE median
TE_p <- matrix(ncol=21,nrow=6)
colnames(TE_p) <- sample
rownames(TE_p) <- c("class1","class2","class3","class4","class5","complex")
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
  TE_p["class1",k] <- median(TE2[TE2$class=="class1", ]$value)
  TE_p["class2",k] <- median(TE2[TE2$class=="class2", ]$value)
  TE_p["class3",k] <- median(TE2[TE2$class=="class3", ]$value)
  TE_p["class4",k] <- median(TE2[TE2$class=="class4", ]$value)
  TE_p["class5",k] <- median(TE2[TE2$class=="class5", ]$value)
  TE_p["complex",k] <- median(TE2[TE2$class=="complex", ]$value)
  
}

TE_p2=as.data.frame(TE_p)
TE_p2$compare=row.names(TE_p)
fwrite(TE_p2,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/5.result_shift1/overlaptype_compare/TE_median_class_TEuniq.txt",sep='\t')


#class13=rbind(TE_p3[TE_p3$compare=="class1_3",],TE_p2[TE_p2$compare=="class1" | TE_p2$compare=="class3",]) #class1<=class3
#class23=rbind(TE_p3[TE_p3$compare=="class2_3",],TE_p2[TE_p2$compare=="class2" | TE_p2$compare=="class3",]) #class2<class3
#class43=rbind(TE_p3[TE_p3$compare=="class3_4",],TE_p2[TE_p2$compare=="class4" | TE_p2$compare=="class3",]) #?
#class53=rbind(TE_p3[TE_p3$compare=="class3_5",],TE_p2[TE_p2$compare=="class5" | TE_p2$compare=="class3",]) #?
#classc3=rbind(TE_p3[TE_p3$compare=="class3_c",],TE_p2[TE_p2$compare=="complex" | TE_p2$compare=="class3",]) #complex<=class3


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

multiplot(p1, p2, p3, p4, p5, p6,cols=3) ##fig_kozak_TE_overlaptype_TEuniq, 16*8

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

multiplot(p1, p2, p3, p4, p5, p6,cols=3) ##fig_BLS_TE_overlaptype_TEuniq, 16*8

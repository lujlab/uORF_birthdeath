
library(data.table)
library(ggplot2)
library(plyr)
library(Rmisc)
#uORF_matrix
uORF_matrix <- fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS_postoCDS.csv")
colnames(uORF_matrix)[2]="uorf_align"
uORF_matrix=uORF_matrix[uORF_matrix$dm6>0 |uORF_matrix$PacBioSim>0 , ]
uORF_matrix=uORF_matrix[uORF_matrix$pos_to_cds<=300,]
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
#max(uORF_TE_BLS_p) # 0.7826774

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
  coord_flip() #fig2E_TE_BLS_21sample_300nt
pdf("/results/figS19a_TE_BLS_21sample_300nt.pdf",width=4,height=6)
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
cor.test(x_dmel[[2]],x_dmel[[1]],method="spearman",exact = FALSE) #0.09378252 ,p-value = 0.00000 00000 00000 4603

x_dmel$log2_TE=log2(x_dmel$dmel_male_head_TE)
x_dmel$log2_BLS=log2(x_dmel$BLS)
p=ggplot(x_dmel,aes(x=log2(BLS),y=log2(dmel_male_head_TE)))+
  geom_point(color="#FF4040",alpha=0.3)+
  theme_classic()+
  #coord_cartesian(xlim=c(-5,5),ylim=c(-5,0))+
  labs(x="log2(BLS)",y="log2(TE)",title = "male head")
#4*4 #fig2F_TE_BLS_MH_300nt
pdf("/results/figS19b_TE_BLS_MH_300nt.pdf",width=4,height=4)
print(p)
dev.off()









library(data.table)
library(ggplot2)
library(ggsignif)
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
# uORF TE,with cds overlapped
uORF_TE <- fread("/data/rpkm/dmel_uORF_mRNA_ribo_merge.txt",sep='\t',header=T)
uORF_matrix_TE=merge(uORF_matrix,uORF_TE,by="uorf_id") #rm overlap with cds

uORF_matrix_TE_mel=uORF_matrix_TE[uORF_matrix_TE$dm6==1,]
#unique(uORF_matrix_TE[uORF_matrix_TE$dm6==1,c("seqname","gene_pos")])

#####for each sample, how many uORFs are translated
# CDS TE-FC dmel-dsim
temp_mRNA_ribo_merge <- fread("/data/rpkm/dmel_CDS_mRNA_ribo_merge.txt",sep='\t',header=T)
colnames(temp_mRNA_ribo_merge)[1]="FBtr"

# most abundant tr
dm6_most_tr=fread("/data/rpkm/tr_mRNA_merge_rank_21sample.txt",header=T,sep='\t')


##### all transcripts were used ####
k<-c("Kronja_mature_oocyte","Kronja_activated_egg",
     "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
     "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
     "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva","dmel_pupa",
     "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
     "Dunn_em02")
k0<-c("Kronja_mature_oocyte","Kronja_activated_egg",
      "Samuels_GSC_AHS0","Samuels_GSC_AHS0.5","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
      "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
      "dmel_em_0_2h","dmel_em_2_6h","dmel_em_6_12h","dmel_em_12_24h","dmel_larva","dmel_pupa",
      "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
      "Dunn_em02")
rpkm1_num=matrix(ncol=21,nrow=1)
TE01_num=matrix(ncol=21,nrow=1)
TE05_num=matrix(ncol=21,nrow=1)
rpkm1_con=matrix(ncol=21,nrow=1)
TE01_con=matrix(ncol=21,nrow=1)
TE05_con=matrix(ncol=21,nrow=1)

colnames(rpkm1_num) <- k
rownames(rpkm1_num) <- c("dmel")

colnames(TE01_num) <- k
rownames(TE01_num) <- c("dmel")

colnames(TE05_num) <- k
rownames(TE05_num) <- c("dmel")

colnames(rpkm1_con) <- k
rownames(rpkm1_con) <- c("dmel")

colnames(TE01_con) <- k
rownames(TE01_con) <- c("dmel")

colnames(TE05_con) <- k
rownames(TE05_con) <- c("dmel")

#
for(i in 1:21){
  k<-c("Kronja_mature_oocyte","Kronja_activated_egg",
       "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
       "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
       "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva","dmel_pupa",
       "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
       "Dunn_em02")[i]
  FBtr <- temp_mRNA_ribo_merge$FBtr[temp_mRNA_ribo_merge[[paste0(k,"_CDS_mRNA_rpkm")]]>1 ]
  tmp <- uORF_matrix_TE
  tmp <- tmp[tmp$transcriptID %in% FBtr,]
  #tmp <- tmp[tmp$transcriptID %in% dm6_canon$FBtr,]
  selected_columns <- c(paste0(k, "_uORF_mRNA_rpkm"),paste0( k, "_TE"),"seqname","gene_pos")
  x_dmel <- tmp[tmp$dm6==1,..selected_columns]
  x_dmel=x_dmel[x_dmel[[paste0(k,"_uORF_mRNA_rpkm")]]>1]
  x_dmel2=x_dmel[x_dmel[[paste0(k, "_TE")]]>0.1]
  x_dmel3=x_dmel[x_dmel[[paste0(k, "_TE")]]>0.5]

  rpkm1_num["dmel",k] <- x_dmel[!duplicated(x_dmel[,c("seqname","gene_pos")]),.N]
  TE01_num["dmel",k] <- x_dmel2[!duplicated(x_dmel2[,c("seqname","gene_pos")]),.N]
  TE05_num["dmel",k] <- x_dmel3[!duplicated(x_dmel3[,c("seqname","gene_pos")]),.N]
}

x=cbind(t(rpkm1_num),t(TE01_num),t(TE05_num))

x=data.table(x)
colnames(x)=c("expressed_RPKM1", "translated_TE0.1","translated_TE0.5")
x$percentage01=x$translated_TE0.1/x$expressed_RPKM1
x$percentage05=x$translated_TE0.5/x$expressed_RPKM1

x$sample=c("Kronja_mature_oocyte","Kronja_activated_egg",
           "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
           "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
           "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva","dmel_pupa",
           "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
           "Dunn_em02")
#fwrite(x,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/4.count_rpkm_shift1/translated_uORF_total_percentage_21sample.txt",sep='\t')



#for total in 21 samples

t=data.table()
for(i in 1:21){
  k<-c("Kronja_mature_oocyte","Kronja_activated_egg",
       "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
       "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
       "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva","dmel_pupa",
       "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
       "Dunn_em02")[i]
  FBtr <- temp_mRNA_ribo_merge$FBtr[temp_mRNA_ribo_merge[[paste0(k,"_CDS_mRNA_rpkm")]]>1 ]
  tmp <- uORF_matrix_TE
  tmp <- tmp[tmp$transcriptID %in% FBtr,]
  #tmp <- tmp[tmp$transcriptID %in% dm6_canon$FBtr,]
  selected_columns <- c(paste0(k, "_uORF_mRNA_rpkm"),paste0( k, "_TE"),"seqname","gene_pos")
  x_dmel <- tmp[tmp$dm6==1,..selected_columns]
  colnames(x_dmel)=c("uORF_mRNA_rpkm","uORF_TE","seqname","gene_pos")
  x_dmel=x_dmel[uORF_mRNA_rpkm>1,]
  x_dmel$sample=k
  t=rbind(t,x_dmel)
}


#total rpkm>1 in at least one sample
t[!duplicated(t[,c("seqname","gene_pos")]),.N] #34527

#total rpkm>1 & TE>0.1 in at least one sample
t2=t[uORF_mRNA_rpkm>1 &uORF_TE>0.1,]
t2[!duplicated(t2[,c("seqname","gene_pos")]),.N] #29844

#total rpkm>1 & TE>0.5 in at least one sample
t3=t[uORF_mRNA_rpkm>1 &uORF_TE>0.5,]
t3[!duplicated(t3[,c("seqname","gene_pos")]),.N] #27468

x_total=data.table(`expressed_RPKM1`=t[!duplicated(t[,c("seqname","gene_pos")]),.N],
                   `translated_TE0.1`= t2[!duplicated(t2[,c("seqname","gene_pos")]),.N],
                   `translated_TE0.5`=t3[!duplicated(t3[,c("seqname","gene_pos")]),.N],
                   `sample`="total")
x_total$percentage01=x_total$translated_TE0.1/x_total$expressed_RPKM1
x_total$percentage05=x_total$translated_TE0.5/x_total$expressed_RPKM1

x2=rbind(x,x_total)


fwrite(x2,"/results/tableS2_translated_uORF_total_percentage_21sample.txt",sep='\t')
#fwrite(t,"/results/translated_uORF_21sample.txt",sep='\t')





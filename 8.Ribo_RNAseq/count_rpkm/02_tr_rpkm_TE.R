setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/4.count_rpkm_shift1")
library(data.table)
#mel libsize
libsize1<- fread('/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.map_new/map2genome/mrna_libsize',header=F)
colnames(libsize1)=c("lib","sample","mapped")
libsize2<- fread('/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.map_new/lib_size/ribo_libsize',header=F)
colnames(libsize2)=c("lib","sample","mapped")
libsize=rbind(libsize1,libsize2)
libsize=libsize[,c("sample","mapped")]
libsize_merge=libsize[, .(mapped = sum(mapped)), by = .(sample)]
libsize_merge[libsize_merge$sample=="Douka_S2_mrna",]$sample="Douka_S2_UTD2_mrna"
libsize_merge[libsize_merge$sample=="Luo_S2_cell_normal_starvation_ribo",]$sample="Luo_S2_cell_normal_ribo"

#dmel 
dmel_tr_count <- fread('/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.map3/all_tr_cnt/all_tr_cnt.tab')
colnames(dmel_tr_count)
dmel_tr_count$`dmel_female_body.weighted_tr_cnt.bed`=dmel_tr_count$`dmel_female_body1.weighted_tr_cnt.bed`+dmel_tr_count$`dmel_female_body2.weighted_tr_cnt.bed`
dmel_tr_count$`dmel_male_body.weighted_tr_cnt.bed`=dmel_tr_count$`dmel_male_body1.weighted_tr_cnt.bed`+dmel_tr_count$`dmel_male_body2.weighted_tr_cnt.bed`

dmel_tr_count$`Dunn_em02.weighted_tr_cnt.bed`=dmel_tr_count$`Dunn_em02_rep1.weighted_tr_cnt.bed`+dmel_tr_count$`Dunn_em02_rep2.weighted_tr_cnt.bed`

dmel_tr_count$`Samuels_GSC_AHS0.5.weighted_tr_cnt.bed`=dmel_tr_count$`Samuels_GSC_AHS0.5_rep1.weighted_tr_cnt.bed`+dmel_tr_count$`Samuels_GSC_AHS0.5_rep2.weighted_tr_cnt.bed`

dmel_tr_count$`Samuels_GSC_AHS0.weighted_tr_cnt.bed`=dmel_tr_count$`Samuels_GSC_AHS0_rep1.weighted_tr_cnt.bed`+dmel_tr_count$`Samuels_GSC_AHS0_rep2.weighted_tr_cnt.bed`
dmel_tr_count$`Samuels_GSC_AHS18.weighted_tr_cnt.bed`=dmel_tr_count$`Samuels_GSC_AHS18_rep1.weighted_tr_cnt.bed`+dmel_tr_count$`Samuels_GSC_AHS18_rep2.weighted_tr_cnt.bed`
dmel_tr_count$`Samuels_GSC_AHS28.weighted_tr_cnt.bed`=dmel_tr_count$`Samuels_GSC_AHS28_rep1.weighted_tr_cnt.bed`+dmel_tr_count$`Samuels_GSC_AHS28_rep2.weighted_tr_cnt.bed`
dmel_tr_count$`Samuels_GSC_AHS38.weighted_tr_cnt.bed`=dmel_tr_count$`Samuels_GSC_AHS38_rep1.weighted_tr_cnt.bed`+dmel_tr_count$`Samuels_GSC_AHS38_rep2.weighted_tr_cnt.bed`
dmel_tr_count$`Samuels_GSC_AHS5.weighted_tr_cnt.bed`=dmel_tr_count$`Samuels_GSC_AHS5_rep1.weighted_tr_cnt.bed`+dmel_tr_count$`Samuels_GSC_AHS5_rep2.weighted_tr_cnt.bed`
dmel_tr_count$`Samuels_GSC_AHS9.weighted_tr_cnt.bed`=dmel_tr_count$`Samuels_GSC_AHS9_rep1.weighted_tr_cnt.bed`+dmel_tr_count$`Samuels_GSC_AHS9_rep2.weighted_tr_cnt.bed`

dmel_tr_count=dmel_tr_count[,c("trid","tr_len","Kronja_mature_oocyte.weighted_tr_cnt.bed",
                               "Kronja_activated_egg.weighted_tr_cnt.bed",
                               "Samuels_GSC_AHS0.weighted_tr_cnt.bed",
                               "Samuels_GSC_AHS0.5.weighted_tr_cnt.bed",
                               "Samuels_GSC_AHS5.weighted_tr_cnt.bed",
                               "Samuels_GSC_AHS9.weighted_tr_cnt.bed",
                               "Samuels_GSC_AHS18.weighted_tr_cnt.bed",
                               "Samuels_GSC_AHS28.weighted_tr_cnt.bed",
                               "Samuels_GSC_AHS38.weighted_tr_cnt.bed",
                               "Douka_S2.weighted_tr_cnt.bed",
                               "dmel_em_0_2h.weighted_tr_cnt.bed",
                               "dmel_em_2_6h.weighted_tr_cnt.bed",
                               "dmel_em_6_12h.weighted_tr_cnt.bed",
                               "dmel_em_12_24h.weighted_tr_cnt.bed",
                               "dmel_larva.weighted_tr_cnt.bed",
                               "dmel_pupa.weighted_tr_cnt.bed",
                               "dmel_female_body.weighted_tr_cnt.bed",
                               "dmel_male_body.weighted_tr_cnt.bed",
                               "dmel_female_head.weighted_tr_cnt.bed",
                               "dmel_male_head.weighted_tr_cnt.bed",
                               "S2.weighted_tr_cnt.bed",
                               "Luo_S2_cell_normal.weighted_tr_cnt.bed",
                               "Luo_S2_cell_serum_starvation.weighted_tr_cnt.bed",
                               "Luo_S2_cell_ago2KD.weighted_tr_cnt.bed" ,
                               "Luo_S2_cell_ago2KD_serum_starvation.weighted_tr_cnt.bed",
                               "Luo_S2_cell_DMSO.weighted_tr_cnt.bed",
                               "Luo_S2_cell_dsNC.weighted_tr_cnt.bed",
                               "Luo_S2_cell_rapamycin.weighted_tr_cnt.bed",
                               "Luo_S2_cell_ssNC.weighted_tr_cnt.bed",
                               "Luo_S2_cell_T3.weighted_tr_cnt.bed",
                               "Luo_S2_cell_T6.weighted_tr_cnt.bed",
                               "Luo_S2_cell_T10.weighted_tr_cnt.bed",
                               "Dunn_em02.weighted_tr_cnt.bed")]

colnames(dmel_tr_count)=c("trid","dmel_tr_len","Kronja_mature_oocyte_tr_mRNA",
                          "Kronja_activated_egg_tr_mRNA",
                          "Samuels_GSC_AHS0_tr_mRNA",
                          "Samuels_GSC_AHS05_tr_mRNA",
                          "Samuels_GSC_AHS5_tr_mRNA",
                          "Samuels_GSC_AHS9_tr_mRNA",
                          "Samuels_GSC_AHS18_tr_mRNA",
                          "Samuels_GSC_AHS28_tr_mRNA",
                          "Samuels_GSC_AHS38_tr_mRNA",
                          "Douka_S2_UTD2_tr_mRNA",
                          "dmel_em02_tr_mRNA",
                          "dmel_em26_tr_mRNA",
                          "dmel_em612_tr_mRNA",
                          "dmel_em1224_tr_mRNA",
                          "dmel_larva_tr_mRNA",
                          "dmel_pupa_tr_mRNA",
                          "dmel_female_body_tr_mRNA",
                          "dmel_male_body_tr_mRNA",
                          "dmel_female_head_tr_mRNA",
                          "dmel_male_head_tr_mRNA",
                          "S2_tr_mRNA",
                          "Luo_S2_cell_normal_tr_mRNA",
                          "Luo_S2_cell_serum_starvation_tr_mRNA",
                          "Luo_S2_cell_ago2KD_tr_mRNA" ,
                          "Luo_S2_cell_ago2KD_serum_starvation_tr_mRNA",
                          "Luo_S2_cell_DMSO_tr_mRNA",
                          "Luo_S2_cell_dsNC_tr_mRNA",
                          "Luo_S2_cell_rapamycin_tr_mRNA",
                          "Luo_S2_cell_ssNC_tr_mRNA",
                          "Luo_S2_cell_T3_tr_mRNA",
                          "Luo_S2_cell_T6_tr_mRNA",
                          "Luo_S2_cell_T10_tr_mRNA",
                          "Dunn_em02_tr_mRNA")

k<-c("Kronja_mature_oocyte","Kronja_activated_egg",
     "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
     "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
     "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva","dmel_pupa",
     "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head","S2",
     "Luo_S2_cell_normal","Luo_S2_cell_serum_starvation",
     "Luo_S2_cell_ago2KD" ,"Luo_S2_cell_ago2KD_serum_starvation",
     "Luo_S2_cell_DMSO","Luo_S2_cell_dsNC",
     "Luo_S2_cell_rapamycin","Luo_S2_cell_ssNC",
     "Luo_S2_cell_T3","Luo_S2_cell_T6","Luo_S2_cell_T10",
     "Dunn_em02")
#libsize_merge$sample
k0<-c("Kronja_mature_oocyte","Kronja_activated_egg",
      "Samuels_GSC_AHS0","Samuels_GSC_AHS0.5","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
      "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
      "dmel_em_0_2h","dmel_em_2_6h","dmel_em_6_12h","dmel_em_12_24h","dmel_larva","dmel_pupa",
      "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head","S2",
      "Luo_S2_cell_normal","Luo_S2_cell_serum_starvation",
      "Luo_S2_cell_ago2KD" ,"Luo_S2_cell_ago2KD_serum_starvation",
      "Luo_S2_cell_DMSO","Luo_S2_cell_dsNC",
      "Luo_S2_cell_rapamycin","Luo_S2_cell_ssNC",
      "Luo_S2_cell_T3","Luo_S2_cell_T6","Luo_S2_cell_T10",
      "Dunn_em02")


#dmel tr mRNA

for(i in 1:length(k)){
  sample_id=paste0(k[i],"_tr_mRNA")
  total=libsize_merge[libsize_merge$sample==paste0(k0[i],"_mrna"),mapped]
  sample_rpkm_id=paste0(k[i],"_tr_mRNA_rpkm")
  dmel_tr_count[,sample_rpkm_id]=dmel_tr_count[,..sample_id]*1e9/total/dmel_tr_count$dmel_tr_len
}

gene=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/0.ann/r6.04.gene_tr",header=F)
colnames(gene)=c("gene","trid")
tr_mRNA_merge=merge(gene,dmel_tr_count,by="trid")

fwrite(tr_mRNA_merge,"tr_mRNA_merge_33sample.txt",sep='\t')
tr_mRNA_merge21sample=tr_mRNA_merge[,c("trid","gene","dmel_tr_len","Kronja_mature_oocyte_tr_mRNA",                    
                                       "Kronja_activated_egg_tr_mRNA","Samuels_GSC_AHS0_tr_mRNA",                        
                                       "Samuels_GSC_AHS05_tr_mRNA","Samuels_GSC_AHS5_tr_mRNA",                        
                                       "Samuels_GSC_AHS9_tr_mRNA","Samuels_GSC_AHS18_tr_mRNA" ,                     
                                       "Samuels_GSC_AHS28_tr_mRNA","Samuels_GSC_AHS38_tr_mRNA",                       
                                       "Douka_S2_UTD2_tr_mRNA","dmel_em02_tr_mRNA",                             
                                       "dmel_em26_tr_mRNA","dmel_em612_tr_mRNA",                              
                                       "dmel_em1224_tr_mRNA","dmel_larva_tr_mRNA",                              
                                       "dmel_pupa_tr_mRNA","dmel_female_body_tr_mRNA",                 
                                       "dmel_male_body_tr_mRNA" ,"dmel_female_head_tr_mRNA",                       
                                       "dmel_male_head_tr_mRNA",                         
                                       "Dunn_em02_tr_mRNA","Kronja_mature_oocyte_tr_mRNA_rpkm",  
                                       "Kronja_activated_egg_tr_mRNA_rpkm","Samuels_GSC_AHS0_tr_mRNA_rpkm",                   
                                       "Samuels_GSC_AHS05_tr_mRNA_rpkm","Samuels_GSC_AHS5_tr_mRNA_rpkm",                   
                                       "Samuels_GSC_AHS9_tr_mRNA_rpkm","Samuels_GSC_AHS18_tr_mRNA_rpkm",                  
                                       "Samuels_GSC_AHS28_tr_mRNA_rpkm" ,"Samuels_GSC_AHS38_tr_mRNA_rpkm",                  
                                       "Douka_S2_UTD2_tr_mRNA_rpkm","dmel_em02_tr_mRNA_rpkm",                          
                                       "dmel_em26_tr_mRNA_rpkm","dmel_em612_tr_mRNA_rpkm",                         
                                       "dmel_em1224_tr_mRNA_rpkm","dmel_larva_tr_mRNA_rpkm",                         
                                       "dmel_pupa_tr_mRNA_rpkm","dmel_female_body_tr_mRNA_rpkm",                   
                                       "dmel_male_body_tr_mRNA_rpkm","dmel_female_head_tr_mRNA_rpkm",                  
                                       "dmel_male_head_tr_mRNA_rpkm","Dunn_em02_tr_mRNA_rpkm")]
fwrite(tr_mRNA_merge21sample,"tr_mRNA_merge_21sample.txt",sep='\t')

#find the most abundant tr
tr_mRNA_merge=fread("tr_mRNA_merge_33sample.txt",header=T,sep='\t')
tr_mRNA_merge_rank=tr_mRNA_merge

for(i in 1:33){
  k<-c("Kronja_mature_oocyte","Kronja_activated_egg",
       "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
       "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
       "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva","dmel_pupa",
       "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head","S2",
       "Luo_S2_cell_normal","Luo_S2_cell_serum_starvation",
       "Luo_S2_cell_ago2KD" ,"Luo_S2_cell_ago2KD_serum_starvation",
       "Luo_S2_cell_DMSO","Luo_S2_cell_dsNC",
       "Luo_S2_cell_rapamycin","Luo_S2_cell_ssNC",
       "Luo_S2_cell_T3","Luo_S2_cell_T6","Luo_S2_cell_T10",
       "Dunn_em02")[i]
  selected_columns <- c("gene","trid","dmel_tr_len",paste0(k, "_tr_mRNA_rpkm"))
  x_dmel <- tr_mRNA_merge[,..selected_columns]
  x_dmel=x_dmel[order(-get(paste0(k, "_tr_mRNA_rpkm")),-get("dmel_tr_len")),]
  rank=paste0(k,"_tr_rank")
  x_dmel[,rank]=0
  for (j in unique(x_dmel$gene)){
    x_dmel[gene==j,rank]=c(1:length(x_dmel[gene==j,]$gene))
  }
  selected_columns2 <- c("trid",rank)
  x_dmel2=x_dmel[,..selected_columns2]
  tr_mRNA_merge_rank=merge(tr_mRNA_merge_rank,x_dmel2,by="trid",sort =F)
}

tr_mRNA_merge_rank_total=tr_mRNA_merge_rank
tr_mRNA_merge_rank_total[, dmel_mean_rpkm := rowMeans(.SD, na.rm = TRUE), .SDcols = 37:69]
selected_columns <- c("gene","trid","dmel_tr_len","dmel_mean_rpkm")
x_dmel <- tr_mRNA_merge_rank_total[,..selected_columns]
x_dmel=x_dmel[order(-get("dmel_mean_rpkm"),-get("dmel_tr_len")),]
rank=paste0("dmel_mean_tr_rank")
x_dmel[,rank]=0
for (j in unique(x_dmel$gene)){
  x_dmel[gene==j,rank]=c(1:length(x_dmel[gene==j,]$gene))
}
selected_columns2 <- c("trid",rank)
x_dmel2=x_dmel[,..selected_columns2]
tr_mRNA_merge_rank_total2=merge(tr_mRNA_merge_rank_total,x_dmel2,by="trid",sort =F)


fwrite(tr_mRNA_merge_rank_total2,"tr_mRNA_merge_rank_33sample.txt",sep='\t')


#### 21sample
#find the most abundant tr
tr_mRNA_merge=fread("tr_mRNA_merge_21sample.txt",header=T,sep='\t')
tr_mRNA_merge_rank=tr_mRNA_merge

for(i in 1:21){
  k<-c("Kronja_mature_oocyte","Kronja_activated_egg",
       "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
       "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
       "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva","dmel_pupa",
       "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
       "Dunn_em02")[i]
  selected_columns <- c("gene","trid","dmel_tr_len",paste0(k, "_tr_mRNA_rpkm"))
  x_dmel <- tr_mRNA_merge[,..selected_columns]
  x_dmel=x_dmel[order(-get(paste0(k, "_tr_mRNA_rpkm")),-get("dmel_tr_len")),]
  rank=paste0(k,"_tr_rank")
  x_dmel[,rank]=0
  for (j in unique(x_dmel$gene)){
    x_dmel[gene==j,rank]=c(1:length(x_dmel[gene==j,]$gene))
  }
  selected_columns2 <- c("trid",rank)
  x_dmel2=x_dmel[,..selected_columns2]
  tr_mRNA_merge_rank=merge(tr_mRNA_merge_rank,x_dmel2,by="trid",sort =F)
}

tr_mRNA_merge_rank_total=tr_mRNA_merge_rank
tr_mRNA_merge_rank_total[, dmel_mean_rpkm := rowMeans(.SD, na.rm = TRUE), .SDcols = 25:45]
selected_columns <- c("gene","trid","dmel_tr_len","dmel_mean_rpkm")
x_dmel <- tr_mRNA_merge_rank_total[,..selected_columns]
x_dmel=x_dmel[order(-get("dmel_mean_rpkm"),-get("dmel_tr_len")),]
rank=paste0("dmel_mean_tr_rank")
x_dmel[,rank]=0
for (j in unique(x_dmel$gene)){
  x_dmel[gene==j,rank]=c(1:length(x_dmel[gene==j,]$gene))
}
selected_columns2 <- c("trid",rank)
x_dmel2=x_dmel[,..selected_columns2]
tr_mRNA_merge_rank_total2=merge(tr_mRNA_merge_rank_total,x_dmel2,by="trid",sort =F)


fwrite(tr_mRNA_merge_rank_total2,"tr_mRNA_merge_rank_21sample.txt",sep='\t')

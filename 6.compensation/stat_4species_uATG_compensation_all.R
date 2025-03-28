library(data.table)
#new
#PacBioSim
uORF_matrix=fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",sep=',',header=T)
#number for fig1B
#remove duplicated
uORF_matrix_uniq=uORF_matrix[!duplicated(uORF_matrix, by = c("seqname", "gene_pos"))]

uORF_matrix_uniq[uORF_matrix_uniq$dm6==1&uORF_matrix_uniq$PacBioSim==1,]->PacBioSim_conserved #30975
uORF_matrix_uniq[uORF_matrix_uniq$dm6==1&uORF_matrix_uniq$PacBioSim==0,]->PacBioSim_mel #5590
uORF_matrix_uniq[uORF_matrix_uniq$dm6==0&uORF_matrix_uniq$PacBioSim==1,]->PacBioSim_sim  #4752

#class1*
uORF_matrix_uniq[uORF_matrix_uniq$dm6==1&uORF_matrix_uniq$PacBioSim==1,]->z #30975
length(unique(z$geneID)) #7053

#class6 lost in b2
uORF_matrix_uniq[uORF_matrix_uniq$dm6==1 & uORF_matrix_uniq$PacBioSim==0 & uORF_matrix_uniq$droYak3==1,]->z #1075
length(unique(z$geneID)) #829
#class5 gained in b2
uORF_matrix_uniq[uORF_matrix_uniq$dm6==0 & uORF_matrix_uniq$PacBioSim==1 & uORF_matrix_uniq$droYak3==0,]->z #3538
length(unique(z$geneID)) #2209
#class4 lost in b1
uORF_matrix_uniq[uORF_matrix_uniq$dm6==0 & uORF_matrix_uniq$PacBioSim==1 & uORF_matrix_uniq$droYak3==1,]->z #1214
length(unique(z$geneID)) #989
class4_list=unique(z$geneID)
class4_tr_list=unique(z$transcriptID)
#class3
uORF_matrix_uniq[uORF_matrix_uniq$dm6==1 & uORF_matrix_uniq$PacBioSim==0 & uORF_matrix_uniq$droYak3==0,]->z #4515
length(unique(z$geneID)) #2503
class3_list=unique(z$geneID)
class3_tr_list=unique(z$transcriptID)
#class2 lost in b3
uORF_matrix_uniq[uORF_matrix_uniq$dm6==0 & uORF_matrix_uniq$PacBioSim==0 & uORF_matrix_uniq$droYak3==1 & uORF_matrix_uniq$droAna3==1,]->z #372
length(unique(z$geneID)) #335
class2_list=unique(z$geneID)
class2_tr_list=unique(z$transcriptID)
#class1
uORF_matrix_uniq[uORF_matrix_uniq$dm6==1 & uORF_matrix_uniq$PacBioSim==1 & uORF_matrix_uniq$droYak3==0 & uORF_matrix_uniq$droAna3==0,]->z #5882
length(unique(z$geneID)) #3003
class1_list=unique(z$geneID)
class1_tr_list=unique(z$transcriptID)

#++++
uORF_matrix_uniq[uORF_matrix_uniq$dm6==1 & uORF_matrix_uniq$PacBioSim==1 & uORF_matrix_uniq$droYak3==1 & uORF_matrix_uniq$droAna3==1,]->z #10870
length(unique(z$geneID)) #3889
#+-++
uORF_matrix_uniq[uORF_matrix_uniq$dm6==1 & uORF_matrix_uniq$PacBioSim==1 & uORF_matrix_uniq$droYak3==0 & uORF_matrix_uniq$droAna3==1,]->z #702
length(unique(z$geneID)) #532
#-+++
uORF_matrix_uniq[uORF_matrix_uniq$dm6==1 & uORF_matrix_uniq$PacBioSim==1 & uORF_matrix_uniq$droYak3==1 & uORF_matrix_uniq$droAna3==0,]->z #13521
length(unique(z$geneID)) #4854



####gain-loss compensation for canonical transcripts
library(data.table)

canonical=fread("/data/dmel-all-r6.04.canonical.transcripts.txt",header=F,sep='\t')
#y=x[x$transcriptID %in% canonical$V3,]
#y[y$dm6==1&y$droSim1==1,]->z #17567
#y[y$dm6==1&y$droSim1==0,]->z #5274
#y[y$dm6==0&y$droSim1==1,]->z #2358

#new
#PacBioSim
uORF_matrix=fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",sep=',',header=T)
uORF_matrix_canonical=uORF_matrix[uORF_matrix$transcriptID %in% canonical$V3,]
uORF_matrix_canonical[uORF_matrix_canonical$dm6==1&uORF_matrix_canonical$PacBioSim==1,]->PacBioSim_conserved #18426
uORF_matrix_canonical[uORF_matrix_canonical$dm6==1&uORF_matrix_canonical$PacBioSim==0,]->PacBioSim_mel #2765
uORF_matrix_canonical[uORF_matrix_canonical$dm6==0&uORF_matrix_canonical$PacBioSim==1,]->PacBioSim_sim  #2450

#class1*
uORF_matrix_canonical[uORF_matrix_canonical$dm6==1&uORF_matrix_canonical$PacBioSim==1,]->z #18426
length(unique(z$geneID)) #6223

#class6 lost in b2
uORF_matrix_canonical[uORF_matrix_canonical$dm6==1 & uORF_matrix_canonical$PacBioSim==0 & uORF_matrix_canonical$droYak3==1,]->z #515
length(unique(z$geneID)) #451
#class5 gained in b2
uORF_matrix_canonical[uORF_matrix_canonical$dm6==0 & uORF_matrix_canonical$PacBioSim==1 & uORF_matrix_canonical$droYak3==0,]->z #1822
length(unique(z$geneID)) #1364
#class4 lost in b1
uORF_matrix_canonical[uORF_matrix_canonical$dm6==0 & uORF_matrix_canonical$PacBioSim==1 & uORF_matrix_canonical$droYak3==1,]->z #628
length(unique(z$geneID)) #553
class4_list=unique(z$geneID)
class4_tr_list=unique(z$transcriptID)
#class3
uORF_matrix_canonical[uORF_matrix_canonical$dm6==1 & uORF_matrix_canonical$PacBioSim==0 & uORF_matrix_canonical$droYak3==0,]->z #2250
length(unique(z$geneID)) #1604
class3_list=unique(z$geneID)
class3_tr_list=unique(z$transcriptID)
#class2 lost in b3
uORF_matrix_canonical[uORF_matrix_canonical$dm6==0 & uORF_matrix_canonical$PacBioSim==0 & uORF_matrix_canonical$droYak3==1 & uORF_matrix_canonical$droAna3==1,]->z #217
length(unique(z$geneID)) #203
class2_list=unique(z$geneID)
class2_tr_list=unique(z$transcriptID)
#class1
uORF_matrix_canonical[uORF_matrix_canonical$dm6==1 & uORF_matrix_canonical$PacBioSim==1 & uORF_matrix_canonical$droYak3==0 & uORF_matrix_canonical$droAna3==0,]->z #3009
length(unique(z$geneID)) #1970
class1_list=unique(z$geneID)
class1_tr_list=unique(z$transcriptID)

#++++
uORF_matrix_canonical[uORF_matrix_canonical$dm6==1 & uORF_matrix_canonical$PacBioSim==1 & uORF_matrix_canonical$droYak3==1 & uORF_matrix_canonical$droAna3==1,]->z #7046
length(unique(z$geneID)) #3332
#+-++
uORF_matrix_canonical[uORF_matrix_canonical$dm6==1 & uORF_matrix_canonical$PacBioSim==1 & uORF_matrix_canonical$droYak3==0 & uORF_matrix_canonical$droAna3==1,]->z #413
length(unique(z$geneID)) #329
#-+++
uORF_matrix_canonical[uORF_matrix_canonical$dm6==1 & uORF_matrix_canonical$PacBioSim==1 & uORF_matrix_canonical$droYak3==1 & uORF_matrix_canonical$droAna3==0,]->z #7958
length(unique(z$geneID)) #3894
#gain-loss compensation

#transcript
gain_loss_tr=intersect(class1_tr_list,class4_tr_list) #231
loss_gain_tr=intersect(class2_tr_list,class3_tr_list) #82
anc_gain_tr=class1_tr_list
anc_loss_tr=class2_tr_list
mel_gain_tr=class3_tr_list
mel_loss_tr=class4_tr_list
compensation_tr_list=list(anc_gain=anc_gain_tr,anc_loss=anc_loss_tr,mel_gain=mel_gain_tr,mel_loss=mel_loss_tr,gain_loss=gain_loss_tr,loss_gain=loss_gain_tr)
saveRDS(compensation_tr_list, file = "/results/uORF_gain_loss_trlist.rds")


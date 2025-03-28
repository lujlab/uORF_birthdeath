library(data.table)
shortntron_ATG_matrix=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/14.shortintron_ATG/ATG_matrix_dm6_PacBioSim_27species_shortintron_update.csv",sep=',',header=T)

#sp4_1111
sp4_1111=nrow(shortntron_ATG_matrix[shortntron_ATG_matrix$dm6==1 & shortntron_ATG_matrix $PacBioSim==1 &
                        shortntron_ATG_matrix$droYak3==1 & shortntron_ATG_matrix$droAna3==1,]) #136
sp4_1011=nrow(shortntron_ATG_matrix[shortntron_ATG_matrix$dm6==1 & shortntron_ATG_matrix $PacBioSim==1 &
                                 shortntron_ATG_matrix$droYak3==0 & shortntron_ATG_matrix$droAna3==1,]) #89
sp4_0111=nrow(shortntron_ATG_matrix[shortntron_ATG_matrix$dm6==1 & shortntron_ATG_matrix $PacBioSim==1 &
                                 shortntron_ATG_matrix$droYak3==1 & shortntron_ATG_matrix$droAna3==0,]) #1567
gain_in_b3=nrow(shortntron_ATG_matrix[shortntron_ATG_matrix$dm6==1 & shortntron_ATG_matrix $PacBioSim==1 &
                                 shortntron_ATG_matrix$droYak3==0 & shortntron_ATG_matrix$droAna3==0,]) #2422
lost_in_b3=nrow(shortntron_ATG_matrix[shortntron_ATG_matrix$dm6==0 & shortntron_ATG_matrix $PacBioSim==0 &
                                   shortntron_ATG_matrix$droYak3==1 & shortntron_ATG_matrix$droAna3==1,]) #84
gain_in_b1=nrow(shortntron_ATG_matrix[shortntron_ATG_matrix$dm6==1 & shortntron_ATG_matrix $PacBioSim==0 &
                                   shortntron_ATG_matrix$droYak3==0 ,]) #1836
lost_in_b1=nrow(shortntron_ATG_matrix[shortntron_ATG_matrix$dm6==0 & shortntron_ATG_matrix $PacBioSim==1 &
                                   shortntron_ATG_matrix$droYak3==1 ,]) #362
gain_in_b2=nrow(shortntron_ATG_matrix[shortntron_ATG_matrix$dm6==0 & shortntron_ATG_matrix $PacBioSim==1 &
                                   shortntron_ATG_matrix$droYak3==0 ,]) #1851
lost_in_b2=nrow(shortntron_ATG_matrix[shortntron_ATG_matrix$dm6==1 & shortntron_ATG_matrix $PacBioSim==0 &
                                   shortntron_ATG_matrix$droYak3==1 ,]) #394

total_dm6=nrow(shortntron_ATG_matrix[shortntron_ATG_matrix$dm6==1,]) #6444
total_sim=nrow(shortntron_ATG_matrix[shortntron_ATG_matrix$PacBioSim==1,]) #6427

f_shortintron=data.table(sp4_1111=sp4_1111,sp4_1011=sp4_1011,sp4_0111=sp4_0111,
                         gain_in_b3=gain_in_b3,lost_in_b3=lost_in_b3,
                         gain_in_b1=gain_in_b1,lost_in_b1=lost_in_b1,
                         gain_in_b2=gain_in_b2,lost_in_b2=lost_in_b2)
fwrite(f_shortintron,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/14.shortintron_ATG/shortintron_gainloss_number_4sp.txt")

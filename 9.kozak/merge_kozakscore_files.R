library(data.table)
dm6_kozak=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/3.kozak/dm6_addkozakscore_tmp",head=F)
dm6_kozak=dm6_kozak[,c(1,7,8,11,12,14,15)]
#dsim_kozak=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/3.kozak/dsim_uORF_nogap.bed.corrected.addkozakscore",head=F)
dsim_kozak=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/3.kozak/dsim_uORF_nogap.bed.corrected.addkozakscore2",head=F)
dsim_kozak=dsim_kozak[,c(4,5,8,9,11,12)]

colnames(dm6_kozak)=c("uorf_align","pos_withgap","dm6_nogap_pos","dm6","PacBioSim","kozak_seq_dm6","kozak_score_dm6")
colnames(dsim_kozak)=c("pos_withgap","dsim_nogap_pos","dm6","PacBioSim","kozak_seq_PacBioSim","kozak_score_PacBioSim")
m=merge(dm6_kozak,dsim_kozak,by=c("pos_withgap","dm6","PacBioSim"),all=TRUE)

#fwrite(m,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/3.kozak/uORF_kozakscore.txt",sep='\t')
fwrite(m,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/3.kozak/uORF_kozakscore.txt2",sep='\t')
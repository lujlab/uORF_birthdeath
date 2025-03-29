setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/4.count_rpkm_shift1")
library(data.table)
libsize1<- fread('/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.map3/map2genome/mrna_libsize',header=F)
colnames(libsize1)=c("lib","sample","mapped")
libsize2<- fread('/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.map3/lib_size/ribo_libsize',header=F)
colnames(libsize2)=c("lib","sample","mapped")
libsize=rbind(libsize1,libsize2)
libsize=libsize[,c("sample","mapped")]
libsize_merge=libsize[, .(mapped = sum(mapped)), by = .(sample)]
libsize_merge[libsize_merge$sample=="Douka_S2_mrna",]$sample="Douka_S2_UTD2_mrna"
libsize_merge[libsize_merge$sample=="Luo_S2_cell_normal_starvation_ribo",]$sample="Luo_S2_cell_normal_ribo"


#calculate Ribo rpkm and TE of potential ORF start codons NATGN, S2 harringtonine
dmel_uORF_count <- fread('/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/2.map3/all_potential_orf_start_plus1_cnt/all_potential_orf_start_plus1_cnt.tab')
dmel_uORF_count=dmel_uORF_count[,c("trid","uatg_tstart","uatg_tend","id","class","codon", "name","S2_harr.Psite_potential_orf_start_plus1_cnt.txt","S2.weighted_potential_orf_start_plus1_cnt.txt")]
colnames(dmel_uORF_count)=c("trid","uatg_tstart","uatg_tend","id","class","codon", "name",
                            "S2_NXXXN_Ribo","S2_NXXXN_mRNA")
dmel_uORF_count$len=5
dmel_uORF_count[dmel_uORF_count$uatg_tstart==0,]$len=4
total=libsize_merge[libsize_merge$sample=="S2_mrna",mapped]
dmel_uORF_count[,"S2_NXXXN_mRNA_rpkm"]=dmel_uORF_count[,"S2_NXXXN_mRNA"]*1e9/total/dmel_uORF_count$len
total=libsize_merge[libsize_merge$sample=="S2_ribo",mapped]
dmel_uORF_count[,"S2_NXXXN_Ribo_rpkm"]=dmel_uORF_count[,"S2_NXXXN_Ribo"]*1e9/total/dmel_uORF_count$len
dmel_uORF_count[,"S2_initiationsignal"]=dmel_uORF_count$S2_NXXXN_Ribo_rpkm/dmel_uORF_count$S2_NXXXN_mRNA_rpkm

fwrite(dmel_uORF_count,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak_add/4.count_rpkm_shift1/dmel_potential_uORF_start_signal_S2harr_rpkm.txt",sep='\t')

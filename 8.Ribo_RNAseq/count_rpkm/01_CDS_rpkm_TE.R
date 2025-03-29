
library(data.table)
library(plyr)
library(dplyr)
library(Rmisc)
#mel libsize
libsize_mel<- fread('/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/1.STAR_new2/mel/map2genome/mel.libsize',header=F)
colnames(libsize_mel)=c("lib","mapped")
libsize_mel2=data.table(lib=c("dmel_female_body_mrna","dmel_female_body_ribo","dmel_male_body_mrna","dmel_male_body_ribo"),
                        mapped=c(libsize_mel[libsize_mel$lib=="dmel_female_body1_mrna",]$mapped+libsize_mel[libsize_mel$lib=="dmel_female_body2_mrna",]$mapped,
                                 libsize_mel[libsize_mel$lib=="dmel_female_body1_ribo",]$mapped+libsize_mel[libsize_mel$lib=="dmel_female_body2_ribo",]$mapped,
                                 libsize_mel[libsize_mel$lib=="dmel_male_body1_mrna",]$mapped+libsize_mel[libsize_mel$lib=="dmel_male_body2_mrna",]$mapped,
                                 libsize_mel[libsize_mel$lib=="dmel_male_body1_ribo",]$mapped+libsize_mel[libsize_mel$lib=="dmel_male_body2_ribo",]$mapped))
libsize_mel=rbind(libsize_mel,libsize_mel2)

libsize_mel2=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/1.STAR_new2/mel2/lib_size/mel_ribo_lib",header=F)

colnames(libsize_mel2)=c("lib","mapped")
libsize_mel=rbind(libsize_mel,libsize_mel2)
libsize_mel2=data.table(lib=c("dmel_female_body","dmel_male_body"),
                        mapped=c(libsize_mel[libsize_mel$lib=="dmel_female_body1",]$mapped+libsize_mel[libsize_mel$lib=="dmel_female_body2",]$mapped,
                                 libsize_mel[libsize_mel$lib=="dmel_male_body1",]$mapped+libsize_mel[libsize_mel$lib=="dmel_male_body2",]$mapped))
libsize_mel=rbind(libsize_mel,libsize_mel2)
#dmel CDS
dmel_CDS_count <- fread('/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/1.STAR_new2/mel2/all_cds_cnt/all_cds_cnt.tab')
colnames(dmel_CDS_count)
dmel_CDS_count$`dmel_female_body.weighted`=dmel_CDS_count$`dmel_female_body1.weighted`+dmel_CDS_count$`dmel_female_body2.weighted`
dmel_CDS_count$`dmel_female_body.Psite`=dmel_CDS_count$`dmel_female_body1.Psite`+dmel_CDS_count$`dmel_female_body2.Psite`
dmel_CDS_count$`dmel_male_body.weighted`=dmel_CDS_count$`dmel_male_body1.weighted`+dmel_CDS_count$`dmel_male_body2.weighted`
dmel_CDS_count$`dmel_male_body.Psite`=dmel_CDS_count$`dmel_male_body1.Psite`+dmel_CDS_count$`dmel_male_body2.Psite`

dmel_CDS_count=dmel_CDS_count[,c("trid","cds_len","dmel_em_0_2h.weighted","dmel_em_0_2h.Psite","dmel_em_2_6h.weighted","dmel_em_2_6h.Psite",
           "dmel_em_6_12h.weighted", "dmel_em_6_12h.Psite","dmel_em_12_24h.weighted","dmel_em_12_24h.Psite","dmel_larva.weighted",
           "dmel_larva.Psite","dmel_pupa.weighted","dmel_pupa.Psite","dmel_female_body.weighted","dmel_female_body.Psite",
           "dmel_male_body.weighted","dmel_male_body.Psite","dmel_female_head.weighted","dmel_female_head.Psite",
           "dmel_male_head.weighted","dmel_male_head.Psite")]

colnames(dmel_CDS_count)=c( "trid","dmel_CDS_len","dmel_em02_CDS_mRNA","dmel_em02_CDS_Ribo","dmel_em26_CDS_mRNA","dmel_em26_CDS_Ribo",
                 "dmel_em612_CDS_mRNA","dmel_em612_CDS_Ribo","dmel_em1224_CDS_mRNA","dmel_em1224_CDS_Ribo","dmel_larva_CDS_mRNA","dmel_larva_CDS_Ribo",
                 "dmel_pupa_CDS_mRNA","dmel_pupa_CDS_Ribo","dmel_fem_CDS_mRNA","dmel_fem_CDS_Ribo","dmel_mal_CDS_mRNA","dmel_mal_CDS_Ribo",
                 "dmel_FH_CDS_mRNA","dmel_FH_CDS_Ribo","dmel_MH_CDS_mRNA","dmel_MH_CDS_Ribo")


k <- c("em02","em26","em612","em1224","larva","pupa","fem","mal","FH","MH")
k0 <- c("em_0_2h","em_2_6h","em_6_12h","em_12_24h","larva","pupa","female_body","male_body","female_head","male_head")


for(i in 1:length(k)){
  sample_id=paste0("dmel_",k[i],"_CDS_mRNA")
  total=libsize_mel[libsize_mel$lib==paste0("dmel_",k0[i],"_mrna"),mapped]
  sample_rpkm_id=paste0("dmel_",k[i],"_CDS_mRNA_rpkm")
  dmel_CDS_count[,sample_rpkm_id]=dmel_CDS_count[,..sample_id]*1e9/total/dmel_CDS_count$dmel_CDS_len
}
for(i in 1:length(k)){
  sample_id=paste0("dmel_",k[i],"_CDS_Ribo")
  total=libsize_mel[libsize_mel$lib==paste0("dmel_",k0[i]),mapped]
  sample_rpkm_id=paste0("dmel_",k[i],"_CDS_Ribo_rpkm")
  dmel_CDS_count[,sample_rpkm_id]=dmel_CDS_count[,..sample_id]*1e9/total/dmel_CDS_count$dmel_CDS_len
}


#TE
for(i in 1:length(k)){
  sample_TE_id=paste0("dmel_",k[i],"_TE")
  mRNA=paste0("dmel_",k[i],"_CDS_mRNA_rpkm")
  Ribo=paste0("dmel_",k[i],"_CDS_Ribo_rpkm")
  dmel_CDS_count[,sample_TE_id]=dmel_CDS_count[,..Ribo]/dmel_CDS_count[,..mRNA]
}

fwrite(dmel_CDS_count,"dmel_CDS_mRNA_ribo_merge.txt",sep='\t')

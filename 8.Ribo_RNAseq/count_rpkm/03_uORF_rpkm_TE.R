
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

#dmel uORF
dmel_uORF_count <- fread('/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/11.kozak/1.STAR_new2/mel2/all_uorf_cnt/all_uorf_cnt.tab')
colnames(dmel_uORF_count)
dmel_uORF_count$`dmel_female_body.weighted`=dmel_uORF_count$`dmel_female_body1.weighted`+dmel_uORF_count$`dmel_female_body2.weighted`
dmel_uORF_count$`dmel_female_body.Psite`=dmel_uORF_count$`dmel_female_body1.Psite`+dmel_uORF_count$`dmel_female_body2.Psite`
dmel_uORF_count$`dmel_male_body.weighted`=dmel_uORF_count$`dmel_male_body1.weighted`+dmel_uORF_count$`dmel_male_body2.weighted`
dmel_uORF_count$`dmel_male_body.Psite`=dmel_uORF_count$`dmel_male_body1.Psite`+dmel_uORF_count$`dmel_male_body2.Psite`

dmel_uORF_count=dmel_uORF_count[,c("trid","uorf_id","uorf_len","dmel_em_0_2h.weighted","dmel_em_0_2h.Psite","dmel_em_2_6h.weighted","dmel_em_2_6h.Psite",
                                 "dmel_em_6_12h.weighted", "dmel_em_6_12h.Psite","dmel_em_12_24h.weighted","dmel_em_12_24h.Psite","dmel_larva.weighted",
                                 "dmel_larva.Psite","dmel_pupa.weighted","dmel_pupa.Psite","dmel_female_body.weighted","dmel_female_body.Psite",
                                 "dmel_male_body.weighted","dmel_male_body.Psite","dmel_female_head.weighted","dmel_female_head.Psite",
                                 "dmel_male_head.weighted","dmel_male_head.Psite")]

colnames(dmel_uORF_count)=c( "trid","uorf_id","dmel_uORF_len","dmel_em02_uORF_mRNA","dmel_em02_uORF_Ribo","dmel_em26_uORF_mRNA","dmel_em26_uORF_Ribo",
                            "dmel_em612_uORF_mRNA","dmel_em612_uORF_Ribo","dmel_em1224_uORF_mRNA","dmel_em1224_uORF_Ribo","dmel_larva_uORF_mRNA","dmel_larva_uORF_Ribo",
                            "dmel_pupa_uORF_mRNA","dmel_pupa_uORF_Ribo","dmel_fem_uORF_mRNA","dmel_fem_uORF_Ribo","dmel_mal_uORF_mRNA","dmel_mal_uORF_Ribo",
                            "dmel_FH_uORF_mRNA","dmel_FH_uORF_Ribo","dmel_MH_uORF_mRNA","dmel_MH_uORF_Ribo")


k <- c("em02","em26","em612","em1224","larva","pupa","fem","mal","FH","MH")
k0 <- c("em_0_2h","em_2_6h","em_6_12h","em_12_24h","larva","pupa","female_body","male_body","female_head","male_head")
#dmel uORF mRNA

for(i in 1:length(k)){
  sample_id=paste0("dmel_",k[i],"_uORF_mRNA")
  total=libsize_mel[libsize_mel$lib==paste0("dmel_",k0[i],"_mrna"),mapped]
  sample_rpkm_id=paste0("dmel_",k[i],"_uORF_mRNA_rpkm")
  dmel_uORF_count[,sample_rpkm_id]=dmel_uORF_count[,..sample_id]*1e9/total/dmel_uORF_count$dmel_uORF_len
}
for(i in 1:length(k)){
  sample_id=paste0("dmel_",k[i],"_uORF_Ribo")
  total=libsize_mel[libsize_mel$lib==paste0("dmel_",k0[i]),mapped]
  sample_rpkm_id=paste0("dmel_",k[i],"_uORF_Ribo_rpkm")
  dmel_uORF_count[,sample_rpkm_id]=dmel_uORF_count[,..sample_id]*1e9/total/dmel_uORF_count$dmel_uORF_len
}



#TE
for(i in 1:length(k)){
  sample_TE_id=paste0("dmel_",k[i],"_TE")
  mRNA=paste0("dmel_",k[i],"_uORF_mRNA_rpkm")
  Ribo=paste0("dmel_",k[i],"_uORF_Ribo_rpkm")
  tmp=as.data.frame(dmel_uORF_count[,c(..mRNA,..Ribo)])
  colnames(tmp)=c("mRNA","Ribo")
  tmp[is.na(tmp)] <- 0
  tmp[tmp$Ribo>0 & tmp$mRNA==0,]$Ribo=tmp[tmp$Ribo>0 & tmp$mRNA==0,]$Ribo+0.1
  tmp[tmp$Ribo>0 & tmp$mRNA==0,]$mRNA=tmp[tmp$Ribo>0 & tmp$mRNA==0,]$mRNA+0.1
  dmel_uORF_count[,sample_TE_id]=tmp$Ribo/tmp$mRNA
  
}
fwrite(dmel_uORF_count,"dmel_uORF_mRNA_ribo_merge.txt",sep='\t')

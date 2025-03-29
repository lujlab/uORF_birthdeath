#RPF Psite related start sites

library(data.table)
library(plyr)
library(dplyr)
library(Rmisc)
#dmel potential startcodon
count <- fread('/data/rpkm/all_potential_start_cnt.tab',header=T)
#colnames(count)
f=fread("/data/rpkm/total_potential_start_sort.bed",header=F)
f=cbind(f,count)
f_ATG=f[f$V5=="pos0" & f$V6=="ATG",]

colnames(f_ATG)
f_ATG=f_ATG[,c("V4","Kronja_mature_oocyte.Psite_potential_start_cnt.txt",
               "Kronja_activated_egg.Psite_potential_start_cnt.txt",
               "Samuels_GSC_AHS0.Psite_potential_start_cnt.txt",
               "Samuels_GSC_AHS0.5.Psite_potential_start_cnt.txt",
               "Samuels_GSC_AHS5.Psite_potential_start_cnt.txt",
               "Samuels_GSC_AHS9.Psite_potential_start_cnt.txt",
               "Samuels_GSC_AHS18.Psite_potential_start_cnt.txt",
               "Samuels_GSC_AHS28.Psite_potential_start_cnt.txt",
               "Samuels_GSC_AHS38.Psite_potential_start_cnt.txt",
               "Douka_S2_UTD2.Psite_potential_start_cnt.txt",
               "dmel_em_0_2h.Psite_potential_start_cnt.txt",
               "dmel_em_2_6h.Psite_potential_start_cnt.txt",
               "dmel_em_6_12h.Psite_potential_start_cnt.txt",
               "dmel_em_12_24h.Psite_potential_start_cnt.txt",
               "dmel_larva.Psite_potential_start_cnt.txt",
               "dmel_pupa.Psite_potential_start_cnt.txt",
               "dmel_female_body.Psite_potential_start_cnt.txt",
               "dmel_male_body.Psite_potential_start_cnt.txt",
               "dmel_female_head.Psite_potential_start_cnt.txt",
               "dmel_male_head.Psite_potential_start_cnt.txt",
               "Dunn_em02.Psite_potential_start_cnt.txt")]

colnames(f_ATG)=c("uorf_align","Kronja_mature_oocyte","Kronja_activated_egg",
                  "Samuels_GSC_AHS0","Samuels_GSC_AHS0.5","Samuels_GSC_AHS5","Samuels_GSC_AHS9","Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38",
                  "Douka_S2_UTD2",
                  "dmel_em_0_2h",
                  "dmel_em_2_6h",
                  "dmel_em_6_12h",
                  "dmel_em_12_24h",
                  "dmel_larva",
                  "dmel_pupa",
                  "dmel_female_body",
                  "dmel_male_body",
                  "dmel_female_head",
                  "dmel_male_head",
                  "Dunn_em02")



uORF_matrix <- fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",sep=',',header=T)
colnames(uORF_matrix)[1]="uorf_align"
uORF_matrix$geno_pos=paste0(uORF_matrix$seqname,"_",uORF_matrix$gene_pos)

#merge(f_ATG,uORF_matrix["dm6"==1,],by="uorf_align")
f_ATG_pos=merge(f_ATG,uORF_matrix[dm6==1,c("uorf_align","geno_pos")],by="uorf_align")
f_ATG_pos_only=unique(f_ATG_pos[,2:ncol(f_ATG_pos)])
f_ATG_pos_only[duplicated(f_ATG_pos_only[,"geno_pos"])]

#f_mean=aggregate(n ~ geno_pos, data = f, FUN = max)
sample=c("Kronja_mature_oocyte","Kronja_activated_egg",
         "Samuels_GSC_AHS0","Samuels_GSC_AHS0.5","Samuels_GSC_AHS5","Samuels_GSC_AHS9","Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38",
         "Douka_S2_UTD2",
         "dmel_em_0_2h","dmel_em_2_6h","dmel_em_6_12h","dmel_em_12_24h","dmel_larva","dmel_pupa","dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
         "Dunn_em02")

uORF_matrix_unique=uORF_matrix[!duplicated(uORF_matrix[,"geno_pos"])]
uORF_matrix_unique_dm6=uORF_matrix_unique[dm6==1,]
uORF_matrix_unique_dm6_pos=uORF_matrix_unique_dm6[,c("seqname","gene_pos","geno_pos")]

out=uORF_matrix_unique_dm6_pos
for (i in c(1:21)){
  tmp=aggregate(get(sample[i]) ~ geno_pos, data = f_ATG_pos_only, FUN = max)
  colnames(tmp)=c("geno_pos",sample[i])
  out=merge(out,tmp,by="geno_pos",all.x = T)
}

out[is.na(out)]=0

out[, sum := rowSums(.SD, na.rm = TRUE), .SDcols = 4:24]

#fwrite(out,"./ATG_pos_Psite_cov_21sample.txt",sep='\t')

sample=c("Kronja_mature_oocyte","Kronja_activated_egg",
         "Samuels_GSC_AHS0","Samuels_GSC_AHS0.5","Samuels_GSC_AHS5","Samuels_GSC_AHS9","Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38",
         "Douka_S2_UTD2",
         "dmel_em_0_2h","dmel_em_2_6h","dmel_em_6_12h","dmel_em_12_24h","dmel_larva","dmel_pupa","dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
         "Dunn_em02","sum")

#filter translated uORFs
t=fread("/data/rpkm/translated_uORF_21sample.txt",sep='\t')
t$geno_pos=paste0(t$seqname,"_",t$gene_pos)
out2=out[out$geno_pos %in% t$geno_pos,]

num_ATG_withPsite=data.table(sample=sample,Psite1=0,Psite3=0)
for (i in c(1:22)){
  num_ATG_withPsite[i,2]=nrow(out2[get(sample[i])>=1,])
  num_ATG_withPsite[i,3]=nrow(out2[get(sample[i])>=3,])
  
}

fwrite(num_ATG_withPsite,"/results/tableS3_num_ATG_withPsite_translateduORF.txt",sep='\t')


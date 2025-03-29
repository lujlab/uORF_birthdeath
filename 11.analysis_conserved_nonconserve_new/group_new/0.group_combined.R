library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(Rmisc)
t=fread("/data/rpkm/group/translated_uORF_21sample.txt",sep='\t')
t$geno_pos=paste(t$seqname,t$gene_pos,sep="_")
t=t[t$uORF_TE>0.1,]

uORF_stat_3m=fread("/data/rpkm/group/uATG_genopos_class.txt",sep='\t')[,1:22]
colnames(uORF_stat_3m)=c("geno_pos","Kronja_mature_oocyte","Kronja_activated_egg","Samuels_GSC_AHS0","Samuels_GSC_AHS05",
                         "Samuels_GSC_AHS5","Samuels_GSC_AHS9","Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38",
                         "Douka_S2_UTD2","dmel_em02","dmel_em26","dmel_em612","dmel_em1224",
                         "dmel_larva" ,"dmel_pupa","dmel_female_body","dmel_male_body","dmel_female_head",
                         "dmel_male_head","Dunn_em02")
sample=c("Kronja_mature_oocyte","Kronja_activated_egg","Samuels_GSC_AHS0","Samuels_GSC_AHS05",
         "Samuels_GSC_AHS5","Samuels_GSC_AHS9","Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38",
         "Douka_S2_UTD2","dmel_em02","dmel_em26","dmel_em612","dmel_em1224",
         "dmel_larva" ,"dmel_pupa","dmel_female_body","dmel_male_body","dmel_female_head",
         "dmel_male_head","Dunn_em02")



####merge TE>0.1 and identified by threemethod
uORF_stat_combined=data.table(matrix(ncol=22,nrow=36565))
colnames(uORF_stat_combined) <- c("geno_pos","Kronja_mature_oocyte","Kronja_activated_egg",
                         "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
                         "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
                         "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva",
                         "dmel_pupa","dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
                         "Dunn_em02")
uORF_stat_combined$geno_pos <- uORF_stat_3m$geno_pos
uORF_stat_combined[is.na(uORF_stat_combined)]=0

sample=c("Kronja_mature_oocyte","Kronja_activated_egg",
          "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
          "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
          "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva",
          "dmel_pupa","dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
          "Dunn_em02")


for (i in c(1:21)){
  k=sample[i]
  selected_columns=c("geno_pos",k)
  temp1=uORF_stat_3m[,..selected_columns] #identified by 3 methods
  colnames(temp1)=c("geno_pos","3m")
  temp1_pos=temp1[`3m`==1, ]$geno_pos
  index <- uORF_stat_combined$geno_pos %in% temp1_pos
  uORF_stat_combined[index, (k) := 1]
  temp2_pos=t[sample==k,]$geno_pos #identified by 3 methods
  index <- uORF_stat_combined$geno_pos %in% temp2_pos
  uORF_stat_combined[index, (k) := 1]
}

uORF_stat_combined[, sum := rowSums(.SD, na.rm = TRUE), .SDcols = c("Kronja_mature_oocyte","Kronja_activated_egg",
                                                                    "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
                                                                    "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
                                                                    "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva",
                                                                    "dmel_pupa","dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
                                                                    "Dunn_em02")]
#nrow(uORF_stat_combined[sum>0,]) #30037
#30037/36565=82.15%


class=data.frame(table(uORF_stat_combined$sum))
colnames(class)=c("sample","uORF")
class$sample=as.numeric(as.character(class$sample))
class$uORF=as.numeric(class$uORF)
#class2=class[class$sample>0,]

pdf("/results/figS17.pdf",width=6,height=4)
ggplot(class)+
  geom_bar(aes(x=sample,y=uORF),stat="identity",width=0.8)+
  theme_classic() #sample_uORF_group_combined
dev.off()

sum(class[class$sample>14,]$uORF) #8837
sum(class[class$sample>=6 & class$sample<=14,]$uORF) #8585
sum(class[class$sample>=2 & class$sample<=5,]$uORF) #8595
sum(class[class$sample<2,]$uORF) #10548
sum(class[class$sample>=1,]$uORF) #30037

uORF_stat_combined$group="group0"
uORF_stat_combined[sum>=2&sum<=5]$group="group1"
uORF_stat_combined[sum>5&sum<=14]$group="group2"
uORF_stat_combined[sum>14]$group="group3"


fwrite(uORF_stat_combined,"/results/uATG_genopos_TEgroup_combined.txt",sep='\t')


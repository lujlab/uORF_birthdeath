
library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(Rmisc)

##### CDS #####
file=fread("/data/periodicity/framefile_21sample",header=F)
colnames(file)=c("filename","sample0")
file$frame0=0
file$frame1=0
file$frame2=0
# most abundant tr
dm6_most_tr=fread("/data/rpkm/tr_mRNA_merge_rank_21sample.txt",header=T,sep='\t')

s=data.table(sample0<-c("Kronja_mature_oocyte","Kronja_activated_egg",
                        "Samuels_GSC_AHS0","Samuels_GSC_AHS0.5","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
                        "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
                        "dmel_em_0_2h","dmel_em_2_6h","dmel_em_6_12h","dmel_em_12_24h","dmel_larva","dmel_pupa",
                        "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
                        "Dunn_em02"),
             sample<-c("Kronja_mature_oocyte","Kronja_activated_egg",
                       "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
                       "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
                       "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva","dmel_pupa",
                       "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
                       "Dunn_em02"))
colnames(s)=c("sample0","sample")
file=merge(file,s,by="sample0")

for (i in c(1:21)){
  j=file$filename[i]
  k=file$sample[i]
  f=fread(paste0("/data/periodicity/",j),header=F)
  f_most=f[f$V1 %in% dm6_most_tr[get(paste0(k,"_tr_rank"))==1,]$trid,]
  
  file[filename==j,]$frame0=sum(f_most[f_most$V9==1,]$V4)
  file[filename==j,]$frame1=sum(f_most[f_most$V9==2,]$V4)
  file[filename==j,]$frame2=sum(f_most[f_most$V9==0,]$V4)
}


file2=file[,c("filename","sample0","frame0","frame1","frame2")]

file2$total=file2$frame0+file2$frame1+file2$frame2
file2$frame0_proportion=file2$frame0/file2$total
file2$frame1_proportion=file2$frame1/file2$total
file2$frame2_proportion=file2$frame2/file2$total

f4=file2[,c("sample0","frame0_proportion","frame1_proportion", "frame2_proportion")]
colnames(f4)=c("sample","frame0","frame1", "frame2")
f5 <- melt(f4, id.vars = "sample", measure.vars = c("frame0","frame1", "frame2"))
samplelevels=c("Kronja_activated_egg","Kronja_mature_oocyte","Samuels_GSC_AHS0","Samuels_GSC_AHS0.5","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
               "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
               "dmel_em_0_2h","dmel_em_2_6h","dmel_em_6_12h","dmel_em_12_24h","dmel_larva","dmel_pupa",
               "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
               "Dunn_em02")
f5$sample=factor(f5$sample,levels=samplelevels)                              
colnames(f5)=c("sample","frame","proportion")


pdf("/results/FigS1A_periodicity_CDS.pdf",width=8,height=6)
ggplot(data=f5,aes(x=sample,y=proportion,fill=frame))+
  geom_bar(stat="identity",position= position_dodge(0.8))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #periodicity_CDS_21sample_mosttr.pdf 8*6
dev.off()

##### CDS #####
file=fread("/data/periodicity/framefile_uORF_21sample",header=F)
colnames(file)=c("filename","sample0")
file$frame0=0
file$frame1=0
file$frame2=0
# most abundant tr
dm6_most_tr=fread("/data/rpkm/tr_mRNA_merge_rank_21sample.txt",header=T,sep='\t')

s=data.table(sample0<-c("Kronja_mature_oocyte","Kronja_activated_egg",
                        "Samuels_GSC_AHS0","Samuels_GSC_AHS0.5","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
                        "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
                        "dmel_em_0_2h","dmel_em_2_6h","dmel_em_6_12h","dmel_em_12_24h","dmel_larva","dmel_pupa",
                        "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
                        "Dunn_em02"),
             sample<-c("Kronja_mature_oocyte","Kronja_activated_egg",
                       "Samuels_GSC_AHS0","Samuels_GSC_AHS05","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
                       "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
                       "dmel_em02","dmel_em26","dmel_em612","dmel_em1224","dmel_larva","dmel_pupa",
                       "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
                       "Dunn_em02"))
colnames(s)=c("sample0","sample")
file=merge(file,s,by="sample0")

for (i in c(1:21)){
  j=file$filename[i]
  k=file$sample[i]
  f=fread(paste0("/data/periodicity/",j),header=F)
  f_most=f[f$V1 %in% dm6_most_tr[get(paste0(k,"_tr_rank"))==1,]$trid,]
  
  file[filename==j,]$frame0=sum(f_most[f_most$V9==1,]$V4)
  file[filename==j,]$frame1=sum(f_most[f_most$V9==2,]$V4)
  file[filename==j,]$frame2=sum(f_most[f_most$V9==0,]$V4)
}


file2=file[,c("filename","sample0","frame0","frame1","frame2")]

file2$total=file2$frame0+file2$frame1+file2$frame2
file2$frame0_proportion=file2$frame0/file2$total
file2$frame1_proportion=file2$frame1/file2$total
file2$frame2_proportion=file2$frame2/file2$total

f4=file2[,c("sample0","frame0_proportion","frame1_proportion", "frame2_proportion")]
colnames(f4)=c("sample","frame0","frame1", "frame2")
f5 <- melt(f4, id.vars = "sample", measure.vars = c("frame0","frame1", "frame2"))
samplelevels=c("Kronja_activated_egg","Kronja_mature_oocyte","Samuels_GSC_AHS0","Samuels_GSC_AHS0.5","Samuels_GSC_AHS5","Samuels_GSC_AHS9",
               "Samuels_GSC_AHS18","Samuels_GSC_AHS28","Samuels_GSC_AHS38","Douka_S2_UTD2",
               "dmel_em_0_2h","dmel_em_2_6h","dmel_em_6_12h","dmel_em_12_24h","dmel_larva","dmel_pupa",
               "dmel_female_body","dmel_male_body","dmel_female_head","dmel_male_head",
               "Dunn_em02")
f5$sample=factor(f5$sample,levels=samplelevels)                              
colnames(f5)=c("sample","frame","proportion")


pdf("/results/FigS1B_periodicity_uORF.pdf",width=8,height=6)
ggplot(data=f5,aes(x=sample,y=proportion,fill=frame))+
  geom_bar(stat="identity",position= position_dodge(0.8))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #periodicity_uORF_21sample_mosttr.pdf 8*6
dev.off()
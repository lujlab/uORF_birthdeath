library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)


fisher.test(matrix(c(4530696-36565,36565,393875-6444,6444),nrow=2)) #p-value < 2.2e-16
fisher.test(matrix(c(4530696-4515,4515,393875-1836,1836),nrow=2)) #p-value < 2.2e-16
fisher.test(matrix(c(4530696-1214,1214,393875-362,362),nrow=2)) #p-value < 2.2e-16
fisher.test(matrix(c(4515,1214,1836,362),nrow=2)) #p-value = 0.000002021



#CG content of intron (71243  + 62726)/393875=0.3401308
#CG content of intron (1006813 + 970267)/4530696=0.4363745

##### for5UTR
#ATG number for gene
f1<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv")
f1_unique=f1[!duplicated(f1[,c("seqname","gene_pos")]),]
#f1_unique_dm6=f1_unique[dm6==1,]
f1_unique$gain=0
f1_unique$loss=0
f1_unique[f1_unique$dm6==1&f1_unique$PacBioSim==0&f1_unique$droYak3==0,]$gain=1
f1_unique[f1_unique$dm6==0&f1_unique$PacBioSim==1&f1_unique$droYak3==1,]$loss=1
#sum(f1_unique$gain) #4515
#sum(f1_unique$loss) # 1214
uATG_number=data.frame(table(f1_unique[dm6==1,]$geneID))
colnames(uATG_number)=c("FBid","uATG_number")
uATG_gain=data.frame(table(f1_unique[gain==1,]$geneID))
colnames(uATG_gain)=c("FBid","uATG_gain")
uATG_loss=data.frame(table(f1_unique[loss==1,]$geneID))
colnames(uATG_loss)=c("FBid","uATG_loss")
f2=merge(uATG_number,uATG_gain,all=T,by="FBid")
f2=merge(f2,uATG_loss,all=T,by="FBid")
f2[is.na(f2)]=0

#dinuc stat
f3=fread("/data/tr_len_27sp/dm6_merge_gene_stat_rmcds.txt",header=T)

f4=merge(f3,f2,by="FBid",all=T)
f4[is.na(f4)]=0
f4$type="5UTR"
f_5utr=f4

##### for shortintron
f5<-fread("/data/shortintron/ATG_matrix_dm6_PacBioSim_27species_shortintron_update.csv")
f5_unique=f5[!duplicated(f5[,c("seqname","gene_pos")]),]
f5_unique_dm6=f5_unique[dm6==1,] #6444

f5_unique$gain=0
f5_unique$loss=0
f5_unique[f5_unique$dm6==1&f5_unique$PacBioSim==0&f5_unique$droYak3==0,]$gain=1
f5_unique[f5_unique$dm6==0&f5_unique$PacBioSim==1&f5_unique$droYak3==1,]$loss=1
#sum(f5_unique$gain) #1836
#sum(f5_unique$loss) #360

uATG_number=data.frame(table(f5_unique[dm6==1,]$intronid))
colnames(uATG_number)=c("FBid","uATG_number")
uATG_gain=data.frame(table(f5_unique[gain==1,]$intronid))
colnames(uATG_gain)=c("FBid","uATG_gain")
uATG_loss=data.frame(table(f5_unique[loss==1,]$intronid))
colnames(uATG_loss)=c("FBid","uATG_loss")

f6=merge(uATG_number,uATG_gain,all=T,by="FBid")
f6=merge(f6,uATG_loss,all=T,by="FBid")
f6[is.na(f6)]=0
#dinuc stat
#grep -E "dm6|#" dm6_PacBioSim_27species_shortintron.dinuc.txt|sed 's/dm6.//g' >dm6_shortintron.dinuc.txt
f7=fread("/data/shortintron/dm6_shortintron.dinuc.txt",header=T)
colnames(f7)[1]="FBid"
f8=merge(f7,f6,by="FBid",all=T)
f8[is.na(f8)]=0
f8$type="shortintron"


f_shortintron=f8


#bin

f_5utr$GC_content=(f_5utr$C+f_5utr$G)/(f_5utr$C+f_5utr$G+f_5utr$A+f_5utr$T)
f_5utr$num_by_len=f_5utr$uATG_number/f_5utr$len
f_5utr <- f_5utr[order(f_5utr$GC_content), ]
f_5utr$gain_by_len=f_5utr$uATG_gain/f_5utr$len
f_5utr$loss_by_len=f_5utr$uATG_loss/f_5utr$len
f_5utr$gain_by_loss=f_5utr$uATG_gain/f_5utr$uATG_loss


f_shortintron$GC_content=(f_shortintron$C+f_shortintron$G)/(f_shortintron$C+f_shortintron$G+f_shortintron$A+f_shortintron$T)
f_shortintron$num_by_len=f_shortintron$uATG_number/f_shortintron$len
f_shortintron <- f_shortintron[order(f_shortintron$GC_content), ]
f_shortintron$gain_by_len=f_shortintron$uATG_gain/f_shortintron$len
f_shortintron$loss_by_len=f_shortintron$uATG_loss/f_shortintron$len
f_shortintron$gain_by_loss=f_shortintron$uATG_gain/f_shortintron$uATG_loss

###8 bins for shortintron, same breaks for 5utr
n_bins <- 8
f_shortintron$Bin <- cut(1:nrow(f_shortintron), 
                         breaks = n_bins, 
                         labels = FALSE)
bin_sum_shortintron <- f_shortintron %>%
  group_by(Bin) %>%
  summarise(
    Sum_C = sum(C),
    Sum_G = sum(G),
    Sum_A = sum(A),
    Sum_T = sum(T),
    Sum_len = sum(len),
    Sum_gain = sum(uATG_gain),
    Sum_loss = sum(uATG_loss),
    Sum_uATG = sum(uATG_number),
    max_GC=max(GC_content),
    min_GC=min(GC_content)
  )
bin_sum_shortintron$GC_content=(bin_sum_shortintron$Sum_C+bin_sum_shortintron$Sum_G)/(bin_sum_shortintron$Sum_C+bin_sum_shortintron$Sum_G+bin_sum_shortintron$Sum_A+bin_sum_shortintron$Sum_T)

bin_breaks <- c(bin_sum_shortintron$min_GC,1)
#bin for 5UTR, according to GC_content
f_5utr$Bin <- cut(f_5utr$GC_content, 
                  breaks = bin_breaks, 
                  labels = FALSE, 
                  include.lowest = TRUE, 
                  right = TRUE)
bin_sum_5utr <- f_5utr %>%
  group_by(Bin) %>%
  summarise(
    Sum_C = sum(C),
    Sum_G = sum(G),
    Sum_A = sum(A),
    Sum_T = sum(T),
    Sum_len = sum(len),
    Sum_gain = sum(uATG_gain),
    Sum_loss = sum(uATG_loss),
    Sum_uATG = sum(uATG_number),
    max_GC=max(GC_content),
    min_GC=min(GC_content)
  )
bin_sum_5utr$GC_content=(bin_sum_5utr$Sum_C+bin_sum_5utr$Sum_G)/(bin_sum_5utr$Sum_C+bin_sum_5utr$Sum_G+bin_sum_5utr$Sum_A+bin_sum_5utr$Sum_T)

pdf("/results/figS17a_uATG_GC_8bin.pdf",width=6,height=4)
ggplot()+
  geom_point(data=bin_sum_shortintron,aes(x=GC_content,y=Sum_uATG/Sum_len),color="red",alpha=0.5)+
  geom_point(data=bin_sum_5utr,aes(x=GC_content,y=Sum_uATG/Sum_len),color="blue",alpha=0.5)+
  geom_smooth(method = 'lm', se = FALSE) +
  geom_vline(xintercept = bin_breaks,color="grey",linetype="dashed")+
  theme_classic()+
  #coord_cartesian(xlim=c(-5,5),ylim=c(-5,0))+
  labs(x="Mean_GCcontent",y="Sum_uATG/Sum_len") #uATG_GC_8bin
dev.off()

#bin_sum_5utr and bin_sum_shortintron should be in the same order for the paired test
wilcox.test(bin_sum_5utr$Sum_uATG/bin_sum_5utr$Sum_len,bin_sum_shortintron$Sum_uATG/bin_sum_shortintron$Sum_len,paired = T) #p-value = 0.007813

pdf("/results/figS17b_uATG_gain_GC_8bin.pdf",width=6,height=4)
ggplot()+
  geom_point(data=bin_sum_shortintron,aes(x=GC_content,y=Sum_gain/Sum_len),color="red",alpha=0.5)+
  geom_point(data=bin_sum_5utr,aes(x=GC_content,y=Sum_gain/Sum_len),color="blue",alpha=0.5)+
  geom_smooth(method = 'lm', se = FALSE) +
  geom_vline(xintercept = bin_breaks,color="grey",linetype="dashed")+
  theme_classic()+
  #coord_cartesian(xlim=c(-5,5),ylim=c(-5,0))+
  labs(x="Mean_GCcontent",y="Sum_gain/Sum_len") #uATG_gain_GC_8bin
dev.off()

#bin_sum_5utr and bin_sum_shortintron should be in the same order for the paired test
wilcox.test(bin_sum_5utr$Sum_gain/bin_sum_5utr$Sum_len,bin_sum_shortintron$Sum_gain/bin_sum_shortintron$Sum_len,paired = T) #p-value = 0.007813

pdf("/results/figS17c_uATGloss_GC_8bin.pdf",width=6,height=4)
ggplot()+
  geom_point(data=bin_sum_shortintron,aes(x=GC_content,y=Sum_loss/Sum_len),color="red",alpha=0.5)+
  geom_point(data=bin_sum_5utr,aes(x=GC_content,y=Sum_loss/Sum_len),color="blue",alpha=0.5)+
  geom_smooth(method = 'lm', se = FALSE) +
  geom_vline(xintercept = bin_breaks,color="grey",linetype="dashed")+
  theme_classic()+
  #coord_cartesian(xlim=c(-5,5),ylim=c(-5,0))+
  labs(x="Mean_GCcontent",y="Sum_loss/Sum_len") #uATGloss_GC_8bin
dev.off()

#bin_sum_5utr and bin_sum_shortintron should be in the same order for the paired test
wilcox.test(bin_sum_5utr$Sum_loss/bin_sum_5utr$Sum_len,bin_sum_shortintron$Sum_loss/bin_sum_shortintron$Sum_len,paired = T) #p-value = 0.007813

pdf("/results/figS17d_uATGgainloss_GC_8bin.pdf",width=6,height=4)
ggplot()+
  geom_point(data=bin_sum_shortintron,aes(x=GC_content,y=Sum_gain/Sum_loss),color="red",alpha=0.5)+
  geom_point(data=bin_sum_5utr,aes(x=GC_content,y=Sum_gain/Sum_loss),color="blue",alpha=0.5)+
  geom_smooth(method = 'lm', se = FALSE) +
  geom_vline(xintercept = bin_breaks,color="grey",linetype="dashed")+
  theme_classic()+
  #coord_cartesian(xlim=c(-5,5),ylim=c(-5,0))+
  labs(x="Mean_GCcontent",y="Sum_gain/Sum_loss") #uATGgainloss_GC_8bin
dev.off()

#bin_sum_5utr and bin_sum_shortintron should be in the same order for the paired test
wilcox.test(bin_sum_5utr$Sum_gain/bin_sum_5utr$Sum_loss,bin_sum_shortintron$Sum_gain/bin_sum_shortintron$Sum_loss,paired = T) #p-value = 0.1484
wilcox.test((bin_sum_5utr$Sum_gain/bin_sum_5utr$Sum_loss)[c(1,3:8)],(bin_sum_shortintron$Sum_gain/bin_sum_shortintron$Sum_loss)[c(1,3:8)],paired = T) #p-value = 0.01563




#compare 5utr and intron in each bins, p value
p_table=data.table(Bin=1:8,uATG_number_p=NA,uATGgain_p=NA,uATGloss_p=NA,uATGgainlossratio_p=NA)
for (i in c(1:8)){
  p_table$uATG_number_p[i]=fisher.test(matrix(c(bin_sum_5utr[bin_sum_5utr$Bin==i,]$Sum_uATG,
                                                bin_sum_5utr[bin_sum_5utr$Bin==i,]$Sum_len,
                                                bin_sum_shortintron[bin_sum_shortintron$Bin==i,]$Sum_uATG,
                                                bin_sum_shortintron[bin_sum_shortintron$Bin==i,]$Sum_len),nrow=2),alternative = "less")$p.value
  p_table$uATGgain_p[i]=fisher.test(matrix(c(bin_sum_5utr[bin_sum_5utr$Bin==i,]$Sum_gain,
                                                bin_sum_5utr[bin_sum_5utr$Bin==i,]$Sum_len,
                                                bin_sum_shortintron[bin_sum_shortintron$Bin==i,]$Sum_gain,
                                                bin_sum_shortintron[bin_sum_shortintron$Bin==i,]$Sum_len),nrow=2),alternative = "less")$p.value
  p_table$uATGloss_p[i]=fisher.test(matrix(c(bin_sum_5utr[bin_sum_5utr$Bin==i,]$Sum_loss,
                                                bin_sum_5utr[bin_sum_5utr$Bin==i,]$Sum_len,
                                                bin_sum_shortintron[bin_sum_shortintron$Bin==i,]$Sum_loss,
                                                bin_sum_shortintron[bin_sum_shortintron$Bin==i,]$Sum_len),nrow=2),alternative = "less")$p.value
  p_table$uATGgainlossratio_p[i]=fisher.test(matrix(c(bin_sum_5utr[bin_sum_5utr$Bin==i,]$Sum_gain,
                                                bin_sum_5utr[bin_sum_5utr$Bin==i,]$Sum_loss,
                                                bin_sum_shortintron[bin_sum_shortintron$Bin==i,]$Sum_gain,
                                                bin_sum_shortintron[bin_sum_shortintron$Bin==i,]$Sum_loss),nrow=2),alternative = "less")$p.value
  
}
#fisher's method to combine
t1=-2 * sum(log(p_table$uATG_number_p)) 
pchisq(t1,df=16,lower.tail = FALSE) #p=0,p<4.940656e-324
# pchisq(1567,df=16,lower.tail = FALSE)
#[1] 0
# pchisq(1566,df=16,lower.tail = FALSE)
#[1] 4.940656e-324

t1=-2 * sum(log(p_table$uATGgain_p))
pchisq(t1,df=16,lower.tail = FALSE) #p=4.519243e-310
t1=-2 * sum(log(p_table$uATGloss_p))
pchisq(t1,df=16,lower.tail = FALSE) #p=2.398629e-46
t1=-2 * sum(log(p_table$uATGgainlossratio_p))
pchisq(t1,df=16,lower.tail = FALSE) #p=0.003054655

#test_statistic <- -2 * sum(log(row))
#combined_p_value <- pchisq(test_statistic, df = 2 * length(row), lower.tail = FALSE)
















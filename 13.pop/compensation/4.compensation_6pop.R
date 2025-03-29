library(data.table)
library(ggplot2)
library(ggpubr)
library(UpSetR)
CN=fread("/data/6pop_compensation/CN_gain_loss_p.txt",header=T)
FR=fread("/data/6pop_compensation/FR_gain_loss_p.txt",header=T)
RAL=fread("/data/6pop_compensation/RAL_gain_loss_p.txt",header=T)
EF=fread("/data/6pop_compensation/EF_gain_loss_p.txt",header=T)
SD=fread("/data/6pop_compensation/SD_gain_loss_p.txt",header=T)
ZI=fread("/data/6pop_compensation/ZI_gain_loss_p.txt",header=T)


#daf5
CN=CN[a_freq>0.05 & a_freq<0.95 & b_freq>0.05 & b_freq<0.95]
FR=FR[a_freq>0.05 & a_freq<0.95 & b_freq>0.05 & b_freq<0.95]
RAL=RAL[a_freq>0.05 & a_freq<0.95 & b_freq>0.05 & b_freq<0.95]
EF=EF[a_freq>0.05 & a_freq<0.95 & b_freq>0.05 & b_freq<0.95]
SD=SD[a_freq>0.05 & a_freq<0.95 & b_freq>0.05 & b_freq<0.95]
ZI=ZI[a_freq>0.05 & a_freq<0.95 & b_freq>0.05 & b_freq<0.95]

sig_matrix=data.table(pop=c("CN","FR","RAL","EF","SD","ZI"),
                      `D>0_FDR>0.05`=rep(0,6),
                      `D>0_FDR<0.05`=rep(0,6),
                      `D<0_FDR>0.05`=rep(0,6),
                      `D<0_FDR<0.05`=rep(0,6))

#rm dup
CN2=CN[!duplicated(CN[,c("gain_chr","gain_pos","loss_chr","loss_pos")]),] #73 pair
CN2$q=p.adjust(CN2$p,method = "BH")
CN2$ab_OEratio=CN2$ab_new/CN2$ab_Expected
CN2$Ab_OEratio=CN2$Ab_new/CN2$Ab_Expected
CN2$aB_OEratio=CN2$aB_new/CN2$aB_Expected
CN2$AB_OEratio=CN2$AB_new/CN2$AB_Expected
#D=P(AB)-P(A)*P(B)
#r2=D*D/(P(A)P(a)P(B)P(b))
CN2$D=CN2$ab_new/(CN2$AB_new+CN2$aB_new+CN2$Ab_new+CN2$ab_new)-CN2$a_freq*CN2$b_freq
CN2$r2=CN2$D*CN2$D/(CN2$A_freq*CN2$B_freq*CN2$a_freq*CN2$b_freq)
sig_matrix[pop=="CN"]$`D>0_FDR>0.05`=nrow(CN2[CN2$D>0 & CN2$q>0.05])
sig_matrix[pop=="CN"]$`D>0_FDR<0.05`=nrow(CN2[CN2$D>0 & CN2$q<0.05])
sig_matrix[pop=="CN"]$`D<0_FDR>0.05`=nrow(CN2[CN2$D<0 & CN2$q>0.05])
sig_matrix[pop=="CN"]$`D<0_FDR<0.05`=nrow(CN2[CN2$D<0 & CN2$q<0.05])
#fwrite(CN2,"/results/CN_gain_loss_maf5_uniq_q.txt",sep='\t')

FR2=FR[!duplicated(FR[,c("gain_chr","gain_pos","loss_chr","loss_pos")]),] #123 pair
FR2$q=p.adjust(FR2$p,method = "BH")
FR2$ab_OEratio=FR2$ab_new/FR2$ab_Expected
FR2$Ab_OEratio=FR2$Ab_new/FR2$Ab_Expected
FR2$aB_OEratio=FR2$aB_new/FR2$aB_Expected
FR2$AB_OEratio=FR2$AB_new/FR2$AB_Expected
#D=P(AB)-P(A)*P(B)
#r2=D*D/(P(A)P(a)P(B)P(b))
FR2$D=FR2$ab_new/(FR2$AB_new+FR2$aB_new+FR2$Ab_new+FR2$ab_new)-FR2$a_freq*FR2$b_freq
FR2$r2=FR2$D*FR2$D/(FR2$A_freq*FR2$B_freq*FR2$a_freq*FR2$b_freq)
sig_matrix[pop=="FR"]$`D>0_FDR>0.05`=nrow(FR2[FR2$D>0 & FR2$q>0.05])
sig_matrix[pop=="FR"]$`D>0_FDR<0.05`=nrow(FR2[FR2$D>0 & FR2$q<0.05])
sig_matrix[pop=="FR"]$`D<0_FDR>0.05`=nrow(FR2[FR2$D<0 & FR2$q>0.05])
sig_matrix[pop=="FR"]$`D<0_FDR<0.05`=nrow(FR2[FR2$D<0 & FR2$q<0.05])
#fwrite(FR2,"/results/FR_gain_loss_maf5_uniq_q.txt",sep='\t')

RAL2=RAL[!duplicated(RAL[,c("gain_chr","gain_pos","loss_chr","loss_pos")]),] #211 pair
RAL2$q=p.adjust(RAL2$p,method = "BH")
RAL2$ab_OEratio=RAL2$ab_new/RAL2$ab_Expected
RAL2$Ab_OEratio=RAL2$Ab_new/RAL2$Ab_Expected
RAL2$aB_OEratio=RAL2$aB_new/RAL2$aB_Expected
RAL2$AB_OEratio=RAL2$AB_new/RAL2$AB_Expected
#D=P(AB)-P(A)*P(B)
#r2=D*D/(P(A)P(a)P(B)P(b))
RAL2$D=RAL2$ab_new/(RAL2$AB_new+RAL2$aB_new+RAL2$Ab_new+RAL2$ab_new)-RAL2$a_freq*RAL2$b_freq
RAL2$r2=RAL2$D*RAL2$D/(RAL2$A_freq*RAL2$B_freq*RAL2$a_freq*RAL2$b_freq)
sig_matrix[pop=="RAL"]$`D>0_FDR>0.05`=nrow(RAL2[RAL2$D>0 & RAL2$q>0.05])
sig_matrix[pop=="RAL"]$`D>0_FDR<0.05`=nrow(RAL2[RAL2$D>0 & RAL2$q<0.05])
sig_matrix[pop=="RAL"]$`D<0_FDR>0.05`=nrow(RAL2[RAL2$D<0 & RAL2$q>0.05])
sig_matrix[pop=="RAL"]$`D<0_FDR<0.05`=nrow(RAL2[RAL2$D<0 & RAL2$q<0.05])
#fwrite(RAL2,"/results/RAL_gain_loss_maf5_uniq_q.txt",sep='\t')


EF2=EF[!duplicated(EF[,c("gain_chr","gain_pos","loss_chr","loss_pos")]),] #227 pair
EF2$q=p.adjust(EF2$p,method = "BH")
EF2$ab_OEratio=EF2$ab_new/EF2$ab_Expected
EF2$Ab_OEratio=EF2$Ab_new/EF2$Ab_Expected
EF2$aB_OEratio=EF2$aB_new/EF2$aB_Expected
EF2$AB_OEratio=EF2$AB_new/EF2$AB_Expected
#D=P(AB)-P(A)*P(B)
#r2=D*D/(P(A)P(a)P(B)P(b))
EF2$D=EF2$ab_new/(EF2$AB_new+EF2$aB_new+EF2$Ab_new+EF2$ab_new)-EF2$a_freq*EF2$b_freq
EF2$r2=EF2$D*EF2$D/(EF2$A_freq*EF2$B_freq*EF2$a_freq*EF2$b_freq)
sig_matrix[pop=="EF"]$`D>0_FDR>0.05`=nrow(EF2[EF2$D>0 & EF2$q>0.05])
sig_matrix[pop=="EF"]$`D>0_FDR<0.05`=nrow(EF2[EF2$D>0 & EF2$q<0.05])
sig_matrix[pop=="EF"]$`D<0_FDR>0.05`=nrow(EF2[EF2$D<0 & EF2$q>0.05])
sig_matrix[pop=="EF"]$`D<0_FDR<0.05`=nrow(EF2[EF2$D<0 & EF2$q<0.05])
#fwrite(EF2,"/results/EF_gain_loss_maf5_uniq_q.txt",sep='\t')


SD2=SD[!duplicated(SD[,c("gain_chr","gain_pos","loss_chr","loss_pos")]),] #242 pair
SD2$q=p.adjust(SD2$p,method = "BH")
SD2$ab_OEratio=SD2$ab_new/SD2$ab_Expected
SD2$Ab_OEratio=SD2$Ab_new/SD2$Ab_Expected
SD2$aB_OEratio=SD2$aB_new/SD2$aB_Expected
SD2$AB_OEratio=SD2$AB_new/SD2$AB_Expected
#D=P(AB)-P(A)*P(B)
#r2=D*D/(P(A)P(a)P(B)P(b))
SD2$D=SD2$ab_new/(SD2$AB_new+SD2$aB_new+SD2$Ab_new+SD2$ab_new)-SD2$a_freq*SD2$b_freq
SD2$r2=SD2$D*SD2$D/(SD2$A_freq*SD2$B_freq*SD2$a_freq*SD2$b_freq)
sig_matrix[pop=="SD"]$`D>0_FDR>0.05`=nrow(SD2[SD2$D>0 & SD2$q>0.05])
sig_matrix[pop=="SD"]$`D>0_FDR<0.05`=nrow(SD2[SD2$D>0 & SD2$q<0.05])
sig_matrix[pop=="SD"]$`D<0_FDR>0.05`=nrow(SD2[SD2$D<0 & SD2$q>0.05])
sig_matrix[pop=="SD"]$`D<0_FDR<0.05`=nrow(SD2[SD2$D<0 & SD2$q<0.05])
#fwrite(SD2,"/results/SD_gain_loss_maf5_uniq_q.txt",sep='\t')



ZI2=ZI[!duplicated(ZI[,c("gain_chr","gain_pos","loss_chr","loss_pos")]),] 
ZI2$q=p.adjust(ZI2$p,method = "BH")
ZI2$ab_OEratio=ZI2$ab_new/ZI2$ab_Expected
ZI2$Ab_OEratio=ZI2$Ab_new/ZI2$Ab_Expected
ZI2$aB_OEratio=ZI2$aB_new/ZI2$aB_Expected
ZI2$AB_OEratio=ZI2$AB_new/ZI2$AB_Expected
#D=P(AB)-P(A)*P(B)
#r2=D*D/(P(A)P(a)P(B)P(b))
ZI2$D=ZI2$ab_new/(ZI2$AB_new+ZI2$aB_new+ZI2$Ab_new+ZI2$ab_new)-ZI2$a_freq*ZI2$b_freq
ZI2$r2=ZI2$D*ZI2$D/(ZI2$A_freq*ZI2$B_freq*ZI2$a_freq*ZI2$b_freq)
sig_matrix[pop=="ZI"]$`D>0_FDR>0.05`=nrow(ZI2[ZI2$D>0 & ZI2$q>0.05])
sig_matrix[pop=="ZI"]$`D>0_FDR<0.05`=nrow(ZI2[ZI2$D>0 & ZI2$q<0.05])
sig_matrix[pop=="ZI"]$`D<0_FDR>0.05`=nrow(ZI2[ZI2$D<0 & ZI2$q>0.05])
sig_matrix[pop=="ZI"]$`D<0_FDR<0.05`=nrow(ZI2[ZI2$D<0 & ZI2$q<0.05])
#fwrite(ZI2,"/results/ZI_gain_loss_maf5_uniq_q.txt",sep='\t')

###compare positively correlated and negatively correlated
#sig_matrix$oddsratio=NA
#sig_matrix$p=NA

#for (i in c(1:6)){
#  
#  compare=fisher.test(matrix(c(sig_matrix$`D>0_FDR<0.05`[i],
#                               sig_matrix$`D>0_FDR>0.05`[i],
#                               sig_matrix$`D<0_FDR<0.05`[i],
#                               sig_matrix$`D<0_FDR>0.05`[i]),nrow=2))
#  sig_matrix$oddsratio[i]=compare$estimate
#  sig_matrix$p[i]=compare$p
#}

#t1=-2 * sum(log(sig_matrix$p)) 
#pchisq(t1,df=16,lower.tail = FALSE) #0.5677312



###upset
CN3=CN2[CN2$q<0.05 & CN2$D>0]
FR3=FR2[FR2$q<0.05 & FR2$D>0]
RAL3=RAL2[RAL2$q<0.05 & RAL2$D>0]
EF3=EF2[EF2$q<0.05 & EF2$D>0]
SD3=SD2[SD2$q<0.05 & SD2$D>0]
ZI3=ZI2[ZI2$q<0.05 & ZI2$D>0]

CN3$id=paste(CN3$gain_chr,CN3$gain_pos,CN3$loss_chr,CN3$loss_pos,sep="_")
FR3$id=paste(FR3$gain_chr,FR3$gain_pos,FR3$loss_chr,FR3$loss_pos,sep="_")
RAL3$id=paste(RAL3$gain_chr,RAL3$gain_pos,RAL3$loss_chr,RAL3$loss_pos,sep="_")
EF3$id=paste(EF3$gain_chr,EF3$gain_pos,EF3$loss_chr,EF3$loss_pos,sep="_")
SD3$id=paste(SD3$gain_chr,SD3$gain_pos,SD3$loss_chr,SD3$loss_pos,sep="_")
ZI3$id=paste(ZI3$gain_chr,ZI3$gain_pos,ZI3$loss_chr,ZI3$loss_pos,sep="_")

CN3$pop="CN"
FR3$pop="FR"
RAL3$pop="RAL"
EF3$pop="EF"
SD3$pop="SD"
ZI3$pop="ZI"

f=rbind(CN3[,c("id","pop")],FR3[,c("id","pop")],RAL3[,c("id","pop")],EF3[,c("id","pop")],SD3[,c("id","pop")],ZI3[,c("id","pop")])
f2=data.table(dcast(data=f, id~pop,fun.aggregate=length))
f3=f2[,c("CN","FR","RAL","EF","SD","ZI")]
#f3[f3 > 0] <- 1

pdf("/results/figS24a_upset_6pop_uATG_compensation.pdf",width=10,height=6)
upset(f3, sets=rev(c("CN","FR","RAL","EF","SD","ZI")),nsets = 6, nintersects = 40,keep.order = T,mb.ratio = c(0.5, 0.5)) #upset_6pop_uATG_compensation
dev.off()

#f2[f2$CN==1 & f2$FR==1 & f2$RAL==1 & f2$EF==1 &f2$SD==1 &f2$ZI==1]
#id    CN    EF    FR   RAL    SD    ZI
#1: 2R_10303372_2R_10303346     1     1     1     1     1     1
#2: 2R_20979257_2R_20979196     1     1     1     1     1     1
#3:   X_14866339_X_14866345     1     1     1     1     1     1

#2R_10303372_2R_10303346
CN_eg1=data.table(t(CN3[CN3$id=="2R_10303372_2R_10303346",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(CN_eg1)=c("count")
CN_eg1$allele=c("AB","Ab","aB","ab")
CN_eg1$freq=CN_eg1$count/(sum(CN_eg1$count))
CN_eg1$allele=factor(CN_eg1$allele,levels=c("AB","Ab","aB","ab"))
CN_eg1$pop="CN"

FR_eg1=data.table(t(FR3[FR3$id=="2R_10303372_2R_10303346",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(FR_eg1)=c("count")
FR_eg1$allele=c("AB","Ab","aB","ab")
FR_eg1$freq=FR_eg1$count/(sum(FR_eg1$count))
FR_eg1$allele=factor(FR_eg1$allele,levels=c("AB","Ab","aB","ab"))
FR_eg1$pop="FR"

RAL_eg1=data.table(t(RAL3[RAL3$id=="2R_10303372_2R_10303346",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(RAL_eg1)=c("count")
RAL_eg1$allele=c("AB","Ab","aB","ab")
RAL_eg1$freq=RAL_eg1$count/(sum(RAL_eg1$count))
RAL_eg1$allele=factor(RAL_eg1$allele,levels=c("AB","Ab","aB","ab"))
RAL_eg1$pop="RAL"

EF_eg1=data.table(t(EF3[EF3$id=="2R_10303372_2R_10303346",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(EF_eg1)=c("count")
EF_eg1$allele=c("AB","Ab","aB","ab")
EF_eg1$freq=EF_eg1$count/(sum(EF_eg1$count))
EF_eg1$allele=factor(EF_eg1$allele,levels=c("AB","Ab","aB","ab"))
EF_eg1$pop="EF"

SD_eg1=data.table(t(SD3[SD3$id=="2R_10303372_2R_10303346",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(SD_eg1)=c("count")
SD_eg1$allele=c("AB","Ab","aB","ab")
SD_eg1$freq=SD_eg1$count/(sum(SD_eg1$count))
SD_eg1$allele=factor(SD_eg1$allele,levels=c("AB","Ab","aB","ab"))
SD_eg1$pop="SD"

ZI_eg1=data.table(t(ZI3[ZI3$id=="2R_10303372_2R_10303346",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(ZI_eg1)=c("count")
ZI_eg1$allele=c("AB","Ab","aB","ab")
ZI_eg1$freq=ZI_eg1$count/(sum(ZI_eg1$count))
ZI_eg1$allele=factor(ZI_eg1$allele,levels=c("AB","Ab","aB","ab"))
ZI_eg1$pop="ZI"

x1=rbind(CN_eg1,FR_eg1,RAL_eg1,EF_eg1,SD_eg1,ZI_eg1)
x1$allele=factor(x1$allele,levels=rev(c("AB","Ab","aB","ab")))
x1$pop=factor(x1$pop,levels=c("CN","FR","RAL","EF","SD","ZI"))

p1=ggplot(x1,aes(x=pop,y=freq,fill=allele))+
  geom_bar(position="stack",stat="identity",width=0.8)+
  scale_fill_manual(values=c("#CD3129","#F8DE91","#9CBEE1","#5873B1"))+
  theme_classic() #2R_10303372_2R_10303346_FBtr0299828_freq


#2R_20979257_2R_20979196
CN_eg1=data.table(t(CN3[CN3$id=="2R_20979257_2R_20979196",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(CN_eg1)=c("count")
CN_eg1$allele=c("AB","Ab","aB","ab")
CN_eg1$freq=CN_eg1$count/(sum(CN_eg1$count))
CN_eg1$allele=factor(CN_eg1$allele,levels=c("AB","Ab","aB","ab"))
CN_eg1$pop="CN"

FR_eg1=data.table(t(FR3[FR3$id=="2R_20979257_2R_20979196",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(FR_eg1)=c("count")
FR_eg1$allele=c("AB","Ab","aB","ab")
FR_eg1$freq=FR_eg1$count/(sum(FR_eg1$count))
FR_eg1$allele=factor(FR_eg1$allele,levels=c("AB","Ab","aB","ab"))
FR_eg1$pop="FR"

RAL_eg1=data.table(t(RAL3[RAL3$id=="2R_20979257_2R_20979196",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(RAL_eg1)=c("count")
RAL_eg1$allele=c("AB","Ab","aB","ab")
RAL_eg1$freq=RAL_eg1$count/(sum(RAL_eg1$count))
RAL_eg1$allele=factor(RAL_eg1$allele,levels=c("AB","Ab","aB","ab"))
RAL_eg1$pop="RAL"

EF_eg1=data.table(t(EF3[EF3$id=="2R_20979257_2R_20979196",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(EF_eg1)=c("count")
EF_eg1$allele=c("AB","Ab","aB","ab")
EF_eg1$freq=EF_eg1$count/(sum(EF_eg1$count))
EF_eg1$allele=factor(EF_eg1$allele,levels=c("AB","Ab","aB","ab"))
EF_eg1$pop="EF"

SD_eg1=data.table(t(SD3[SD3$id=="2R_20979257_2R_20979196",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(SD_eg1)=c("count")
SD_eg1$allele=c("AB","Ab","aB","ab")
SD_eg1$freq=SD_eg1$count/(sum(SD_eg1$count))
SD_eg1$allele=factor(SD_eg1$allele,levels=c("AB","Ab","aB","ab"))
SD_eg1$pop="SD"

ZI_eg1=data.table(t(ZI3[ZI3$id=="2R_20979257_2R_20979196",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(ZI_eg1)=c("count")
ZI_eg1$allele=c("AB","Ab","aB","ab")
ZI_eg1$freq=ZI_eg1$count/(sum(ZI_eg1$count))
ZI_eg1$allele=factor(ZI_eg1$allele,levels=c("AB","Ab","aB","ab"))
ZI_eg1$pop="ZI"

x1=rbind(CN_eg1,FR_eg1,RAL_eg1,EF_eg1,SD_eg1,ZI_eg1)
x1$allele=factor(x1$allele,levels=rev(c("AB","Ab","aB","ab")))
x1$pop=factor(x1$pop,levels=c("CN","FR","RAL","EF","SD","ZI"))

p2=ggplot(x1,aes(x=pop,y=freq,fill=allele))+
  geom_bar(position="stack",stat="identity",width=0.8)+
  scale_fill_manual(values=c("#CD3129","#F8DE91","#9CBEE1","#5873B1"))+
  theme_classic() #2R_20979257_2R_20979196_FBtr0345666_freq



#X_14866339_X_14866345
CN_eg1=data.table(t(CN3[CN3$id=="X_14866339_X_14866345",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(CN_eg1)=c("count")
CN_eg1$allele=c("AB","Ab","aB","ab")
CN_eg1$freq=CN_eg1$count/(sum(CN_eg1$count))
CN_eg1$allele=factor(CN_eg1$allele,levels=c("AB","Ab","aB","ab"))
CN_eg1$pop="CN"

FR_eg1=data.table(t(FR3[FR3$id=="X_14866339_X_14866345",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(FR_eg1)=c("count")
FR_eg1$allele=c("AB","Ab","aB","ab")
FR_eg1$freq=FR_eg1$count/(sum(FR_eg1$count))
FR_eg1$allele=factor(FR_eg1$allele,levels=c("AB","Ab","aB","ab"))
FR_eg1$pop="FR"

RAL_eg1=data.table(t(RAL3[RAL3$id=="X_14866339_X_14866345",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(RAL_eg1)=c("count")
RAL_eg1$allele=c("AB","Ab","aB","ab")
RAL_eg1$freq=RAL_eg1$count/(sum(RAL_eg1$count))
RAL_eg1$allele=factor(RAL_eg1$allele,levels=c("AB","Ab","aB","ab"))
RAL_eg1$pop="RAL"

EF_eg1=data.table(t(EF3[EF3$id=="X_14866339_X_14866345",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(EF_eg1)=c("count")
EF_eg1$allele=c("AB","Ab","aB","ab")
EF_eg1$freq=EF_eg1$count/(sum(EF_eg1$count))
EF_eg1$allele=factor(EF_eg1$allele,levels=c("AB","Ab","aB","ab"))
EF_eg1$pop="EF"

SD_eg1=data.table(t(SD3[SD3$id=="X_14866339_X_14866345",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(SD_eg1)=c("count")
SD_eg1$allele=c("AB","Ab","aB","ab")
SD_eg1$freq=SD_eg1$count/(sum(SD_eg1$count))
SD_eg1$allele=factor(SD_eg1$allele,levels=c("AB","Ab","aB","ab"))
SD_eg1$pop="SD"

ZI_eg1=data.table(t(ZI3[ZI3$id=="X_14866339_X_14866345",c("AB_new","Ab_new","aB_new","ab_new")]))
colnames(ZI_eg1)=c("count")
ZI_eg1$allele=c("AB","Ab","aB","ab")
ZI_eg1$freq=ZI_eg1$count/(sum(ZI_eg1$count))
ZI_eg1$allele=factor(ZI_eg1$allele,levels=c("AB","Ab","aB","ab"))
ZI_eg1$pop="ZI"

x1=rbind(CN_eg1,FR_eg1,RAL_eg1,EF_eg1,SD_eg1,ZI_eg1)
x1$allele=factor(x1$allele,levels=rev(c("AB","Ab","aB","ab")))
x1$pop=factor(x1$pop,levels=c("CN","FR","RAL","EF","SD","ZI"))

p3=ggplot(x1,aes(x=pop,y=freq,fill=allele))+
  geom_bar(position="stack",stat="identity",width=0.8)+
  scale_fill_manual(values=c("#CD3129","#F8DE91","#9CBEE1","#5873B1"))+
  theme_classic() #X_14866339_X_14866345_FBtr0073983_freq


pdf("/results/figS24bcd_pop_compensation_eg.pdf",width=6,height=3)
ggarrange(p1,p2,p3,ncol=3,common.legend = T)
dev.off()
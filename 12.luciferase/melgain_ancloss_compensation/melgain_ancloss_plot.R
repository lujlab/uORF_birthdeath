library(data.table)
library(ggplot2)
library(ggpubr)
f=fread("/code/12.luciferase/melgain_ancloss_compensation/melgain_ancloss_protein.txt",header=T)
f$ratio=f$renilla/f$firefly
f$mut=factor(f$mut,levels=c("mel","sim","yak"))

#qPCR result
q=fread("/code/12.luciferase/melgain_ancloss_compensation/melgain_ancloss_mRNA.txt",header=T)

#q_ct=aggregate(q$ct,by=list(q$ID,q$`5UTR`,q$mut,q$treatment,q$target),mean)
q_ct=data.table(aggregate(q$ct,by=list(q$ID,q$`5UTR`,q$mut,q$treatment,q$target),mean))
colnames(q_ct)=c("ID","5UTR","mut","treatment","target","ct_mean")
#delta ct
q_ct=dcast(q_ct,ID+`5UTR`+mut+treatment~target,value.var="ct_mean")
q_ct$delta_ct=q_ct$Rluc-q_ct$fluc

wilcox.test(q_ct[q_ct$`5UTR`=="Art4"&q_ct$mut=="mel"&q_ct$treatment=="H2O",]$delta_ct,
            q_ct[q_ct$`5UTR`=="Art4"&q_ct$mut=="sim"&q_ct$treatment=="H2O",]$delta_ct)
wilcox.test(q_ct[q_ct$`5UTR`=="CG9933"&q_ct$mut=="mel"&q_ct$treatment=="H2O",]$delta_ct,
            q_ct[q_ct$`5UTR`=="CG9933"&q_ct$mut=="sim"&q_ct$treatment=="H2O",]$delta_ct)
wilcox.test(q_ct[q_ct$`5UTR`=="nrv3"&q_ct$mut=="mel"&q_ct$treatment=="H2O",]$delta_ct,
            q_ct[q_ct$`5UTR`=="nrv3"&q_ct$mut=="sim"&q_ct$treatment=="H2O",]$delta_ct)

wilcox.test(q_ct[q_ct$`5UTR`=="Art4"&q_ct$mut=="sim"&q_ct$treatment=="H2O",]$delta_ct,
            q_ct[q_ct$`5UTR`=="Art4"&q_ct$mut=="yak"&q_ct$treatment=="H2O",]$delta_ct)
wilcox.test(q_ct[q_ct$`5UTR`=="CG9933"&q_ct$mut=="sim"&q_ct$treatment=="H2O",]$delta_ct,
            q_ct[q_ct$`5UTR`=="CG9933"&q_ct$mut=="yak"&q_ct$treatment=="H2O",]$delta_ct)
wilcox.test(q_ct[q_ct$`5UTR`=="nrv3"&q_ct$mut=="sim"&q_ct$treatment=="H2O",]$delta_ct,
            q_ct[q_ct$`5UTR`=="nrv3"&q_ct$mut=="yak"&q_ct$treatment=="H2O",]$delta_ct)


q_ct2=aggregate(q_ct$delta_ct,by=list(q_ct$`5UTR`,q_ct$mut,q_ct$treatment),median)
colnames(q_ct2)=c("5UTR","mut","treatment","delta_ct")
q_ct2$`2^delta_ct`=2^(-q_ct2$delta_ct)

#Correction of expression
f2=merge(f,q_ct2[,c("5UTR","mut","treatment","2^delta_ct")])

f2$ratio2=f2$ratio/f2$`2^delta_ct`
f2$mut=factor(f2$mut,levels=c("mel","sim","yak"))

mean_r=aggregate(f2$ratio2,by=list(f2$`5UTR`,f2$mut,f2$treatment),mean)
colnames(mean_r)=c("5UTR","mut","treatment","n")
mean_r=mean_r[mean_r$mut=="mel",]
f3=merge(f2,mean_r[,c("5UTR","treatment","n")],by=c("5UTR","treatment"))
f3$ratio3=f3$ratio2/f3$n


i="Art4"
p1=ggplot(f3[f3$`5UTR`==i,],aes(x=treatment,y=ratio3,fill=mut))+
  geom_bar(stat="summary",fun="mean",width=0.5,color="black",position=position_dodge())+
  stat_summary(fun.data="mean_sd",geom="errorbar",color="black",width=0.2,position=position_dodge(width = 0.5))+
  geom_jitter(aes(x=treatment,y=ratio3,fill=mut),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),size=0.6)+
  #ylim(0,1.3)+
  ylab("relative translational efficiency")+
  xlab(i)+
  theme_classic()
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="H2O",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="H2O",]$ratio3) #p-value = 0.008658
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="H2O",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="H2O",]$ratio3) #p-value = 0.09307
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="H2O",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="H2O",]$ratio3) #p-value = 0.002165

wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="NaAsO2",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="NaAsO2",]$ratio3) #p-value = 0.04113
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="NaAsO2",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="NaAsO2",]$ratio3) #p-value = 0.1797
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="NaAsO2",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="NaAsO2",]$ratio3) #p-value = 0.02597

mean(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="H2O",]$ratio3)
mean(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="H2O",]$ratio3) #1.276158
mean(f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="H2O",]$ratio3) #1.370875
mean(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="NaAsO2",]$ratio3) #1.192864
mean(f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="NaAsO2",]$ratio3) #1.234105

i="CG9933"
p2=ggplot(f3[f3$`5UTR`==i,],aes(x=treatment,y=ratio3,fill=mut))+
  geom_bar(stat="summary",fun="mean",width=0.5,color="black",position=position_dodge())+
  stat_summary(fun.data="mean_sd",geom="errorbar",color="black",width=0.2,position=position_dodge(width = 0.5))+
  geom_jitter(aes(x=treatment,y=ratio3,fill=mut),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),size=0.6)+
  ylab("relative translational efficiency")+
  xlab(i)+
  #ylim(0,1.2)+
  theme_classic()
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="H2O",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="H2O",]$ratio3) #p-value = 0.008658
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="H2O",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="H2O",]$ratio3) #p-value = 0.09307
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="H2O",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="H2O",]$ratio3) #p-value = 0.002165

wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="NaAsO2",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="NaAsO2",]$ratio3) #p-value = 0.04113
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="NaAsO2",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="NaAsO2",]$ratio3) #p-value = 0.1797
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="NaAsO2",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="NaAsO2",]$ratio3) #p-value = 0.02597

mean(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="H2O",]$ratio3)
mean(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="H2O",]$ratio3) #1.1183
mean(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="NaAsO2",]$ratio3) #1
mean(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="NaAsO2",]$ratio3) #1.442753

i="nrv3"
p3=ggplot(f3[f3$`5UTR`==i,],aes(x=treatment,y=ratio3,fill=mut))+
  geom_bar(stat="summary",fun="mean",width=0.5,color="black",position=position_dodge())+
  stat_summary(fun.data="mean_sd",geom="errorbar",color="black",width=0.2,position=position_dodge(width = 0.5))+
  geom_jitter(aes(x=treatment,y=ratio3,fill=mut),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),size=0.6)+
  ylab("relative translational efficiency")+
  xlab(i)+
  #coord_cartesian(ylim = c(0.6, 1.5))+
  theme_classic()
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="H2O",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="H2O",]$ratio3) #p-value = 0.008658
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="H2O",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="H2O",]$ratio3) #p-value = 0.09307
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="H2O",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="H2O",]$ratio3) #p-value = 0.002165

wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="NaAsO2",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="NaAsO2",]$ratio3) #p-value =  0.004329
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="NaAsO2",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="NaAsO2",]$ratio3) #p-value = 0.04113
wilcox.test(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="NaAsO2",]$ratio3,
            f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="NaAsO2",]$ratio3) #p-value = 0.02597

mean(f3[f3$`5UTR`==i&f3$mut=="mel"&f3$treatment=="H2O",]$ratio3)
mean(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="H2O",]$ratio3) #1.889986
mean(f3[f3$`5UTR`==i&f3$mut=="sim"&f3$treatment=="NaAsO2",]$ratio3) #2.19468
mean(f3[f3$`5UTR`==i&f3$mut=="yak"&f3$treatment=="NaAsO2",]$ratio3) #1.72855

#ggarrange(p1,p2,p3,ncol=3,common.legend = T)

pdf("/results/fig5abc.pdf",width=9,height=3)
ggarrange(p1,p2,p3,ncol=3,common.legend = T)
dev.off()

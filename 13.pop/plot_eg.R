library(data.table)
library(ggplot2)
p=fread("/data/PBS_6pop_AXsep_specific_AF_DAFfilter_27.bed")

s1=data.table(t(p[p$V5=="FBgn0286567",23:28]))
s_gry=data.table(freq=s1$V1,pop=c("CN","FR","RAL","EF","SD","ZI"))
s_gry$pop=factor(s_gry$pop,levels=c("CN","FR","RAL","EF","SD","ZI"))
pdf("/results/fig6d.pdf",width=4,height=3)
ggplot(data=s_gry,aes(x=pop,y=freq,fill=pop))+
  geom_col()+
  ylim(0,1)+
  theme_classic()
dev.off()

s2=data.table(t(p[p$V5=="FBgn0033784",23:28]))
s_SCCRO3=data.table(freq=s2$V1,pop=c("CN","FR","RAL","EF","SD","ZI"))
s_SCCRO3$pop=factor(s_SCCRO3$pop,levels=c("CN","FR","RAL","EF","SD","ZI"))
pdf("/results/fig6e.pdf",width=4,height=3)
ggplot(data=s_SCCRO3,aes(x=pop,y=freq,fill=pop))+
  geom_col()+
  ylim(0,0.8)+
  theme_classic()
dev.off()

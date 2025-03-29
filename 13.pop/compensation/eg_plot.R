library(data.table)
library(ggplot2)
f1=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/eg/2R_10303372_10303346.txt",sep='\t')
f2=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/eg/X_14866339_14866345.txt",sep='\t')

colnames(f1)=c("freq","GT")
colnames(f2)=c("freq","GT")

f1$GT=factor(f1$GT,levels=c("AB","Ab","aB","ab"))
f2$GT=factor(f2$GT,levels=c("AB","Ab","aB","ab"))

ggplot(f1,aes(x=GT,y=freq,fill=GT))+
  geom_bar(stat="identity",width=0.8)+
  theme_classic()

ggplot(f2,aes(x=GT,y=freq,fill=GT))+
  geom_bar(stat="identity",width=0.8)+
  theme_classic()

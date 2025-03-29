## plot for alpha_asyptotic ####

library(ggplot2)
library(data.table)

mytheme <- theme_classic(base_size = 14) + theme(
  axis.text = element_text(color = 'black'),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  legend.text=element_text(size=9,face = "italic")
  #,axis.text.x = element_text(angle = 45,vjust = 0.9, hjust=1)
)

hh<-read.table("/data/MKtest_300bpDistance2CDS_result.txt",header = T,stringsAsFactors = F)
as.data.table(hh)->hh
hh$GL<-c("uORF gain","uORF gain","uORF gain","uORF loss","uORF loss","uORF loss")
hh$description<-c("Genes with uORF D. mel_gain without anc_loss",
                  "Genes with uORF D. mel_gain and anc_loss",
                  "All genes with uORF D. mel_gain",
                  "Genes with uORF D. mel_loss without anc_gain",
                  "Genes with uORF D. mel_loss and anc_gain",
                  "All genes with uORF D. mel_loss")

description_order=c("All genes with uORF D. mel_gain",
                    "Genes with uORF D. mel_gain and anc_loss",
                    "Genes with uORF D. mel_gain without anc_loss",
                    "All genes with uORF D. mel_loss",
                    "Genes with uORF D. mel_loss and anc_gain",
                    "Genes with uORF D. mel_loss without anc_gain")
hh$GL<-factor(hh$GL,levels = c("uORF gain","uORF loss"))
hh$description<-factor(hh$description,levels = description_order)

pdf("/results/figS22MK_300nt.pdf",width=12,height=4)
ggplot(hh, aes(x=GL, y=alpha_asymptotic,fill=description)) +
  geom_bar(aes(color = description,group=description,fill=description,),stat="identity",position = position_dodge(),width=0.8)+ 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),color="black",width = 0.2, position = position_dodge(0.8))+
  labs(x = '', y = 'alpha_asymptotic',color = NULL)+
  coord_cartesian(ylim = c(0, 1))+
  mytheme
dev.off()

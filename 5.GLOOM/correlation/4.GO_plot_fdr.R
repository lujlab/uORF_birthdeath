#metascape result plot
#q value

library(data.table)
library(ggplot2)
library(dplyr)
#all gainloss top500 
go=fread("/data/GO/metascape_input_all_gainloss_500.result",header=T,sep='\t')
go[, c("Group_num", "Group_summary") := tstrsplit(GroupID, "_", fixed = TRUE)]
go_summary=go[go$Group_summary=="Summary",]
go_summary$GOid=paste(go_summary$Term,go_summary$Description,sep=" ")
go_summary$FDR=-go_summary$`Log(q-value)`
go_summary$GOid=factor(go_summary$GOid,levels=rev(go_summary$GOid))
go_summary=go_summary[go_summary$`Log(q-value)`<log10(0.05)]
pdf(file="/results/figS7a_metascape_input_all_gainloss_500_go.pdf" , height= 4,width = 7, onefile = F)
ggplot(go_summary, aes(GOid, FDR)) +
  geom_bar(aes(fill = FDR), stat="identity", position="dodge", width=0.8, color = "#77776F")+
  theme_bw()+
  theme(panel.background=element_blank(),
        #        panel.grid=element_blank(),
        #        panel.border=element_blank(),
        axis.line=element_line(linewidth=0.3,colour="black"),
        axis.text.x = element_text(color="black", size=11),
        axis.text.y = element_text(color="black", size=11),
        panel.grid.major.y = element_blank(),
        legend.position = 'none')+
  labs(x=NULL, y="-log10(FDR)")+ 
  coord_flip() +
  scale_x_discrete(position = "top")+
  scale_fill_gradient2(low = "#FEE391",mid = "#FE9929", high = "#D95F0E",
                       limit = c(1,20),
                       midpoint = 10) 
dev.off()

#all netgain top500 
go=fread("/data/GO/metascape_input_all_netgain_500.result",header=T,sep='\t')
go[, c("Group_num", "Group_summary") := tstrsplit(GroupID, "_", fixed = TRUE)]
go_summary=go[go$Group_summary=="Summary",]
go_summary$GOid=paste(go_summary$Term,go_summary$Description,sep=" ")
go_summary$FDR=-go_summary$`Log(q-value)`
go_summary$GOid=factor(go_summary$GOid,levels=rev(go_summary$GOid))
go_summary=go_summary[go_summary$`Log(q-value)`<log10(0.05)]

pdf(file="/results/figS7b_metascape_input_all_netgain_500_go.pdf" , height= 4,width = 7, onefile = F)
ggplot(go_summary, aes(GOid, FDR)) +
  geom_bar(aes(fill = FDR), stat="identity", position="dodge", width=0.8, color = "#77776F")+
  theme_bw()+
  theme(panel.background=element_blank(),
        #        panel.grid=element_blank(),
        #        panel.border=element_blank(),
        axis.line=element_line(linewidth=0.3,colour="black"),
        axis.text.x = element_text(color="black", size=11),
        axis.text.y = element_text(color="black", size=11),
        panel.grid.major.y = element_blank(),
        legend.position = 'none')+
  labs(x=NULL, y="-log10(FDR)")+ 
  coord_flip() +
  scale_x_discrete(position = "top")+
  scale_fill_gradient2(low = "#FEE391",mid = "#FE9929", high = "#D95F0E",
                       limit = c(1,20),
                       midpoint = 10) 
dev.off()




#all canonical300nt gainloss top500 
go=fread("/data/GO/metascape_input_canonical300nt_gainloss_500.result",header=T,sep='\t')
go[, c("Group_num", "Group_summary") := tstrsplit(GroupID, "_", fixed = TRUE)]
go_summary=go[go$Group_summary=="Summary",]
go_summary$GOid=paste(go_summary$Term,go_summary$Description,sep=" ")
go_summary$FDR=-go_summary$`Log(q-value)`
go_summary$GOid=factor(go_summary$GOid,levels=rev(go_summary$GOid))
go_summary=go_summary[go_summary$`Log(q-value)`<log10(0.05)]

pdf(file="/results/figS8a_metascape_input_canonical300nt_gainloss_500.pdf" , height= 1,width = 5, onefile = F)
ggplot(go_summary, aes(GOid, FDR)) +
  geom_bar(aes(fill = FDR), stat="identity", position="dodge", width=0.8, color = "#77776F")+
  theme_bw()+
  theme(panel.background=element_blank(),
        #        panel.grid=element_blank(),
        #        panel.border=element_blank(),
        axis.line=element_line(linewidth=0.3,colour="black"),
        axis.text.x = element_text(color="black", size=11),
        axis.text.y = element_text(color="black", size=11),
        panel.grid.major.y = element_blank(),
        legend.position = 'none')+
  labs(x=NULL, y="-log10(FDR)")+ 
  coord_flip() +
  scale_x_discrete(position = "top")+
  scale_fill_gradient2(low = "#FEE391",mid = "#FE9929", high = "#D95F0E",
                       limit = c(1,20),
                       midpoint = 10) 
dev.off()

#all canonical300nt netgain top500 
go=fread("/data/GO/metascape_input_canonical300nt_netgain_500.result",header=T,sep='\t')
go[, c("Group_num", "Group_summary") := tstrsplit(GroupID, "_", fixed = TRUE)]
go_summary=go[go$Group_summary=="Summary",]
go_summary$GOid=paste(go_summary$Term,go_summary$Description,sep=" ")
go_summary$FDR=-go_summary$`Log(q-value)`
go_summary$GOid=factor(go_summary$GOid,levels=rev(go_summary$GOid))
go_summary=go_summary[go_summary$`Log(q-value)`<log10(0.05)]

pdf(file="/results/figS8b_metascape_input_canonical300nt_netgain_500_go.pdf" , height= 1,width = 6, onefile = F)
ggplot(go_summary, aes(GOid, FDR)) +
  geom_bar(aes(fill = FDR), stat="identity", position="dodge", width=0.8, color = "#77776F")+
  theme_bw()+
  theme(panel.background=element_blank(),
        #        panel.grid=element_blank(),
        #        panel.border=element_blank(),
        axis.line=element_line(linewidth=0.3,colour="black"),
        axis.text.x = element_text(color="black", size=11),
        axis.text.y = element_text(color="black", size=11),
        panel.grid.major.y = element_blank(),
        legend.position = 'none')+
  labs(x=NULL, y="-log10(FDR)")+ 
  coord_flip() +
  scale_x_discrete(position = "top")+
  scale_fill_gradient2(low = "#FEE391",mid = "#FE9929", high = "#D95F0E",
                       limit = c(1,20),
                       midpoint = 10) 
dev.off()

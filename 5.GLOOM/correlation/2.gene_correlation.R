library(ggplot2)
library(data.table)
#ATG number for gene
f1<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv")
f1_unique=f1[!duplicated(f1[,c("seqname","gene_pos")]),]
f1_unique_dm6=f1_unique[dm6==1,]
uATG_number_per_gene=data.frame(table(f1_unique_dm6$geneID))
colnames(uATG_number_per_gene)=c("FBid","uATG_number")

#dinuc stat
f2=fread("/data/correlation/dm6_merge_gene_stat_rmcds.txt",header=T)

gene_dinuc_uATG_number=merge(f2,uATG_number_per_gene,by="FBid",all=T)
gene_dinuc_uATG_number[is.na(gene_dinuc_uATG_number)]=0

##gain_loss_27sp


##### add gain loss number of each gene
jj3<-fread("/data/uORF_matrix_triCas2_ATG_GLOOME_result/gainLossMP.1.PerPosPerBranch.triCas2.unique.txt",header = T)
jj3[,.(num=.N),by=.(geneID,GL)]->GL_num_gene
dcast(GL_num_gene,geneID~GL,value.var ="num")->GL_num_gene2
colnames(GL_num_gene2)=c("FBid","gain","loss")
GL_num_gene2[is.na(GL_num_gene2)]=0

f3=merge(gene_dinuc_uATG_number,GL_num_gene2,by="FBid")

p=ggplot(f3, aes(x=log10(len), y=log10(loss+gain+1))) +
  geom_point(size=3,shape=16,fill="black",alpha=0.05)+
  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
  labs(x = 'log10(length of 5\' UTR)', y = 'log10(Total number of gain and loss events)', color = NULL)+
  #geom_abline(intercept = 0, slope = 1,color="red",linetype="dashed", size=0.8)+
  theme_classic() #fig4a_uATGgainloss_len_pergene
pdf("/results/figS5a_uATGgainloss_len_pergene.pdf",width=4,height=4)
print(p)
dev.off()
cor.test(log10(f3$len),log10(f3$loss+f3$gain+1),method = "spearman") #p-value < 2.2e-16   rho 0.8153781



p=ggplot(f3, aes(x=log10(len), y=log10(loss+1))) +
  geom_point(size=3,shape=16,fill="black",alpha=0.05)+
  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
  labs(x = 'log10(length of 5\' UTR)', y = 'log10(loss occurrence)', color = NULL)+
  #geom_abline(intercept = 0, slope = 1,color="red",linetype="dashed", size=0.8)+
  theme_classic() #fig_uATGloss_len_pergene
pdf("/results/figS5b_uATGloss_len_pergene.pdf",width=4,height=4)
print(p)
dev.off()
cor.test(log10(f3$len),log10(f3$loss+1),method = "spearman") #p-value < 2.2e-16   rho 0.6547831 

p=ggplot(f3, aes(x=log10(len), y=log10(gain+1))) +
  geom_point(size=3,shape=16,fill="black",alpha=0.05)+
  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
  labs(x = 'log10(length of 5\' UTR)', y = 'log10(gain occurrence)', color = NULL)+
  #geom_abline(intercept = 0, slope = 1,color="red",linetype="dashed", size=0.8)+
  theme_classic() #fig_uATGgain_len_pergene
pdf("/results/figS5c_uATGgain_len_pergene.pdf",width=4,height=4)
print(p)
dev.off()
cor.test(log10(f3$len),log10(f3$gain+1),method = "spearman") #p-value < 2.2e-16   rho 0.8124413 

p=ggplot(f3, aes(x=log10(gain+1), y=log10(loss+1))) +
  geom_point(size=3,shape=16,fill="black",alpha=0.05)+
  coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
  labs(x = 'log10(length of 5\' UTR)', y = 'log10(gain occurrence)', color = NULL)+
  geom_abline(intercept = 0, slope = 1,color="red",linetype="dashed", size=0.8)+
  theme_classic() #figS4d_uATGgain_loss_pergene
pdf("/results/figS5d_uATGgain_loss_pergene.pdf",width=4,height=4)
print(p)
dev.off()
cor.test(log10(f3$gain+1),log10(f3$loss+1),method = "spearman") #p-value < 2.2e-16   rho 0.7155545 


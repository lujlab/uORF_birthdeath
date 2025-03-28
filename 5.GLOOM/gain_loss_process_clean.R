#setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/10.GLOOME")

library(data.table)
library(ggplot2)
library(dplyr)
#library(pheatmap)

mytheme <- theme_classic(base_size = 16) + theme(
  axis.text = element_text(color = 'black',size=3),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  legend.text=element_text(size=9,face = "italic")
  #,axis.text.x = element_text(angle = 45,vjust = 0.9, hjust=1)
)

myorder<-c("dm6","PacBioSim","droSec1","N18","N17","droYak3","droEre2","N19","N16","droBia2","droSuz1","N20","N15","droAna3","droBip2","N21","N14","droEug2","N13","droEle2","N12","droKik2","N11","droTak2","N10","droRho2","N9","droFic2","N8",
           "droPse3","droPer1","N23","droMir2","N22","N7","droWil2","N6","droVir3","droMoj3","N26","droAlb1","N25","droGri2","N24","N5","musDom2","N4","anoGam1","N3","apiMel4","N2","triCas2","N1")
myorder2<-c("D.mel","D.sim","D.sec","N18","N17","D.yak","D.ere","N19","N16","D.bia","D.suz","N20","N15","D.ana","D.bip","N21","N14","D.eug","N13","D.ele","N12","D.kik","N11","D.tak","N10","D.rho","N9","D.fic","N8",
            "D.pse","D.per","N23","D.mir","N22","N7","D.wil","N6","D.vir","D.moj","N26","D.alb","N25","D.gri","N24","N5","M.dom","N4","A.gam","N3","A.mel","N2","T.cas","N1")
branch_name<-data.table(myorder,myorder2)
names(branch_name)<-c("branch","branch2")

## remove redundant uATG by unique their genome coordinates
hh<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv")

hh_rm=hh[hh$dm6+hh$PacBioSim+hh$droSec1+hh$droYak3+hh$droEre2+hh$droBia2+hh$droSuz1+hh$droAna3+hh$droBip2+hh$droEug2+hh$droEle2+hh$droKik2+hh$droTak2+hh$droRho2+hh$droFic2+hh$droPse3+hh$droPer1+hh$droMir2+hh$droWil2+hh$droVir3+hh$droMoj3+hh$droAlb1+hh$droGri2+hh$musDom2+hh$anoGam1+hh$apiMel4+hh$triCas2==0,]
hh[,c(1,2,29:34)]->hh
names(hh)[1]<-"ID"
hh2 <- tidyr::extract(hh, col = 'ID', into = c('tr', 'tr_pos_Gap','tr_pos_NoGap'),
                      regex = '(FBtr\\d+)_(\\d+)[(](\\d+)[)]',remove=F) %>% as.data.table()
hh2[,geno_pos:=paste(seqname,gene_pos,sep = "_")]
hh2[,c(1,5,2:4,9,12,11)]->hh2
hh2[,POS:=1:.N]
hh2[order(geno_pos,POS)]->hh2 #total 435425(including drosim1),434416 uATG
hh2[!duplicated(geno_pos)]->hh3 #total 224879(including drosim1),224358 unique uATG

## uATG number identified in each branch
nn<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv")
nn=nn[nn$dm6+nn$PacBioSim+nn$droSec1+nn$droYak3+nn$droEre2+nn$droBia2+nn$droSuz1+nn$droAna3+nn$droBip2+nn$droEug2+nn$droEle2+nn$droKik2+nn$droTak2+nn$droRho2+nn$droFic2+nn$droPse3+nn$droPer1+nn$droMir2+nn$droWil2+nn$droVir3+nn$droMoj3+nn$droAlb1+nn$droGri2+nn$musDom2+nn$anoGam1+nn$apiMel4+nn$triCas2>0,]

names(nn)[1]<-"ID"
nn[,geno_pos:=paste(seqname,gene_pos,sep = "_")]
nn[,POS:=1:.N]
nn[order(geno_pos,POS)]->nn
nn[!duplicated(geno_pos)]->nn2
#nn[!duplicated(dm6,geno_pos)]->nn2
#nn2=nn[!duplicated(nn, by = c("dm6", "geno_pos"))]


nn2[,2:28]->nn3
apply(nn3, 2, sum)->tmp
sum(tmp) # total 668375 uAUGs
names(tmp)<-myorder2[grep("N",myorder2,invert = T)]
as.data.frame(tmp)->tmp1
tmp1$species<-rownames(tmp1)
as.data.table(tmp1)->tmp1
tmp1$species<-factor(tmp1$species,levels = myorder2)
names(tmp1)[1]<-"num"
p=ggplot(tmp1, aes(x=species, y=num)) +
  geom_bar(stat="identity",fill="gray",position = position_dodge(),width=0.8)+
  coord_flip()+
  labs(x = '', y = 'uATG number',color = NULL)+
  mytheme+
  theme(axis.text.y = element_text(face = "italic"))
pdf("/results/figS1b.pdf",width=4,height=4)
print(p)
dev.off()

## add uATG gain/loss events with detailed information
jj<-fread("/code/5.GLOOM/uORF_matrix_triCas2_ATG_GLOOME_result/gainLossMP.1.PerPosPerBranch.txt")
names(jj)[1]<-"GL"
merge(jj,hh2,by="POS")->jj2
fwrite(jj2,"/results/gainLossMP.1.PerPosPerBranch.triCas2.all.txt",sep="\t")
merge(jj,hh3,by="POS")->jj3
fwrite(jj3,"/results/gainLossMP.1.PerPosPerBranch.triCas2.unique.txt",sep="\t")
nrow(jj3) # 335686 GL events
length(jj3[jj3$GL=="gain",]$POS) #302849
length(jj3[jj3$GL=="loss",]$POS) #32837

length(unique(jj3$POS)) # 224319 uATG with gain or loss events happened
x=hh3[!(hh3$ID %in% hh_rm$V1),]
## gain/loss event number for each uATG
tmp<-jj3[,.(num=.N),by=.(POS)]
table(tmp$num)
#1      2      3      4      5      6      7      8      9     10 
#165846  30262  13048   8482   4514   1613    454     83     16      1 


jj3<-fread("/data/gainLossMP.1.PerPosPerBranch.triCas2.unique.txt",header = T)
merge(jj3,branch_name,by="branch")->jj3
jj3[,.(num=.N),by=.(branch2,GL)]->GL_num_branch #number for figS1c

GL_num_branch$branch2<-factor(GL_num_branch$branch2,levels = myorder2)
GL_num_branch$GL<-factor(GL_num_branch$GL,levels = c("gain","loss")) 

dcast(GL_num_branch,branch2~GL,value.var = 'num')->GL_num_branch2
GL_num_branch2[is.na(GL_num_branch2)]<-0
GL_num_branch2$total=GL_num_branch2$gain+GL_num_branch2$loss
GL_num_branch2$gain_raio=GL_num_branch2$gain/GL_num_branch2$total
summary(GL_num_branch[GL=="gain"]$num) # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 7    2120    4458    5824    7615   16600 
summary(GL_num_branch[GL=="loss"]$num) # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 1.0    76.5   323.0   656.7   851.5  5055.0 


## each uATG gain/loss in each brach
dcast(jj3,POS~branch2,value.var = 'GL')->jj4
jj4[is.na(jj4)]<-0
jj4[jj4=="gain"]<-1
jj4[jj4=="loss"]<- -1
as.data.frame(jj4)->jj4
rownames(jj4)<-as.character(jj4$POS)
jj4[,-1]->jj4
jj4[,myorder2[-53]]->jj4 # reorder the column
# most GL as example: Abd-B
# Abd-B:uATG presence or absence in each branch
#aa<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv")
#aa[,c(1,2,29:34)]->aa
#names(aa)[1]<-"ID"
#aa2 <- tidyr::extract(aa, col = 'ID', into = c('tr', 'tr_pos_Gap','tr_pos_NoGap'),
#                      regex = '(FBtr\\d+)_(\\d+)[(](\\d+)[)]',remove=F) %>% as.data.table
#aa2[,geno_pos:=paste(seqname,gene_pos,sep = "_")]

#aa2[,POS:=1:.N]
#aa2[order(geno_pos,POS)]->aa2 #total 435425 uATG
#aa2[!duplicated(geno_pos)]->aa3 #total 224879
#names(aa3)[5:31]<-myorder2[grep("N",myorder2,invert = T)]
#aa3[geneID=="FBgn0000015",31:5]->aa3_abd
#pabd=pheatmap(aa3_abd,cluster_rows=0,cluster_cols=0,treeheight_row = 0,treeheight_col = 0,show_rownames = 0,show_colnames = 1) # ,col=colorRampPalette(c("blue","gray","red"))(3)
#pdf("/results/figS4.pdf",width=4,height=4)
#print(pabd)
#dev.off()

## uATG gain/loss of each gene in all branch
jj3=as.data.table(jj3)
jj3[,.(num=.N),by=.(geneID,GL)]->GL_num_gene
dcast(GL_num_gene,geneID~GL,value.var ="num")->GL_num_gene2
GL_num_gene2=as.data.table(GL_num_gene2)
#GL_num_gene2[!is.na(geneID)]->GL_num_gene2
GL_num_gene2[is.na(GL_num_gene2)]<-0


p2=ggplot(GL_num_gene2, aes(x=log10(gain+1), y=log10(loss+1))) +
  geom_point(size=3,shape=16,fill="black",alpha=0.05)+
  coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
  labs(x = 'log10(gain occurrence)', y = 'log10(loss occurrence)', color = NULL)+
  geom_abline(intercept = 0, slope = 1,color="red",linetype="dashed", size=0.8)+
  mytheme+theme_classic(base_size = 26,base_line_size = 1 )
pdf("/results/figS2b.pdf",width=4,height=4)
print(p2)
dev.off()

#cor.test(log10(GL_num_gene2$gain+1),log10(GL_num_gene2$loss+1),method = "pearson")
cor.test(log10(GL_num_gene2$gain+1),log10(GL_num_gene2$loss+1),method = "spearman")
## ratio of uATG gain-to-loss of each gene in all branch
GL_num_gene2=as.data.table(GL_num_gene2)
GL_num_gene2[,c("gain_ratio","total","net_gain","net_loss"):= .(gain/(gain+loss),loss+gain,gain - loss,loss - gain)]
gene_num<-length(unique(GL_num_gene2$geneID))
GL_num_gene2->tmp
tmp[gain_ratio>1]$gain_ratio <- 1.01
hist(tmp$gain_ratio,breaks = 50,col = "gray",las=1,xlab="Gain ratio",ylab="Count",main="The distribution of percent of gain events in each gene")
## total events of  uATG gain and loss of each gene in all branch correlated with 5UTR length
utr_len<-fread("/data/dmel-all-r6.04-gene-merged-5utr.length",header = T)
GL_num_gene2->tmp
merge(tmp,utr_len[,c(1,5)],by.x="geneID",by.y="gene")->tmp2

p3=ggplot(tmp2, aes(x=log10(merged+1), y=log10(total+1))) +
  geom_point(size=3,shape=16,fill="black",alpha=0.03)+
  coord_cartesian(xlim = c(0, 4),ylim = c(0, 3))+
  labs(x = 'log10(5UTR length)', y = 'log10(Total number of gain and loss)', color = NULL)+
  geom_smooth(method='lm',color="red",linetype="dashed", size=0.8, se=FALSE)+
  mytheme+theme_classic(base_size = 26,base_line_size = 1 )
pdf("/results/figS2a.pdf",width=4,height=4)
print(p3)
dev.off()

cor.test(log10(tmp2$total+1),log10(tmp2$merged+1),method = "spearman")



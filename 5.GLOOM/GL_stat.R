
library(data.table)
library(ggplot2)
library(dplyr)
library(pheatmap)

myorder<-c("dm6","PacBioSim","droSec1","N18","N17","droYak3","droEre2","N19","N16","droBia2","droSuz1","N20","N15","droAna3","droBip2","N21","N14","droEug2","N13","droEle2","N12","droKik2","N11","droTak2","N10","droRho2","N9","droFic2","N8",
           "droPse3","droPer1","N23","droMir2","N22","N7","droWil2","N6","droVir3","droMoj3","N26","droAlb1","N25","droGri2","N24","N5","musDom2","N4","anoGam1","N3","apiMel4","N2","triCas2","N1")
myorder2<-c("D.mel","D.sim","D.sec","N18","N17","D.yak","D.ere","N19","N16","D.bia","D.suz","N20","N15","D.ana","D.bip","N21","N14","D.eug","N13","D.ele","N12","D.kik","N11","D.tak","N10","D.rho","N9","D.fic","N8",
            "D.pse","D.per","N23","D.mir","N22","N7","D.wil","N6","D.vir","D.moj","N26","D.alb","N25","D.gri","N24","N5","M.dom","N4","A.gam","N3","A.mel","N2","T.cas","N1")
branch_name<-data.table(myorder,myorder2)
names(branch_name)<-c("branch","branch2")

## remove redundant uATG by unique their genome coordinates
hh<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv")
#hh$total=hh$dm6+hh$PacBioSim+hh$droSec1+hh$droYak3+hh$droEre2+hh$droBia2+hh$droSuz1+hh$droAna3+hh$droBip2+hh$droEug2+hh$droEle2+hh$droKik2+hh$droTak2+hh$droRho2+hh$droFic2+hh$droPse3+hh$droPer1+hh$droMir2+hh$droWil2+hh$droVir3+hh$droMoj3+hh$droAlb1+hh$droGri2+hh$musDom2+hh$anoGam1+hh$apiMel4+hh$triCas2
hh_rm=hh[hh$dm6+hh$PacBioSim+hh$droSec1+hh$droYak3+hh$droEre2+hh$droBia2+hh$droSuz1+hh$droAna3+hh$droBip2+hh$droEug2+hh$droEle2+hh$droKik2+hh$droTak2+hh$droRho2+hh$droFic2+hh$droPse3+hh$droPer1+hh$droMir2+hh$droWil2+hh$droVir3+hh$droMoj3+hh$droAlb1+hh$droGri2+hh$musDom2+hh$anoGam1+hh$apiMel4+hh$triCas2==0,]
hh[,c(1,2,29:34)]->hh
names(hh)[1]<-"ID"
hh2 <- tidyr::extract(hh, col = 'ID', into = c('tr', 'tr_pos_Gap','tr_pos_NoGap'),
                      regex = '(FBtr\\d+)_(\\d+)[(](\\d+)[)]',remove=F) %>% as.data.table()
hh2[,geno_pos:=paste(seqname,gene_pos,sep = "_")]
hh2[,c(1,5,2:4,9,12,11)]->hh2
hh2[,POS:=1:.N]
hh2[order(geno_pos,POS)]->hh2 #total 435425
hh2[!duplicated(geno_pos)]->hh3 #total 224879


#GainLoss
jj<-fread("/data/uORF_matrix_triCas2_ATG_GLOOME_result/gainLossMP.1.PerPosPerBranch.txt")
names(jj)[1]<-"GL"
#merge(jj,hh2,by="POS")->jj2
#fwrite(jj2,"gainLossMP.1.PerPosPerBranch.triCas2.all.txt",sep="\t")
merge(jj,hh3,by="POS")->jj3
#fwrite(jj3,"gainLossMP.1.PerPosPerBranch.triCas2.unique.txt",sep="\t")
nrow(jj3) # 335686 GL events
length(jj3[jj3$GL=="gain",]$POS) #302849
length(jj3[jj3$GL=="loss",]$POS) #32837

length(unique(jj3$POS)) # 224319 uATG with gain or loss events happened

jj3<-fread("/data/uORF_matrix_triCas2_ATG_GLOOME_result/gainLossMP.1.PerPosPerBranch.triCas2.unique.txt",header = T)
merge(jj3,branch_name,by="branch")->jj3
jj3[,.(num=.N),by=.(branch2,GL)]->GL_num_branch

GL_num_branch$branch2<-factor(GL_num_branch$branch2,levels = myorder2)
GL_num_branch$GL<-factor(GL_num_branch$GL,levels = c("gain","loss"))
#ggplot(GL_num_branch, aes(x = branch2, y = num,fill=GL))+
#  geom_bar(position="dodge", stat="identity")+
#  labs(x = 'Branch', y = 'Number of Gain/Loss events', color = NULL)+
#  mytheme

dcast(GL_num_branch,branch2~GL,value.var = 'num')->GL_num_branch2
GL_num_branch2[is.na(GL_num_branch2)]<-0
#GL_num_branch2[,c("total","gain_raio"):=.(gain+loss,gain/(gain+loss))]
GL_num_branch2$total=GL_num_branch2$gain+GL_num_branch2$loss
GL_num_branch2$gain_raio=GL_num_branch2$gain/GL_num_branch2$total


fwrite(GL_num_branch2,"/results/figS2c_GL_events.txt",sep='\t') #data for figS2c



#####Abd-B

####Abd-B
jj3=data.table(jj3)
jj3_abdb=jj3[geneID=="FBgn0000015",]
jj3_abdb[,.(num=.N),by=.(branch2,GL)]->GL_num_branch_abdb
GL_num_branch_abdb$branch2<-factor(GL_num_branch_abdb$branch2,levels = myorder2)
GL_num_branch_abdb$GL<-factor(GL_num_branch_abdb$GL,levels = c("gain","loss"))
dcast(GL_num_branch_abdb,branch2~GL,value.var = 'num')->GL_num_branch_abdb2
GL_num_branch_abdb2[is.na(GL_num_branch_abdb2)]<-0
GL_num_branch_abdb2$total=GL_num_branch_abdb2$gain+GL_num_branch_abdb2$loss
GL_num_branch_abdb2$gain_raio=GL_num_branch_abdb2$gain/GL_num_branch_abdb2$total
fwrite(GL_num_branch_abdb2,"/results/figS8b_GL_abdb_events.txt",sep='\t') #data for figS8b

###heatmap
hh<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv")
hh=hh[hh$dm6+hh$PacBioSim+hh$droSec1+hh$droYak3+hh$droEre2+hh$droBia2+hh$droSuz1+hh$droAna3+hh$droBip2+hh$droEug2+hh$droEle2+hh$droKik2+hh$droTak2+hh$droRho2+hh$droFic2+hh$droPse3+hh$droPer1+hh$droMir2+hh$droWil2+hh$droVir3+hh$droMoj3+hh$droAlb1+hh$droGri2+hh$musDom2+hh$anoGam1+hh$apiMel4+hh$triCas2>0,]

names(hh)[1]<-"ID"
hh2 <- tidyr::extract(hh, col = 'ID', into = c('tr', 'tr_pos_Gap','tr_pos_NoGap'),
                      regex = '(FBtr\\d+)_(\\d+)[(](\\d+)[)]',remove=F) %>% as.data.table()
hh2=as.data.table(hh2)
hh2[,POS:=1:.N]

hh2[,geno_pos:=paste(seqname,gene_pos,sep = "_")]
hh2[!duplicated(geno_pos)]->hh3
hh3[order(geno_pos,POS)]->hh3
hh3[hh3$dm6==1,]->hh4
gene_uORF_num=as.data.table(table(hh4$geneID))

abd_uORF=hh3[hh3$geneID=="FBgn0000015",c(5:31)] 
sum(abd_uORF$dm6) #112
#reverse
abd_uORF[,POS:=1:.N]
abd_uORF[order(-POS)]->abd_uORF
abd_uORF=abd_uORF[,c(1:27)]
tmp<-sapply(abd_uORF, as.numeric)


pdf("/results/figS8A.pdf",width=3.8,height=4.6)
pheatmap(tmp,cluster_rows=0,cluster_cols=0,treeheight_row = 0,treeheight_col = 0,show_rownames = 0,show_colnames = 1,legend = F) 
dev.off()



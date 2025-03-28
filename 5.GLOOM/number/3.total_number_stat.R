library(data.table)
library(ggplot2)
library(dplyr)
myorder<-c("dm6","PacBioSim","droSim1","droSec1","N18","N17","droYak3","droEre2","N19","N16","droBia2","droSuz1","N20","N15","droAna3","droBip2","N21","N14","droEug2","N13","droEle2","N12","droKik2","N11","droTak2","N10","droRho2","N9","droFic2","N8",
           "droPse3","droPer1","N23","droMir2","N22","N7","droWil2","N6","droVir3","droMoj3","N26","droAlb1","N25","droGri2","N24","N5","musDom2","N4","anoGam1","N3","apiMel4","N2","triCas2","N1")
myorder2<-c("D.mel","D.sim_PacBio","D.sim_droSim1","D.sec","N18","N17","D.yak","D.ere","N19","N16","D.bia","D.suz","N20","N15","D.ana","D.bip","N21","N14","D.eug","N13","D.ele","N12","D.kik","N11","D.tak","N10","D.rho","N9","D.fic","N8",
            "D.pse","D.per","N23","D.mir","N22","N7","D.wil","N6","D.vir","D.moj","N26","D.alb","N25","D.gri","N24","N5","M.dom","N4","A.gam","N3","A.mel","N2","T.cas","N1")
branch_name<-data.table(myorder,myorder2)
names(branch_name)<-c("branch","branch2")


hh<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv")

hh[,c(1,2,29:34)]->hh
names(hh)[1]<-"ID"
hh2 <- tidyr::extract(hh, col = 'ID', into = c('tr', 'tr_pos_Gap','tr_pos_NoGap'),
                      regex = '(FBtr\\d+)_(\\d+)[(](\\d+)[)]',remove=F) %>% as.data.table()
hh2[,geno_pos:=paste(seqname,gene_pos,sep = "_")]
hh2[,c(1,5,2:4,9,12,11)]->hh2
hh2[,POS:=1:.N]
hh2[order(geno_pos,POS)]->hh2 #total 435425
hh2[!duplicated(geno_pos)]->hh3 #total 224879

## uATG number identified in each branch
nn<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv")

names(nn)[1]<-"ID"
nn[,geno_pos:=paste(seqname,gene_pos,sep = "_")]
nn[,POS:=1:.N]
nn[order(geno_pos,POS)]->nn
nn[!duplicated(geno_pos)]->nn2
#nn[!duplicated(dm6,geno_pos)]->nn2
#nn2=nn[!duplicated(nn, by = c("dm6", "geno_pos"))]

nn2[,2:29]->nn3
apply(nn3, 2, sum)->tmp
#sum(tmp) 
names(tmp)<-myorder2[grep("N",myorder2,invert = T)]
as.data.frame(tmp)->tmp1
tmp1$species<-rownames(tmp1)
as.data.table(tmp1)->tmp1
tmp1$species<-factor(tmp1$species,levels = myorder2)
names(tmp1)[1]<-"num"
fwrite(tmp1,"/results/sp27_uATG_number.txt",sep='\t')

pdf("/results/figS2b.pdf",width=3.2,height=4)

ggplot(tmp1, aes(x=species, y=num)) +
  geom_bar(stat="identity",fill="gray",position = position_dodge(),width=0.8)+
  coord_flip()+
  labs(x = '', y = 'uATG number',color = NULL)+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic")) 

dev.off()
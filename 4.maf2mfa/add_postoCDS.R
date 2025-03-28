library(data.table)
#myorder<-c("dm6","PacBioSim","droSec1","N18","N17","droYak3","droEre2","N19","N16","droBia2","droSuz1","N20","N15","droAna3","droBip2","N21","N14","droEug2","N13","droEle2","N12","droKik2","N11","droTak2","N10","droRho2","N9","droFic2","N8",
#           "droPse3","droPer1","N23","droMir2","N22","N7","droWil2","N6","droVir3","droMoj3","N26","droAlb1","N25","droGri2","N24","N5","musDom2","N4","anoGam1","N3","apiMel4","N2","triCas2","N1")
#myorder2<-c("D.mel","D.sim","D.sec","N18","N17","D.yak","D.ere","N19","N16","D.bia","D.suz","N20","N15","D.ana","D.bip","N21","N14","D.eug","N13","D.ele","N12","D.kik","N11","D.tak","N10","D.rho","N9","D.fic","N8",
#            "D.pse","D.per","N23","D.mir","N22","N7","D.wil","N6","D.vir","D.moj","N26","D.alb","N25","D.gri","N24","N5","M.dom","N4","A.gam","N3","A.mel","N2","T.cas","N1")
hh<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",header=T)
colnames(hh)[1]="ID"
len<-fread("/code/4.maf2mfa/dm6_5UTR_nogap_len.txt",header=F)
colnames(len)=c("transcriptID","utr_len")

h2 <- merge(hh, len, by = "transcriptID")

h2$pos_to_cds=h2$utr_len-h2$position_withoutgap

fwrite(h2,"/results/uORF_matrix_dm6_PacBioSim_27species_update_noCDS_postoCDS.csv",sep=',')

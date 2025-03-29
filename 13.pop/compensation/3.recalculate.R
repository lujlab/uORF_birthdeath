library(data.table)
CN=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/CN_gain_loss_filter_new.txt",sep='\t',na="NA",quote=F)
FR=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/FR_gain_loss_filter_new.txt",sep='\t',na="NA",quote=F)
RAL=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/RAL_gain_loss_filter_new.txt",sep='\t',na="NA",quote=F)
EF=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/EF_gain_loss_filter_new.txt",sep='\t',na="NA",quote=F)
SD=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/SD_gain_loss_filter_new.txt",sep='\t',na="NA",quote=F)
ZI=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/ZI_gain_loss_filter_new.txt",sep='\t',na="NA",quote=F)

#MAF>0.05
#EF=EF[EF$gain_freq>0.05&EF$gain_freq<0.95&EF$loss_freq>0.05&EF$gain_freq<0.95]
#FR=FR[FR$gain_freq>0.05&FR$gain_freq<0.95&FR$loss_freq>0.05&FR$gain_freq<0.95]
#RAL=RAL[RAL$gain_freq>0.05&RAL$gain_freq<0.95&RAL$loss_freq>0.05&RAL$gain_freq<0.95]
#EF=EF[EF$gain_freq>0.05&EF$gain_freq<0.95&EF$loss_freq>0.05&EF$gain_freq<0.95]
#SD=SD[SD$gain_freq>0.05&SD$gain_freq<0.95&SD$loss_freq>0.05&SD$gain_freq<0.95]
#ZI=ZI[ZI$gain_freq>0.05&ZI$gain_freq<0.95&ZI$loss_freq>0.05&ZI$gain_freq<0.95]

#CN
CN$A_freq=(CN$Ab_new+CN$AB_new)/(CN$AB_new+CN$Ab_new+CN$aB_new+CN$ab_new)
CN$a_freq=(CN$ab_new+CN$aB_new)/(CN$AB_new+CN$Ab_new+CN$aB_new+CN$ab_new)
CN$B_freq=(CN$aB_new+CN$AB_new)/(CN$AB_new+CN$Ab_new+CN$aB_new+CN$ab_new)
CN$b_freq=(CN$ab_new+CN$Ab_new)/(CN$AB_new+CN$Ab_new+CN$aB_new+CN$ab_new)

CN$AB_Expected=CN$A_freq*CN$B_freq*(CN$AB_new+CN$Ab_new+CN$aB_new+CN$ab_new)
CN$Ab_Expected=CN$A_freq*CN$b_freq*(CN$AB_new+CN$Ab_new+CN$aB_new+CN$ab_new)
CN$aB_Expected=CN$a_freq*CN$B_freq*(CN$AB_new+CN$Ab_new+CN$aB_new+CN$ab_new)
CN$ab_Expected=CN$a_freq*CN$b_freq*(CN$AB_new+CN$Ab_new+CN$aB_new+CN$ab_new)
CN$p=NA
###re-calculate correlation
for (i in c(1:nrow(CN))){
  CN$p[i]=fisher.test(matrix(c(CN$AB_new[i],CN$Ab_new[i],CN$aB_new[i],CN$ab_new[i]),nrow=2))$p
}
fwrite(CN,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/CN_gain_loss_p.txt")


#FR
FR$A_freq=(FR$Ab_new+FR$AB_new)/(FR$AB_new+FR$Ab_new+FR$aB_new+FR$ab_new)
FR$a_freq=(FR$ab_new+FR$aB_new)/(FR$AB_new+FR$Ab_new+FR$aB_new+FR$ab_new)
FR$B_freq=(FR$aB_new+FR$AB_new)/(FR$AB_new+FR$Ab_new+FR$aB_new+FR$ab_new)
FR$b_freq=(FR$ab_new+FR$Ab_new)/(FR$AB_new+FR$Ab_new+FR$aB_new+FR$ab_new)

FR$AB_Expected=FR$A_freq*FR$B_freq*(FR$AB_new+FR$Ab_new+FR$aB_new+FR$ab_new)
FR$Ab_Expected=FR$A_freq*FR$b_freq*(FR$AB_new+FR$Ab_new+FR$aB_new+FR$ab_new)
FR$aB_Expected=FR$a_freq*FR$B_freq*(FR$AB_new+FR$Ab_new+FR$aB_new+FR$ab_new)
FR$ab_Expected=FR$a_freq*FR$b_freq*(FR$AB_new+FR$Ab_new+FR$aB_new+FR$ab_new)
FR$p=NA
###re-calculate correlation
for (i in c(1:nrow(FR))){
  FR$p[i]=fisher.test(matrix(c(FR$AB_new[i],FR$Ab_new[i],FR$aB_new[i],FR$ab_new[i]),nrow=2))$p
}
fwrite(FR,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/FR_gain_loss_p.txt")


#RAL
RAL$A_freq=(RAL$Ab_new+RAL$AB_new)/(RAL$AB_new+RAL$Ab_new+RAL$aB_new+RAL$ab_new)
RAL$a_freq=(RAL$ab_new+RAL$aB_new)/(RAL$AB_new+RAL$Ab_new+RAL$aB_new+RAL$ab_new)
RAL$B_freq=(RAL$aB_new+RAL$AB_new)/(RAL$AB_new+RAL$Ab_new+RAL$aB_new+RAL$ab_new)
RAL$b_freq=(RAL$ab_new+RAL$Ab_new)/(RAL$AB_new+RAL$Ab_new+RAL$aB_new+RAL$ab_new)

RAL$AB_Expected=RAL$A_freq*RAL$B_freq*(RAL$AB_new+RAL$Ab_new+RAL$aB_new+RAL$ab_new)
RAL$Ab_Expected=RAL$A_freq*RAL$b_freq*(RAL$AB_new+RAL$Ab_new+RAL$aB_new+RAL$ab_new)
RAL$aB_Expected=RAL$a_freq*RAL$B_freq*(RAL$AB_new+RAL$Ab_new+RAL$aB_new+RAL$ab_new)
RAL$ab_Expected=RAL$a_freq*RAL$b_freq*(RAL$AB_new+RAL$Ab_new+RAL$aB_new+RAL$ab_new)
RAL$p=NA
###re-calculate correlation
for (i in c(1:nrow(RAL))){
  RAL$p[i]=fisher.test(matrix(c(RAL$AB_new[i],RAL$Ab_new[i],RAL$aB_new[i],RAL$ab_new[i]),nrow=2))$p
}
fwrite(RAL,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/RAL_gain_loss_p.txt")


#EF
EF$A_freq=(EF$Ab_new+EF$AB_new)/(EF$AB_new+EF$Ab_new+EF$aB_new+EF$ab_new)
EF$a_freq=(EF$ab_new+EF$aB_new)/(EF$AB_new+EF$Ab_new+EF$aB_new+EF$ab_new)
EF$B_freq=(EF$aB_new+EF$AB_new)/(EF$AB_new+EF$Ab_new+EF$aB_new+EF$ab_new)
EF$b_freq=(EF$ab_new+EF$Ab_new)/(EF$AB_new+EF$Ab_new+EF$aB_new+EF$ab_new)

EF$AB_Expected=EF$A_freq*EF$B_freq*(EF$AB_new+EF$Ab_new+EF$aB_new+EF$ab_new)
EF$Ab_Expected=EF$A_freq*EF$b_freq*(EF$AB_new+EF$Ab_new+EF$aB_new+EF$ab_new)
EF$aB_Expected=EF$a_freq*EF$B_freq*(EF$AB_new+EF$Ab_new+EF$aB_new+EF$ab_new)
EF$ab_Expected=EF$a_freq*EF$b_freq*(EF$AB_new+EF$Ab_new+EF$aB_new+EF$ab_new)
EF$p=NA
###re-calculate correlation
for (i in c(1:nrow(EF))){
  EF$p[i]=fisher.test(matrix(c(EF$AB_new[i],EF$Ab_new[i],EF$aB_new[i],EF$ab_new[i]),nrow=2))$p
}
fwrite(EF,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/EF_gain_loss_p.txt")


#SD
SD$A_freq=(SD$Ab_new+SD$AB_new)/(SD$AB_new+SD$Ab_new+SD$aB_new+SD$ab_new)
SD$a_freq=(SD$ab_new+SD$aB_new)/(SD$AB_new+SD$Ab_new+SD$aB_new+SD$ab_new)
SD$B_freq=(SD$aB_new+SD$AB_new)/(SD$AB_new+SD$Ab_new+SD$aB_new+SD$ab_new)
SD$b_freq=(SD$ab_new+SD$Ab_new)/(SD$AB_new+SD$Ab_new+SD$aB_new+SD$ab_new)

SD$AB_Expected=SD$A_freq*SD$B_freq*(SD$AB_new+SD$Ab_new+SD$aB_new+SD$ab_new)
SD$Ab_Expected=SD$A_freq*SD$b_freq*(SD$AB_new+SD$Ab_new+SD$aB_new+SD$ab_new)
SD$aB_Expected=SD$a_freq*SD$B_freq*(SD$AB_new+SD$Ab_new+SD$aB_new+SD$ab_new)
SD$ab_Expected=SD$a_freq*SD$b_freq*(SD$AB_new+SD$Ab_new+SD$aB_new+SD$ab_new)
SD$p=NA
###re-calculate correlation
for (i in c(1:nrow(SD))){
  SD$p[i]=fisher.test(matrix(c(SD$AB_new[i],SD$Ab_new[i],SD$aB_new[i],SD$ab_new[i]),nrow=2))$p
}
fwrite(SD,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/SD_gain_loss_p.txt")


#ZI
ZI$A_freq=(ZI$Ab_new+ZI$AB_new)/(ZI$AB_new+ZI$Ab_new+ZI$aB_new+ZI$ab_new)
ZI$a_freq=(ZI$ab_new+ZI$aB_new)/(ZI$AB_new+ZI$Ab_new+ZI$aB_new+ZI$ab_new)
ZI$B_freq=(ZI$aB_new+ZI$AB_new)/(ZI$AB_new+ZI$Ab_new+ZI$aB_new+ZI$ab_new)
ZI$b_freq=(ZI$ab_new+ZI$Ab_new)/(ZI$AB_new+ZI$Ab_new+ZI$aB_new+ZI$ab_new)

ZI$AB_Expected=ZI$A_freq*ZI$B_freq*(ZI$AB_new+ZI$Ab_new+ZI$aB_new+ZI$ab_new)
ZI$Ab_Expected=ZI$A_freq*ZI$b_freq*(ZI$AB_new+ZI$Ab_new+ZI$aB_new+ZI$ab_new)
ZI$aB_Expected=ZI$a_freq*ZI$B_freq*(ZI$AB_new+ZI$Ab_new+ZI$aB_new+ZI$ab_new)
ZI$ab_Expected=ZI$a_freq*ZI$b_freq*(ZI$AB_new+ZI$Ab_new+ZI$aB_new+ZI$ab_new)
ZI$p=NA
###re-calculate correlation
for (i in c(1:nrow(ZI))){
  ZI$p[i]=fisher.test(matrix(c(ZI$AB_new[i],ZI$Ab_new[i],ZI$aB_new[i],ZI$ab_new[i]),nrow=2))$p
}
fwrite(ZI,"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/6.compensation/6pop/ZI_gain_loss_p.txt")


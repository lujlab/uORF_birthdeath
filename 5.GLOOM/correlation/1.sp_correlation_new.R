library(Rmisc)
library(data.table)
library(ggplot2)
myorder<-c("dm6","PacBioSim","droSec1","N18","N17","droYak3","droEre2","N19","N16","droBia2","droSuz1","N20","N15","droAna3","droBip2","N21","N14","droEug2","N13","droEle2","N12","droKik2","N11","droTak2","N10","droRho2","N9","droFic2","N8",
           "droPse3","droPer1","N23","droMir2","N22","N7","droWil2","N6","droVir3","droMoj3","N26","droAlb1","N25","droGri2","N24","N5","musDom2","N4","anoGam1","N3","apiMel4","N2","triCas2","N1")
myorder2<-c("D.mel","D.sim","D.sec","N18","N17","D.yak","D.ere","N19","N16","D.bia","D.suz","N20","N15","D.ana","D.bip","N21","N14","D.eug","N13","D.ele","N12","D.kik","N11","D.tak","N10","D.rho","N9","D.fic","N8",
            "D.pse","D.per","N23","D.mir","N22","N7","D.wil","N6","D.vir","D.moj","N26","D.alb","N25","D.gri","N24","N5","M.dom","N4","A.gam","N3","A.mel","N2","T.cas","N1")
myorder3<-c("D.mel","D.sim","D.sec","D.yak","D.ere","D.bia","D.suz","D.ana","D.bip","D.eug","D.ele","D.kik","D.tak","D.rho","D.fic",
            "D.pse","D.per","D.mir","D.wil","D.vir","D.moj","D.alb","D.gri","M.dom","A.gam","A.mel","T.cas")
myorder4<-c("D.mel","D.sim","D.sec","D.yak","D.ere","D.bia","D.suz","D.ana","D.bip","D.eug","D.ele","D.kik","D.tak","D.rho","D.fic",
            "D.pse","D.per","D.mir","D.wil","D.vir","D.moj","D.alb","D.gri")


branch_name<-data.table(myorder,myorder2)
names(branch_name)<-c("branch","branch2")

f1=fread("/data/correlation/sp27_uATG_number.txt",header=T)
f1=f1[f1$species!="D.sim_droSim1",]
f1[f1$species=="D.sim_PacBio",]$species="D.sim"
colnames(f1)[2]="branch2"
f1=merge(f1,branch_name,by="branch2")

f2=fread("/data/correlation/stat_merged_27sp.txt",header=T)
f2=f2[f2$sp!="droSim1",]
colnames(f2)[24]="branch"
f2=merge(f2,branch_name,by="branch")

f3=fread("/data/correlation/GL_events.txt",header=T)
f3=merge(f3,branch_name,by="branch2")

f=merge(f1,f2,by=c("branch","branch2"))
f=merge(f,f3,by=c("branch","branch2"))
f$branch2=factor(f$branch2,levels=myorder3)
p=ggplot(f, aes(x=len, y=num)) +
  geom_point(aes(color=branch2,shape=branch2),size=3,alpha=1)+
  scale_shape_manual(values = c(1:10,1:10,1:7))+
  #geom_text(aes(label = branch2))+
  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
  labs(x = 'length of 5\' UTR', y = 'Number of uATGs', color = NULL)+
  guides(color = guide_legend(override.aes = list(shape = c(1:10,1:10,1:7))),shape="none") +
  #geom_smooth(method="lm",formula=y~x,color="black")+
  theme_classic( ) #fig_num_len,6.8*4.6
pdf("/results/figS3a_num_len_27sp.pdf",width=6.8,height=4.6)
print(p)
dev.off()
cor.test(f$len,f$num,method = "spearman")
mean(f$num/f$len) #0.007744674
min(f$num/f$len) #0.005878932
max(f$num/f$len) #0.009464298

#num by len
f$num_by_len=f$num/f$len
p=ggplot(f, aes(x=branch2, y=num_by_len)) +
  geom_bar(stat="identity",fill="gray",position = position_dodge(),width=0.8)+
  coord_flip()+
  labs(x = '', y = 'uATG number/length of 5UTR',color = NULL)+
  theme(axis.text.y = element_text(face = "italic")) +
  theme_classic()#fig_num_by_len,4*4.58
pdf("/results/figs3b_num_by_len.pdf",width=4,height=4.6)
print(p)
dev.off()



#only drosophila
f=f[!(f$branch2 %in% c("M.dom","A.gam","A.mel","T.cas")),]
f$branch2=factor(f$branch2,levels=myorder4)

#len vs uATG number
#p=ggplot(f, aes(x=len, y=num)) +
#  geom_point(aes(color=branch2,shape=branch2),size=3,alpha=1)+
#  scale_shape_manual(values = c(1:10,1:10,1:3))+
#  #geom_text(aes(label = branch2))+
#  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
#  labs(x = 'length of 5\' UTR', y = 'Number of uATGs', color = NULL)+
#  guides(color = guide_legend(override.aes = list(shape = c(1:10,1:10,1:3))),shape="none") +
#  #geom_smooth(method="lm",formula=y~x,color="black")+
#  theme_classic( ) #fig_num_len,6.8*4.6
#pdf("/results/fig_num_len_23sp.pdf",width=6.8,height=4.6)
#print(p)
#dev.off()

#cor.test(log10(GL_num_gene2$gain+1),log10(GL_num_gene2$loss+1),method = "pearson")
#cor.test(f$len,f$num,method = "spearman")
#	Spearman's rank correlation rho
#data:  f$len and f$num
#S = 20, p-value = 1.677e-06
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.9901186 



#num vs GC number
f$GC_content=(f$C+f$G)/(f$C+f$G+f$A+f$T)
p=ggplot(f, aes(x=GC_content, y=num/len)) +
  geom_point(aes(color=branch2,shape=branch2),size=3)+
  scale_shape_manual(values = c(1:10,1:10,1:3))+
  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
  labs(x = 'GC content', y = 'Number of uATGs/length of 5UTR', color = NULL)+
  #geom_abline(intercept = 0, slope = 1,color="red",linetype="dashed", size=0.8)+
  guides(color = guide_legend(override.aes = list(shape = c(1:10,1:10,1:3))),shape="none") +
  theme_classic( ) #figS14_num_GCcontent_23sp,6.8*4.6
pdf("/results/figS14_num_GCcontent_23sp.pdf",width=6.8,height=4.6)
print(p)
dev.off()

#cor.test(log10(GL_num_gene2$gain+1),log10(GL_num_gene2$loss+1),method = "pearson")
cor.test(f$GC_content,f$num/f$len,method = "spearman")
#	Spearman's rank correlation rho
#data:  f$GC_content and f$num
#S = 1780, p-value = 0.07183
#alternative hypothesis: true rho is not equal to 0
#sample estimates: rho -0.3833992 


#num vs GC number
f$A_content=(f$A)/(f$C+f$G+f$A+f$T)
f$T_content=(f$T)/(f$C+f$G+f$A+f$T)
f$C_content=(f$C)/(f$C+f$G+f$A+f$T)
f$G_content=(f$G)/(f$C+f$G+f$A+f$T)

#ggplot(f, aes(x=A_content, y=num/len,color=branch2)) +
#  geom_point(size=3,shape=16,fill="black",alpha=0.9)+
#  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
#  labs(x = 'A content', y = 'Number of uATGs/length of 5UTR', color = NULL)+
#  #geom_abline(intercept = 0, slope = 1,color="red",linetype="dashed", size=0.8)+
#  theme_classic( ) #fig_num_Acontent,6.8*4.6
p1=ggplot(f, aes(x=A_content, y=num/len)) +
  geom_point(aes(color=branch2,shape=branch2),size=3)+
  scale_shape_manual(values = c(1:10,1:10,1:3))+
  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
  labs(x = 'A content', y = 'Number of uATGs/length of 5UTR', color = NULL)+
  guides(color = guide_legend(override.aes = list(shape = c(1:10,1:10,1:3))),shape="none") +
  theme_classic( ) #fig_num_Acontent_23sp,6.8*4.6
cor.test(f$A_content,f$num/f$len,method = "spearman") #p-value = 0.1908,rho=0.2826087

p2=ggplot(f, aes(x=T_content, y=num/len)) +
  geom_point(aes(color=branch2,shape=branch2),size=3)+
  scale_shape_manual(values = c(1:10,1:10,1:3))+
  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
  labs(x = 'T content', y = 'Number of uATGs/length of 5UTR', color = NULL)+
  guides(color = guide_legend(override.aes = list(shape = c(1:10,1:10,1:3))),shape="none") +
  theme_classic( ) #fig_num_Tcontent_23sp,6.8*4.6
cor.test(f$T_content,f$num/f$len,method = "spearman") #p-value = 0.03985,rho=0.4337945

p3=ggplot(f, aes(x=C_content, y=num/len)) +
  geom_point(aes(color=branch2,shape=branch2),size=3)+
  scale_shape_manual(values = c(1:10,1:10,1:3))+
  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
  labs(x = 'C content', y = 'Number of uATGs/length of 5UTR', color = NULL)+
  guides(color = guide_legend(override.aes = list(shape = c(1:10,1:10,1:3))),shape="none") +
  theme_classic( ) #fig_num_Ccontent_23sp,6.8*4.6
cor.test(f$C_content,f$num/f$len,method = "spearman") #p-value = 0.07183,rho=-0.3833992 

p4=ggplot(f, aes(x=G_content, y=num/len)) +
  geom_point(aes(color=branch2,shape=branch2),size=3)+
  scale_shape_manual(values = c(1:10,1:10,1:3))+
  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
  labs(x = 'G content', y = 'Number of uATGs/length of 5UTR', color = NULL)+
  guides(color = guide_legend(override.aes = list(shape = c(1:10,1:10,1:3))),shape="none") +
  theme_classic( ) #fig_num_Gcontent_23sp,6.8*4.6
cor.test(f$G_content,f$num/f$len,method = "spearman") #p-value = 0.1107,rho=-0.3418972


pdf("/results/figS13_num_ATCGcontent_23sp.pdf",width=16,height=8)
multiplot(p1, p2, p3, p4, cols=2) ##fig_num_ATCGcontent_23sp, 8*8
dev.off()

#gain_ratio vs. GC
#p=ggplot(f, aes(x=GC_content, y=gain_raio)) +
#  geom_point(aes(color=branch2,shape=branch2),size=3)+
#  scale_shape_manual(values = c(1:10,1:10,1:3))+
#  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
#  labs(x = 'GC content', y = 'gain/loss', color = NULL)+
#  guides(color = guide_legend(override.aes = list(shape = c(1:10,1:10,1:3))),#shape="none") +
#  theme_classic( ) #fig_gainratio_GCcontent_23sp,6.8*4.6
#pdf("/results/fig_gainratio_GCcontent_23sp.pdf",width=6.8,height=4.6)
#print(p)
#dev.off()
#cor.test(f$gain_raio,f$GC_content,method = "spearman") #p-value = 0.006863,rho=-0.5543478
#cor.test(f$gain/f$loss,f$GC_content,method = "spearman") #p-value = 0.006863,rho=-0.5543478

#ggplot(f[f$gain_loss_ratio<40,], aes(x=GC_content, y=gain_loss_ratio)) +
#  geom_point(aes(color=branch2,shape=branch2),size=3)+
#  scale_shape_manual(values = c(1:10,1:10))+
#  #coord_cartesian(xlim = c(0, 3),ylim = c(0, 3))+ 
#  labs(x = 'GC content', y = 'gain/loss', color = NULL)+
#  guides(color = guide_legend(override.aes = list(shape = c(1:10,1:10,1:1))),shape="none") +
#  theme_classic( ) #fig_gainratio_GCcontent_23sp,6.8*4.6

#dinuc
f$total_dinu=rowSums(f[,11:26])

dinu=colnames(f)[11:26]
dinu_out <- matrix(ncol=16,nrow=2)
colnames(dinu_out) <- dinu
rownames(dinu_out) <- c("rho","p")

for (i in dinu){
  selected_columns <- c("num_by_len",i,"total_dinu")
  tmp <- f[,..selected_columns]
  colnames(tmp)=c("num_by_len","di","total_dinu")
  c=cor.test(tmp$num_by_len,tmp$di/tmp$total_dinu,method = "spearman")
  dinu_out["rho",i] <- c$estimate
  dinu_out["p",i] <- c$p.value
}

dinu_out2=as.data.table(t(dinu_out))
dinu_out2$dinu=colnames(dinu_out)
colnames(dinu_out2)=c("rho","p","dinu")

dinu_out2=dinu_out2[order(dinu_out2$rho),]
dinu_out2$dinu=factor(dinu_out2$dinu,levels=dinu_out2$dinu)
dinu_out2$q=p.adjust(dinu_out2$p,method="BH")
dinu_out2$color="grey"
dinu_out2[dinu_out2$q<0.05,]$color="red"

p=ggplot(dinu_out2, aes(x = dinu, y = rho,fill=color)) +
  geom_bar(stat = "identity", position = "dodge",fill=dinu_out2$color,width = 0.7) +
  labs(x = "Dinucleotide", y = "Rho") +
  #scale_fill_manual(values=c("#FF4040"))+
  theme_classic()+
  coord_flip() #merged_dinu_uATGdensity_rho_sp23.pdf #2.82*4.72

pdf("/results/figS15_merged_dinu_uATGdensity_rho_sp23.pdf",width=2.82,height=4.72)
print(p)
dev.off()
#fwrite(dinu_out2,"dinu_uATGdensity_rho_p_sp23.txt")

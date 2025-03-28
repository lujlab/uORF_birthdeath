library(data.table)
setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/12.pos_to_cds/2.compensation")
#canonical=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/from_dyg/dmel-all-r6.04.canonical.transcripts.txt",header=F,sep='\t')
uORF_matrix_canonical=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/8.new_27_mfa/PacBioSim/position_cds/uORF_matrix_dm6_PacBioSim_27species_update_noCDS_postoCDS_canonical.csv",sep=',',header=T)
h=uORF_matrix_canonical[,c(1:4,7,11,31:37)]
h=h[h$dm6+h$PacBioSim+h$droYak3+h$droAna3>0,]

#1.stat ratio of GL events
#####4bin stat#####
h2=h
h2$distance_class="0-100nt"
h2[h2$pos_to_cds>100,]$distance_class="100-300nt"
h2[h2$pos_to_cds>300,]$distance_class="300-500nt"
h2[h2$pos_to_cds>500,]$distance_class=">500nt"

l=c("0-100nt","100-300nt","300-500nt",">500nt")
x=matrix(data=NA,nrow=9,ncol=8)
for (i in c(1:length(l))){
  t=h2[h2$distance_class==l[i],]
  #++++,ana yak sim mel
  z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
  x[1,i]=length(z$ID)
  x[1,i+4]=length(unique(z$geneID))
  #+-++
  z=t[t$droAna3==1 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
  x[2,i]=length(z$ID)
  x[2,i+4]=length(unique(z$geneID))
  #-+++
  z=t[t$droAna3==0 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
  x[3,i]=length(z$ID)
  x[3,i+4]=length(unique(z$geneID))
  #--++
  z=t[t$droAna3==0 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
  x[4,i]=length(z$ID)
  x[4,i+4]=length(unique(z$geneID))
  #++--
  z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==0 & t$dm6==0,]
  x[5,i]=length(z$ID)
  x[5,i+4]=length(unique(z$geneID))
  #N--+
  z=t[t$droYak3==0 & t$PacBioSim==0 & t$dm6==1,]
  x[6,i]=length(z$ID)
  x[6,i+4]=length(unique(z$geneID))
  #N++-
  z=t[t$droYak3==1 & t$PacBioSim==1 & t$dm6==0,]
  x[7,i]=length(z$ID)
  x[7,i+4]=length(unique(z$geneID))
  #N-+-
  z=t[t$droYak3==0 & t$PacBioSim==1 & t$dm6==0,]
  x[8,i]=length(z$ID)
  x[8,i+4]=length(unique(z$geneID))
  #N+-+
  z=t[t$droYak3==1 & t$PacBioSim==0 & t$dm6==1,]
  x[9,i]=length(z$ID)
  x[9,i+4]=length(unique(z$geneID))
}
y=data.table(x)
y$type=c("++++","+-++","-+++","--++","++--","N--+","N++-","N-+-","N+-+")
colnames(y)=c("0-100nt_uORF","100-300nt_uORF","300-500nt_uORF",">500nt_uORF","0-100nt_gene","100-300nt_gene","300-500nt_gene",">500nt_gene","type")
y$all_uORF=0
y$all_gene=0
t=h2
#++++,ana yak sim mel
z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[1]=length(z$ID)
y$all_gene[1]=length(unique(z$geneID))
#+-++
z=t[t$droAna3==1 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[2]=length(z$ID)
y$all_gene[2]=length(unique(z$geneID))
#-+++
z=t[t$droAna3==0 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[3]=length(z$ID)
y$all_gene[3]=length(unique(z$geneID))
#--++
z=t[t$droAna3==0 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[4]=length(z$ID)
y$all_gene[4]=length(unique(z$geneID))
#++--
z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==0 & t$dm6==0,]
y$all_uORF[5]=length(z$ID)
y$all_gene[5]=length(unique(z$geneID))
#N--+
z=t[t$droYak3==0 & t$PacBioSim==0 & t$dm6==1,]
y$all_uORF[6]=length(z$ID)
y$all_gene[6]=length(unique(z$geneID))
#N++-
z=t[t$droYak3==1 & t$PacBioSim==1 & t$dm6==0,]
y$all_uORF[7]=length(z$ID)
y$all_gene[7]=length(unique(z$geneID))
#N-+-
z=t[t$droYak3==0 & t$PacBioSim==1 & t$dm6==0,]
y$all_uORF[8]=length(z$ID)
y$all_gene[8]=length(unique(z$geneID))
#N+-+
z=t[t$droYak3==1 & t$PacBioSim==0 & t$dm6==1,]
y$all_uORF[9]=length(z$ID)
y$all_gene[9]=length(unique(z$geneID))

y2=y[,c("type","all_uORF","0-100nt_uORF","100-300nt_uORF","300-500nt_uORF",">500nt_uORF","all_gene","0-100nt_gene","100-300nt_gene","300-500nt_gene",">500nt_gene")]
fwrite(y2,"class_uORF_gene_num_DisToCDS.txt")
y2=fread("class_uORF_gene_num_DisToCDS.txt")
y3=t(y2)

colnames(y3)=c(y3[1,])
y3=y3[-1,]
y3=apply(y3,2,as.numeric)
y4=y3
y4=data.table(y4)
y4$class=c("all_uORF","0-100nt_uORF","100-300nt_uORF","300-500nt_uORF",">500nt_uORF","all_gene","0-100nt_gene","100-300nt_gene","300-500nt_gene",">500nt_gene")

y4$gain_before_b3=y4$`++++`+y4$`+-++`+y4$`-+++`
y4$gain_in_b3=y4$`--++`
y4$lost_in_b3=y4$`++--`
y4$gain_in_b1=y4$`N--+`
y4$lost_in_b1=y4$`N++-`
y4$gain_in_b2=y4$`N-+-`
y4$lost_in_b2=y4$`N+-+`
y4$melsim_cons=y4$`++++`+y4$`+-++`+y4$`-+++`+y4$`--++`
y4$mel_specific=y4$`N--+`+y4$`N+-+`
y4$sim_specific=y4$`N-+-`+y4$`N++-`
y4=y4[,10:20]

y4$gain_in_b3_ratio=y4$gain_in_b3/y4$gain_before_b3
y4$lost_in_b3_ratio=y4$lost_in_b3/y4$gain_before_b3
y4$gain_in_b1_ratio=y4$gain_in_b1/y4$gain_before_b3
y4$lost_in_b1_ratio=y4$lost_in_b1/y4$gain_before_b3
y4$gain_in_b2_ratio=y4$gain_in_b2/y4$gain_before_b3
y4$lost_in_b2_ratio=y4$lost_in_b2/y4$gain_before_b3

y4$`lostb1-b3_ratio`=y4$lost_in_b3_ratio+y4$lost_in_b1_ratio+y4$lost_in_b2_ratio
y4$`gainb1-b3_ratio`=y4$gain_in_b3_ratio+y4$gain_in_b1_ratio+y4$gain_in_b2_ratio
y4$`GLb1-b3_ratio`=y4$`lostb1-b3_ratio`+y4$`gainb1-b3_ratio`

y4$`lostb1-b3`=y4$lost_in_b3+y4$lost_in_b1+y4$lost_in_b2
y4$`gainb1-b3`=y4$gain_in_b3+y4$gain_in_b1+y4$gain_in_b2
y4$`GLb1-b3`=y4$`lostb1-b3`+y4$`gainb1-b3`

y4$mel_specific_ratio=y4$mel_specific/y4$melsim_cons
y4$sim_specific_ratio=y4$sim_specific/y4$melsim_cons
y4$specific_ratio=(y4$sim_specific+y4$mel_specific)/y4$melsim_cons

y5=y4[2:5,]
y5$class=factor(y5$class,levels=y5$class)
ggplot(data = y5, mapping = aes(x = class, y = specific_ratio,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("specific(sim+mel)/conserved")+
  theme_classic() #SpecificToConserved_DisToCDS_4bins.pdf
ggplot(data = y5, mapping = aes(x = class, y = `lostb1-b3_ratio`,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("lostb1-b3_ratio/gain_before_b3")+
  theme_classic() #lostb1-b3_DisToCDS_4bins.pdf
ggplot(data = y5, mapping = aes(x = class, y = `gainb1-b3_ratio`,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("gainb1-b3_ratio/gain_before_b3")+
  theme_classic() #gainb1-b3_ratio_DisToCDS_4bins.pdf
ggplot(data = y5, mapping = aes(x = class, y = `GLb1-b3_ratio`,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("GLb1-b3_ratio/gain_before_b3")+
  theme_classic() #GLb1-b3_ratio_DisToCDS_4bins.pdf


#specific to conserved
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="100-300nt_uORF",]$melsim_cons,y5[y5$class=="100-300nt_uORF",]$mel_specific+y5[y5$class=="100-300nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="300-500nt_uORF",]$melsim_cons,y5[y5$class=="300-500nt_uORF",]$mel_specific+y5[y5$class=="300-500nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class==">500nt_uORF",]$melsim_cons,y5[y5$class==">500nt_uORF",]$mel_specific+y5[y5$class==">500nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="100-300nt_uORF",]$melsim_cons,y5[y5$class=="100-300nt_uORF",]$mel_specific+y5[y5$class=="100-300nt_uORF",]$sim_specific,
                     y5[y5$class=="300-500nt_uORF",]$melsim_cons,y5[y5$class=="300-500nt_uORF",]$mel_specific+y5[y5$class=="300-500nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="100-300nt_uORF",]$melsim_cons,y5[y5$class=="100-300nt_uORF",]$mel_specific+y5[y5$class=="100-300nt_uORF",]$sim_specific,
                     y5[y5$class==">500nt_uORF",]$melsim_cons,y5[y5$class==">500nt_uORF",]$mel_specific+y5[y5$class==">500nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="300-500nt_uORF",]$melsim_cons,y5[y5$class=="300-500nt_uORF",]$mel_specific+y5[y5$class=="300-500nt_uORF",]$sim_specific,
                     y5[y5$class==">500nt_uORF",]$melsim_cons,y5[y5$class==">500nt_uORF",]$mel_specific+y5[y5$class==">500nt_uORF",]$sim_specific),nrow=2))

#lost
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$gain_before_b3,y5[y5$class=="0-100nt_uORF",]$`lostb1-b3`,
                     y5[y5$class=="100-300nt_uORF",]$gain_before_b3,y5[y5$class=="100-300nt_uORF",]$`lostb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$gain_before_b3,y5[y5$class=="0-100nt_uORF",]$`lostb1-b3`,
                     y5[y5$class=="300-500nt_uORF",]$gain_before_b3,y5[y5$class=="300-500nt_uORF",]$`lostb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$gain_before_b3,y5[y5$class=="0-100nt_uORF",]$`lostb1-b3`,
                     y5[y5$class==">500nt_uORF",]$gain_before_b3,y5[y5$class==">500nt_uORF",]$`lostb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="100-300nt_uORF",]$gain_before_b3,y5[y5$class=="100-300nt_uORF",]$`lostb1-b3`,
                     y5[y5$class=="300-500nt_uORF",]$gain_before_b3,y5[y5$class=="300-500nt_uORF",]$`lostb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="100-300nt_uORF",]$gain_before_b3,y5[y5$class=="100-300nt_uORF",]$`lostb1-b3`,
                     y5[y5$class==">500nt_uORF",]$gain_before_b3,y5[y5$class==">500nt_uORF",]$`lostb1-b3`),nrow=2))
#gain
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$gain_before_b3,y5[y5$class=="0-100nt_uORF",]$`gainb1-b3`,
                     y5[y5$class=="100-300nt_uORF",]$gain_before_b3,y5[y5$class=="100-300nt_uORF",]$`gainb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$gain_before_b3,y5[y5$class=="0-100nt_uORF",]$`gainb1-b3`,
                     y5[y5$class=="300-500nt_uORF",]$gain_before_b3,y5[y5$class=="300-500nt_uORF",]$`gainb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$gain_before_b3,y5[y5$class=="0-100nt_uORF",]$`gainb1-b3`,
                     y5[y5$class==">500nt_uORF",]$gain_before_b3,y5[y5$class==">500nt_uORF",]$`gainb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="100-300nt_uORF",]$gain_before_b3,y5[y5$class=="100-300nt_uORF",]$`gainb1-b3`,
                     y5[y5$class=="300-500nt_uORF",]$gain_before_b3,y5[y5$class=="300-500nt_uORF",]$`gainb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="100-300nt_uORF",]$gain_before_b3,y5[y5$class=="100-300nt_uORF",]$`gainb1-b3`,
                     y5[y5$class==">500nt_uORF",]$gain_before_b3,y5[y5$class==">500nt_uORF",]$`gainb1-b3`),nrow=2))


#GL
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$gain_before_b3,y5[y5$class=="0-100nt_uORF",]$`GLb1-b3`,
                     y5[y5$class=="100-300nt_uORF",]$gain_before_b3,y5[y5$class=="100-300nt_uORF",]$`GLb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$gain_before_b3,y5[y5$class=="0-100nt_uORF",]$`GLb1-b3`,
                     y5[y5$class=="300-500nt_uORF",]$gain_before_b3,y5[y5$class=="300-500nt_uORF",]$`GLb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$gain_before_b3,y5[y5$class=="0-100nt_uORF",]$`GLb1-b3`,
                     y5[y5$class==">500nt_uORF",]$gain_before_b3,y5[y5$class==">500nt_uORF",]$`GLb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="100-300nt_uORF",]$gain_before_b3,y5[y5$class=="100-300nt_uORF",]$`GLb1-b3`,
                     y5[y5$class=="300-500nt_uORF",]$gain_before_b3,y5[y5$class=="300-500nt_uORF",]$`GLb1-b3`),nrow=2))
fisher.test(matrix(c(y5[y5$class=="100-300nt_uORF",]$gain_before_b3,y5[y5$class=="100-300nt_uORF",]$`GLb1-b3`,
                     y5[y5$class==">500nt_uORF",]$gain_before_b3,y5[y5$class==">500nt_uORF",]$`GLb1-b3`),nrow=2))















####8bin####
h2=h
h2$distance_class="0-100nt"
h2[h2$pos_to_cds>100,]$distance_class="100-200nt"
h2[h2$pos_to_cds>200,]$distance_class="200-300nt"
h2[h2$pos_to_cds>300,]$distance_class="300-400nt"
h2[h2$pos_to_cds>400,]$distance_class="400-500nt"
h2[h2$pos_to_cds>500,]$distance_class="500-600nt"
h2[h2$pos_to_cds>600,]$distance_class="600-700nt"
h2[h2$pos_to_cds>700,]$distance_class=">700nt"


l=c("0-100nt","100-200nt","200-300nt","300-400nt","400-500nt","500-600nt","600-700nt",">700nt")
x=matrix(data=NA,nrow=9,ncol=16)
for (i in c(1:length(l))){
  t=h2[h2$distance_class==l[i],]
  #++++,ana yak sim mel
  z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
  x[1,i]=length(z$ID)
  x[1,i+8]=length(unique(z$geneID))
  #+-++
  z=t[t$droAna3==1 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
  x[2,i]=length(z$ID)
  x[2,i+8]=length(unique(z$geneID))
  #-+++
  z=t[t$droAna3==0 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
  x[3,i]=length(z$ID)
  x[3,i+8]=length(unique(z$geneID))
  #--++
  z=t[t$droAna3==0 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
  x[4,i]=length(z$ID)
  x[4,i+8]=length(unique(z$geneID))
  #++--
  z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==0 & t$dm6==0,]
  x[5,i]=length(z$ID)
  x[5,i+8]=length(unique(z$geneID))
  #N--+
  z=t[t$droYak3==0 & t$PacBioSim==0 & t$dm6==1,]
  x[6,i]=length(z$ID)
  x[6,i+8]=length(unique(z$geneID))
  #N++-
  z=t[t$droYak3==1 & t$PacBioSim==1 & t$dm6==0,]
  x[7,i]=length(z$ID)
  x[7,i+8]=length(unique(z$geneID))
  #N-+-
  z=t[t$droYak3==0 & t$PacBioSim==1 & t$dm6==0,]
  x[8,i]=length(z$ID)
  x[8,i+8]=length(unique(z$geneID))
  #N+-+
  z=t[t$droYak3==1 & t$PacBioSim==0 & t$dm6==1,]
  x[9,i]=length(z$ID)
  x[9,i+8]=length(unique(z$geneID))
}
y=data.table(x)
y$type=c("++++","+-++","-+++","--++","++--","N--+","N++-","N-+-","N+-+")
colnames(y)=c("0-100nt_uORF","100-200nt_uORF","200-300nt_uORF","300-400nt_uORF","400-500nt_uORF","500-600nt_uORF","600-700nt_uORF",">700nt_uORF","0-100nt_gene","100-200nt_gene","200-300nt_gene","300-400nt_gene","400-500nt_gene","500-600nt_gene","600-700nt_gene",">700nt_gene","type")
y$all_uORF=0
y$all_gene=0
t=h2
#++++,ana yak sim mel
z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[1]=length(z$ID)
y$all_gene[1]=length(unique(z$geneID))
#+-++
z=t[t$droAna3==1 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[2]=length(z$ID)
y$all_gene[2]=length(unique(z$geneID))
#-+++
z=t[t$droAna3==0 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[3]=length(z$ID)
y$all_gene[3]=length(unique(z$geneID))
#--++
z=t[t$droAna3==0 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[4]=length(z$ID)
y$all_gene[4]=length(unique(z$geneID))
#++--
z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==0 & t$dm6==0,]
y$all_uORF[5]=length(z$ID)
y$all_gene[5]=length(unique(z$geneID))
#N--+
z=t[t$droYak3==0 & t$PacBioSim==0 & t$dm6==1,]
y$all_uORF[6]=length(z$ID)
y$all_gene[6]=length(unique(z$geneID))
#N++-
z=t[t$droYak3==1 & t$PacBioSim==1 & t$dm6==0,]
y$all_uORF[7]=length(z$ID)
y$all_gene[7]=length(unique(z$geneID))
#N-+-
z=t[t$droYak3==0 & t$PacBioSim==1 & t$dm6==0,]
y$all_uORF[8]=length(z$ID)
y$all_gene[8]=length(unique(z$geneID))
#N+-+
z=t[t$droYak3==1 & t$PacBioSim==0 & t$dm6==1,]
y$all_uORF[9]=length(z$ID)
y$all_gene[9]=length(unique(z$geneID))

y2=y[,c("type","all_uORF","0-100nt_uORF","100-200nt_uORF","200-300nt_uORF","300-400nt_uORF","400-500nt_uORF","500-600nt_uORF","600-700nt_uORF",">700nt_uORF","all_gene","0-100nt_gene","100-200nt_gene","200-300nt_gene","300-400nt_gene","400-500nt_gene","500-600nt_gene","600-700nt_gene",">700nt_gene")]
fwrite(y2,"class_uORF_gene_num_DisToCDS_8bins.txt")
y2=fread("class_uORF_gene_num_DisToCDS_8bins.txt")
y3=t(y2)

colnames(y3)=c(y3[1,])
y3=y3[-1,]
y3=apply(y3,2,as.numeric)
y4=y3
y4=data.table(y4)
y4$class=c("all_uORF","0-100nt_uORF","100-200nt_uORF","200-300nt_uORF","300-400nt_uORF","400-500nt_uORF","500-600nt_uORF","600-700nt_uORF",">700nt_uORF","all_gene","0-100nt_gene","100-200nt_gene","200-300nt_gene","300-400nt_gene","400-500nt_gene","500-600nt_gene","600-700nt_gene",">700nt_gene")

y4$gain_before_b3=y4$`++++`+y4$`+-++`+y4$`-+++`
y4$gain_in_b3=y4$`--++`
y4$lost_in_b3=y4$`++--`
y4$gain_in_b1=y4$`N--+`
y4$lost_in_b1=y4$`N++-`
y4$gain_in_b2=y4$`N-+-`
y4$lost_in_b2=y4$`N+-+`
y4$melsim_cons=y4$`++++`+y4$`+-++`+y4$`-+++`+y4$`--++`
y4$mel_specific=y4$`N--+`+y4$`N+-+`
y4$sim_specific=y4$`N-+-`+y4$`N++-`
y4=y4[,10:20]

y4$gain_in_b3_ratio=y4$gain_in_b3/y4$gain_before_b3
y4$lost_in_b3_ratio=y4$lost_in_b3/y4$gain_before_b3
y4$gain_in_b1_ratio=y4$gain_in_b1/y4$gain_before_b3
y4$lost_in_b1_ratio=y4$lost_in_b1/y4$gain_before_b3
y4$gain_in_b2_ratio=y4$gain_in_b2/y4$gain_before_b3
y4$lost_in_b2_ratio=y4$lost_in_b2/y4$gain_before_b3

y4$`lostb1-b3_ratio`=y4$lost_in_b3_ratio+y4$lost_in_b1_ratio+y4$lost_in_b2_ratio
y4$`gainb1-b3_ratio`=y4$gain_in_b3_ratio+y4$gain_in_b1_ratio+y4$gain_in_b2_ratio
y4$`GLb1-b3_ratio`=y4$`lostb1-b3_ratio`+y4$`gainb1-b3_ratio`

y4$`lostb1-b3`=y4$lost_in_b3+y4$lost_in_b1+y4$lost_in_b2
y4$`gainb1-b3`=y4$gain_in_b3+y4$gain_in_b1+y4$gain_in_b2
y4$`GLb1-b3`=y4$`lostb1-b3`+y4$`gainb1-b3`

y4$mel_specific_ratio=y4$mel_specific/y4$melsim_cons
y4$sim_specific_ratio=y4$sim_specific/y4$melsim_cons
y4$specific_ratio=(y4$sim_specific+y4$mel_specific)/y4$melsim_cons

y5=y4[2:9,]
y5$class=factor(y5$class,levels=y5$class)
ggplot(data = y5, mapping = aes(x = class, y = specific_ratio,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("specific(sim+mel)/conserved")+
  theme_classic() #SpecificToConserved_DisToCDS_8bins.pdf
ggplot(data = y5, mapping = aes(x = class, y = `lostb1-b3_ratio`,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("lostb1-b3_ratio/gain_before_b3")+
  theme_classic() #lostb1-b3_DisToCDS_8bins.pdf
ggplot(data = y5, mapping = aes(x = class, y = `gainb1-b3_ratio`,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("gainb1-b3_ratio/gain_before_b3")+
  theme_classic() #gainb1-b3_ratio_DisToCDS_8bins.pdf
ggplot(data = y5, mapping = aes(x = class, y = `GLb1-b3_ratio`,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("GLb1-b3_ratio/gain_before_b3")+
  theme_classic() #GLb1-b3_ratio_DisToCDS_8bins.pdf


#specific to conserved
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="100-200nt_uORF",]$melsim_cons,y5[y5$class=="100-200nt_uORF",]$mel_specific+y5[y5$class=="100-200nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="200-300nt_uORF",]$melsim_cons,y5[y5$class=="200-300nt_uORF",]$mel_specific+y5[y5$class=="200-300nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="300-400nt_uORF",]$melsim_cons,y5[y5$class=="300-400nt_uORF",]$mel_specific+y5[y5$class=="300-400nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="400-500nt_uORF",]$melsim_cons,y5[y5$class=="400-500nt_uORF",]$mel_specific+y5[y5$class=="400-500nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="500-600nt_uORF",]$melsim_cons,y5[y5$class=="500-600nt_uORF",]$mel_specific+y5[y5$class=="500-600nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="600-700nt_uORF",]$melsim_cons,y5[y5$class=="600-700nt_uORF",]$mel_specific+y5[y5$class=="600-700nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class==">700nt_uORF",]$melsim_cons,y5[y5$class==">700nt_uORF",]$mel_specific+y5[y5$class==">700nt_uORF",]$sim_specific),nrow=2))


















########3bin#####
h2=h
h2$distance_class="0-200nt"
h2[h2$pos_to_cds>200,]$distance_class=">200nt"


l=c("0-200nt",">200nt")
x=matrix(data=NA,nrow=9,ncol=4)
for (i in c(1:length(l))){
  t=h2[h2$distance_class==l[i],]
  #++++,ana yak sim mel
  z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
  x[1,i]=length(z$ID)
  x[1,i+2]=length(unique(z$geneID))
  #+-++
  z=t[t$droAna3==1 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
  x[2,i]=length(z$ID)
  x[2,i+2]=length(unique(z$geneID))
  #-+++
  z=t[t$droAna3==0 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
  x[3,i]=length(z$ID)
  x[3,i+2]=length(unique(z$geneID))
  #--++
  z=t[t$droAna3==0 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
  x[4,i]=length(z$ID)
  x[4,i+2]=length(unique(z$geneID))
  #++--
  z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==0 & t$dm6==0,]
  x[5,i]=length(z$ID)
  x[5,i+2]=length(unique(z$geneID))
  #N--+
  z=t[t$droYak3==0 & t$PacBioSim==0 & t$dm6==1,]
  x[6,i]=length(z$ID)
  x[6,i+2]=length(unique(z$geneID))
  #N++-
  z=t[t$droYak3==1 & t$PacBioSim==1 & t$dm6==0,]
  x[7,i]=length(z$ID)
  x[7,i+2]=length(unique(z$geneID))
  #N-+-
  z=t[t$droYak3==0 & t$PacBioSim==1 & t$dm6==0,]
  x[8,i]=length(z$ID)
  x[8,i+2]=length(unique(z$geneID))
  #N+-+
  z=t[t$droYak3==1 & t$PacBioSim==0 & t$dm6==1,]
  x[9,i]=length(z$ID)
  x[9,i+2]=length(unique(z$geneID))
}
y=data.table(x)
y$type=c("++++","+-++","-+++","--++","++--","N--+","N++-","N-+-","N+-+")
colnames(y)=c("0-200nt_uORF",">200nt_uORF","0-200nt_gene",">200nt_gene","type")
y$all_uORF=0
y$all_gene=0
t=h2
#++++,ana yak sim mel
z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[1]=length(z$ID)
y$all_gene[1]=length(unique(z$geneID))
#+-++
z=t[t$droAna3==1 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[2]=length(z$ID)
y$all_gene[2]=length(unique(z$geneID))
#-+++
z=t[t$droAna3==0 & t$droYak3==1 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[3]=length(z$ID)
y$all_gene[3]=length(unique(z$geneID))
#--++
z=t[t$droAna3==0 & t$droYak3==0 & t$PacBioSim==1 & t$dm6==1,]
y$all_uORF[4]=length(z$ID)
y$all_gene[4]=length(unique(z$geneID))
#++--
z=t[t$droAna3==1 & t$droYak3==1 & t$PacBioSim==0 & t$dm6==0,]
y$all_uORF[5]=length(z$ID)
y$all_gene[5]=length(unique(z$geneID))
#N--+
z=t[t$droYak3==0 & t$PacBioSim==0 & t$dm6==1,]
y$all_uORF[6]=length(z$ID)
y$all_gene[6]=length(unique(z$geneID))
#N++-
z=t[t$droYak3==1 & t$PacBioSim==1 & t$dm6==0,]
y$all_uORF[7]=length(z$ID)
y$all_gene[7]=length(unique(z$geneID))
#N-+-
z=t[t$droYak3==0 & t$PacBioSim==1 & t$dm6==0,]
y$all_uORF[8]=length(z$ID)
y$all_gene[8]=length(unique(z$geneID))
#N+-+
z=t[t$droYak3==1 & t$PacBioSim==0 & t$dm6==1,]
y$all_uORF[9]=length(z$ID)
y$all_gene[9]=length(unique(z$geneID))

y2=y[,c("type","all_uORF","0-200nt_uORF",">200nt_uORF","all_gene","0-200nt_gene",">200nt_gene")]
#fwrite(y2,"class_uORF_gene_num_DisToCDS_3bins.txt")
#y2=fread("class_uORF_gene_num_DisToCDS_3bins.txt")
y3=t(y2)

colnames(y3)=c(y3[1,])
y3=y3[-1,]
y3=apply(y3,2,as.numeric)
y4=y3
y4=data.table(y4)
y4$class=c("all_uORF","0-200nt_uORF",">200nt_uORF","all_gene","0-200nt_gene",">200nt_gene")

y4$gain_before_b3=y4$`++++`+y4$`+-++`+y4$`-+++`
y4$gain_in_b3=y4$`--++`
y4$lost_in_b3=y4$`++--`
y4$gain_in_b1=y4$`N--+`
y4$lost_in_b1=y4$`N++-`
y4$gain_in_b2=y4$`N-+-`
y4$lost_in_b2=y4$`N+-+`
y4$melsim_cons=y4$`++++`+y4$`+-++`+y4$`-+++`+y4$`--++`
y4$mel_specific=y4$`N--+`+y4$`N+-+`
y4$sim_specific=y4$`N-+-`+y4$`N++-`
y4=y4[,10:20]

y4$gain_in_b3_ratio=y4$gain_in_b3/y4$gain_before_b3
y4$lost_in_b3_ratio=y4$lost_in_b3/y4$gain_before_b3
y4$gain_in_b1_ratio=y4$gain_in_b1/y4$gain_before_b3
y4$lost_in_b1_ratio=y4$lost_in_b1/y4$gain_before_b3
y4$gain_in_b2_ratio=y4$gain_in_b2/y4$gain_before_b3
y4$lost_in_b2_ratio=y4$lost_in_b2/y4$gain_before_b3

y4$`lostb1-b3_ratio`=y4$lost_in_b3_ratio+y4$lost_in_b1_ratio+y4$lost_in_b2_ratio
y4$`gainb1-b3_ratio`=y4$gain_in_b3_ratio+y4$gain_in_b1_ratio+y4$gain_in_b2_ratio
y4$`GLb1-b3_ratio`=y4$`lostb1-b3_ratio`+y4$`gainb1-b3_ratio`

y4$`lostb1-b3`=y4$lost_in_b3+y4$lost_in_b1+y4$lost_in_b2
y4$`gainb1-b3`=y4$gain_in_b3+y4$gain_in_b1+y4$gain_in_b2
y4$`GLb1-b3`=y4$`lostb1-b3`+y4$`gainb1-b3`

y4$mel_specific_ratio=y4$mel_specific/y4$melsim_cons
y4$sim_specific_ratio=y4$sim_specific/y4$melsim_cons
y4$specific_ratio=(y4$sim_specific+y4$mel_specific)/y4$melsim_cons

y5=y4[2:3,]
y5$class=factor(y5$class,levels=y5$class)
ggplot(data = y5, mapping = aes(x = class, y = specific_ratio,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("specific(sim+mel)/conserved")+
  theme_classic() #SpecificToConserved_DisToCDS_8bins.pdf
ggplot(data = y5, mapping = aes(x = class, y = `lostb1-b3_ratio`,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("lostb1-b3_ratio/gain_before_b3")+
  theme_classic() #lostb1-b3_DisToCDS_8bins.pdf
ggplot(data = y5, mapping = aes(x = class, y = `gainb1-b3_ratio`,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("gainb1-b3_ratio/gain_before_b3")+
  theme_classic() #gainb1-b3_ratio_DisToCDS_8bins.pdf
ggplot(data = y5, mapping = aes(x = class, y = `GLb1-b3_ratio`,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("GLb1-b3_ratio/gain_before_b3")+
  theme_classic() #GLb1-b3_ratio_DisToCDS_8bins.pdf


#specific to conserved
fisher.test(matrix(c(y5[y5$class=="0-200nt_uORF",]$melsim_cons,y5[y5$class=="0-200nt_uORF",]$mel_specific+y5[y5$class=="0-200nt_uORF",]$sim_specific,
                     y5[y5$class==">200nt_uORF",]$melsim_cons,y5[y5$class==">200nt_uORF",]$mel_specific+y5[y5$class==">200nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="200-300nt_uORF",]$melsim_cons,y5[y5$class=="200-300nt_uORF",]$mel_specific+y5[y5$class=="200-300nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="300-400nt_uORF",]$melsim_cons,y5[y5$class=="300-400nt_uORF",]$mel_specific+y5[y5$class=="300-400nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="400-500nt_uORF",]$melsim_cons,y5[y5$class=="400-500nt_uORF",]$mel_specific+y5[y5$class=="400-500nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="500-600nt_uORF",]$melsim_cons,y5[y5$class=="500-600nt_uORF",]$mel_specific+y5[y5$class=="500-600nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class=="600-700nt_uORF",]$melsim_cons,y5[y5$class=="600-700nt_uORF",]$mel_specific+y5[y5$class=="600-700nt_uORF",]$sim_specific),nrow=2))
fisher.test(matrix(c(y5[y5$class=="0-100nt_uORF",]$melsim_cons,y5[y5$class=="0-100nt_uORF",]$mel_specific+y5[y5$class=="0-100nt_uORF",]$sim_specific,
                     y5[y5$class==">700nt_uORF",]$melsim_cons,y5[y5$class==">700nt_uORF",]$mel_specific+y5[y5$class==">700nt_uORF",]$sim_specific),nrow=2))





########number of uORF in each bin #########
len=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/12.pos_to_cds/all_tr_melsimGL.txt",sep='\t')
len$class="0-100nt_uORF"

len[len$pos_to_cds_class>=2 &len$pos_to_cds_class<=5 ,]$class="100-300nt_uORF"
len[len$pos_to_cds_class>=6 &len$pos_to_cds_class<=9 ,]$class="300-500nt_uORF"
len[len$pos_to_cds_class>=10,]$class=">500nt_uORF"
totallen=aggregate(len$total_len,by=list(len$class),sum) %>% as.data.table
colnames(totallen)=c("class","len")
y6=merge(y5,totallen)
y6$class=factor(y6$class,levels=y6$class)
ggplot(data = y6, mapping = aes(x = class, y = y6$melsim_cons/y6$len,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("melsim_cons/5UTR_len")+
  theme_classic() #melsim_cons_5UTR_len_4bins.pdf
ggplot(data = y6, mapping = aes(x = class, y = y6$gain_before_b3/y6$len,fill=class)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("melsim_cons/5UTR_len")+
  theme_classic() #gain_before_b3_5UTR_len_4bins.pdf

fisher.test(matrix(c(y6[y6$class=="0-100nt_uORF",]$melsim_cons,y6[y6$class=="0-100nt_uORF",]$len,
                     y6[y6$class=="100-300nt_uORF",]$melsim_cons,y6[y6$class=="100-300nt_uORF",]$len),nrow=2))
fisher.test(matrix(c(y6[y6$class=="0-100nt_uORF",]$melsim_cons,y6[y6$class=="0-100nt_uORF",]$len,
                     y6[y6$class=="300-500nt_uORF",]$melsim_cons,y6[y6$class=="300-500nt_uORF",]$len),nrow=2))
fisher.test(matrix(c(y6[y6$class=="0-100nt_uORF",]$melsim_cons,y6[y6$class=="0-100nt_uORF",]$len,
                     y6[y6$class==">500nt_uORF",]$melsim_cons,y6[y6$class==">500nt_uORF",]$len),nrow=2))

fisher.test(matrix(c(y6[y6$class=="100-300nt_uORF",]$melsim_cons,y6[y6$class=="100-300nt_uORF",]$len,
                     y6[y6$class=="300-500nt_uORF",]$melsim_cons,y6[y6$class=="300-500nt_uORF",]$len),nrow=2))
fisher.test(matrix(c(y6[y6$class=="100-300nt_uORF",]$melsim_cons,y6[y6$class=="100-300nt_uORF",]$len,
                     y6[y6$class==">500nt_uORF",]$melsim_cons,y6[y6$class==">500nt_uORF",]$len),nrow=2))
fisher.test(matrix(c(y6[y6$class=="300-500nt_uORF",]$melsim_cons,y6[y6$class=="300-500nt_uORF",]$len,
                     y6[y6$class==">500nt_uORF",]$melsim_cons,y6[y6$class==">500nt_uORF",]$len),nrow=2))


fisher.test(matrix(c(y6[y6$class=="0-100nt_uORF",]$gain_before_b3,y6[y6$class=="0-100nt_uORF",]$len,
                     y6[y6$class=="100-300nt_uORF",]$gain_before_b3,y6[y6$class=="100-300nt_uORF",]$len),nrow=2))
fisher.test(matrix(c(y6[y6$class=="0-100nt_uORF",]$gain_before_b3,y6[y6$class=="0-100nt_uORF",]$len,
                     y6[y6$class=="300-500nt_uORF",]$gain_before_b3,y6[y6$class=="300-500nt_uORF",]$len),nrow=2))
fisher.test(matrix(c(y6[y6$class=="0-100nt_uORF",]$gain_before_b3,y6[y6$class=="0-100nt_uORF",]$len,
                     y6[y6$class==">500nt_uORF",]$gain_before_b3,y6[y6$class==">500nt_uORF",]$len),nrow=2))

fisher.test(matrix(c(y6[y6$class=="100-300nt_uORF",]$gain_before_b3,y6[y6$class=="100-300nt_uORF",]$len,
                     y6[y6$class=="300-500nt_uORF",]$gain_before_b3,y6[y6$class=="300-500nt_uORF",]$len),nrow=2))
fisher.test(matrix(c(y6[y6$class=="100-300nt_uORF",]$gain_before_b3,y6[y6$class=="100-300nt_uORF",]$len,
                     y6[y6$class==">500nt_uORF",]$gain_before_b3,y6[y6$class==">500nt_uORF",]$len),nrow=2))
fisher.test(matrix(c(y6[y6$class=="300-500nt_uORF",]$gain_before_b3,y6[y6$class=="300-500nt_uORF",]$len,
                     y6[y6$class==">500nt_uORF",]$gain_before_b3,y6[y6$class==">500nt_uORF",]$len),nrow=2))





#2.compensation_ratio

#########2.1 compensation for specific regions#########
uORF_matrix_canonical=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/8.new_27_mfa/PacBioSim/position_cds/uORF_matrix_dm6_PacBioSim_27species_update_noCDS_postoCDS_canonical.csv",sep=',',header=T)
h=uORF_matrix_canonical[,c(1:4,7,11,31:37)]
h=h[h$dm6+h$PacBioSim+h$droYak3+h$droAna3>0,]

gene_num=matrix(0,ncol=6,nrow=6) %>% as.data.table
colnames(gene_num)=c("mel_loss","mel_gain","anc_loss","anc_gain","gain_loss","loss_gain")

for (i in c(1:5)){
  s=as.character(i*100)
  x=matrix(0,ncol=2,nrow=4) %>% as.data.table
  colnames(x)=c("uORF","gene")
  x$type=c("mel_loss","mel_gain","anc_loss","anc_gain")
  
  h2=h[h$pos_to_cds<=i*100,]
  #class4
  h2[h2$dm6==0 & h2$PacBioSim==1 & h2$droYak3==1,]->z 
  x[x$type=="mel_loss",]$uORF=length(z$ID)
  gene_num$"mel_loss"[i]<-x[x$type=="mel_loss",]$gene<-length(unique(z$geneID))
  class4_list=unique(z$geneID)
  class4_tr_list=unique(z$transcriptID)
  #class3
  h2[h2$dm6==1 & h2$PacBioSim==0 & h2$droYak3==0,]->z
  x[x$type=="mel_gain",]$uORF=length(z$ID)
  gene_num$"mel_gain"[i]<-x[x$type=="mel_gain",]$gene<-length(unique(z$geneID))
  class3_list=unique(z$geneID)
  class3_tr_list=unique(z$transcriptID)
  #class2
  h2[h2$dm6==0 & h2$PacBioSim==0 & h2$droYak3==1 & h2$droAna3==1,]->z
  x[x$type=="anc_loss",]$uORF=length(z$ID)
  gene_num$"anc_loss"[i]<-x[x$type=="anc_loss",]$gene<-length(unique(z$geneID))
  class2_list=unique(z$geneID)
  class2_tr_list=unique(z$transcriptID)
  #class1
  h2[h2$dm6==1 & h2$PacBioSim==1 & h2$droYak3==0 & h2$droAna3==0,]->z
  x[x$type=="anc_gain",]$uORF=length(z$ID)
  gene_num$"anc_gain"[i]<-x[x$type=="anc_gain",]$gene<-length(unique(z$geneID)) #713
  class1_list=unique(z$geneID)
  class1_tr_list=unique(z$transcriptID)
  
  fwrite(x,paste0("compensation_stat_",s,".txt"),sep='\t')
  #gene
  gain_loss=intersect(class1_list,class4_list)
  gene_num$"gain_loss"[i]=length(gain_loss)
  loss_gain=intersect(class2_list,class3_list)
  gene_num$"loss_gain"[i]=length(loss_gain)
  anc_gain=class1_list
  anc_loss=class2_list
  mel_gain=class3_list
  mel_loss=class4_list
  
  #transcript
  gain_loss_tr=intersect(class1_tr_list,class4_tr_list) #231
  loss_gain_tr=intersect(class2_tr_list,class3_tr_list) #82
  anc_gain_tr=class1_tr_list
  anc_loss_tr=class2_tr_list
  mel_gain_tr=class3_tr_list
  mel_loss_tr=class4_tr_list
  compensation_tr_list=list(anc_gain=anc_gain_tr,anc_loss=anc_loss_tr,mel_gain=mel_gain_tr,mel_loss=mel_loss_tr,gain_loss=gain_loss_tr,loss_gain=loss_gain_tr)
  filename=paste0("uORF_gain_loss_trlist_",s,".rds")
  saveRDS(compensation_tr_list, file = filename)
}

###all
h2=h
i=6
#class4
h2[h2$dm6==0 & h2$PacBioSim==1 & h2$droYak3==1,]->z 
x[x$type=="mel_loss",]$uORF=length(z$ID)
gene_num$"mel_loss"[i]<-x[x$type=="mel_loss",]$gene<-length(unique(z$geneID))
class4_list=unique(z$geneID)
class4_tr_list=unique(z$transcriptID)
#class3
h2[h2$dm6==1 & h2$PacBioSim==0 & h2$droYak3==0,]->z
x[x$type=="mel_gain",]$uORF=length(z$ID)
gene_num$"mel_gain"[i]<-x[x$type=="mel_gain",]$gene<-length(unique(z$geneID))
class3_list=unique(z$geneID)
class3_tr_list=unique(z$transcriptID)
#class2
h2[h2$dm6==0 & h2$PacBioSim==0 & h2$droYak3==1 & h2$droAna3==1,]->z
x[x$type=="anc_loss",]$uORF=length(z$ID)
gene_num$"anc_loss"[i]<-x[x$type=="anc_loss",]$gene<-length(unique(z$geneID))
class2_list=unique(z$geneID)
class2_tr_list=unique(z$transcriptID)
#class1
h2[h2$dm6==1 & h2$PacBioSim==1 & h2$droYak3==0 & h2$droAna3==0,]->z
x[x$type=="anc_gain",]$uORF=length(z$ID)
gene_num$"anc_gain"[i]<-x[x$type=="anc_gain",]$gene<-length(unique(z$geneID)) #713
class1_list=unique(z$geneID)
class1_tr_list=unique(z$transcriptID)

fwrite(x,"compensation_stat_all.txt",sep='\t')
#gene
gain_loss=intersect(class1_list,class4_list)
gene_num$"gain_loss"[i]=length(gain_loss)
loss_gain=intersect(class2_list,class3_list)
gene_num$"loss_gain"[i]=length(loss_gain)
anc_gain=class1_list
anc_loss=class2_list
mel_gain=class3_list
mel_loss=class4_list

#transcript
gain_loss_tr=intersect(class1_tr_list,class4_tr_list) #231
loss_gain_tr=intersect(class2_tr_list,class3_tr_list) #82
anc_gain_tr=class1_tr_list
anc_loss_tr=class2_tr_list
mel_gain_tr=class3_tr_list
mel_loss_tr=class4_tr_list
compensation_tr_list=list(anc_gain=anc_gain_tr,anc_loss=anc_loss_tr,mel_gain=mel_gain_tr,mel_loss=mel_loss_tr,gain_loss=gain_loss_tr,loss_gain=loss_gain_tr)
filename=paste0("uORF_gain_loss_trlist_all.rds")
saveRDS(compensation_tr_list, file = filename)

gene_num$region=c("0-100nt","0-200nt","0-300nt","0-400nt","0-500nt","all")
fwrite(gene_num,"compensation_gene_number_stat.txt",sep='\t')

#gene_num2=t(gene_num) %>% as.data.table
#gene_num2$type=c("mel_loss","mel_gain","anc_loss","anc_gain","gain_loss","loss_gain")
#colnames(gene_num2)=c("0-100nt","0-200nt","0-300nt","0-400nt","0-500nt","all","type")

gene_num$region=c("0-100nt","0-200nt","0-300nt","0-400nt","0-500nt","all")
gene_num$gain_loss2anc_gain=gene_num$gain_loss/gene_num$anc_gain
gene_num$loss_gain2anc_loss=gene_num$loss_gain/gene_num$anc_loss
gene_num$gain_loss2mel_loss=gene_num$gain_loss/gene_num$mel_loss
gene_num$loss_gain2mel_gain=gene_num$loss_gain/gene_num$mel_gain

gene_num$region=factor(gene_num$region,levels=c(gene_num$region))

ggplot(data = gene_num, mapping = aes(x = region, y = gain_loss2anc_gain,fill=region)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("gain_loss/anc_gain")+
  theme_classic() #gain_loss2anc_gain.pdf

fisher.test(matrix(c(gene_num[gene_num$region=="0-100nt",]$anc_gain,gene_num[gene_num$region=="0-100nt",]$gain_loss,
                     gene_num[gene_num$region=="0-300nt",]$anc_gain,gene_num[gene_num$region=="0-300nt",]$gain_loss),nrow=2))

ggplot(data = gene_num, mapping = aes(x = region, y = loss_gain2anc_loss,fill=region)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("loss_gain/anc_loss")+
  theme_classic() #loss_gain2anc_loss.pdf


ggplot(data = gene_num, mapping = aes(x = region, y = gene_num$gain_loss/231,fill=region)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("gain_loss/all_gain_loss")+
  theme_classic() #gain_loss2all_anc_gain.pdf

#normalized by conserved uATG
gene_num$conserved=0
h2=h
h2$distance_class="0-100nt"
h2[h2$pos_to_cds>100,]$distance_class="100-200nt"
h2[h2$pos_to_cds>200,]$distance_class="200-300nt"
h2[h2$pos_to_cds>300,]$distance_class="300-400nt"
h2[h2$pos_to_cds>400,]$distance_class="400-500nt"
h2[h2$pos_to_cds>500,]$distance_class=">500nt"
h2=h2[h2$PacBioSim==1 & h2$dm6==1,]
gene_num$conserved[1]=length(h2[h2$distance_class=="0-100nt",]$ID)
gene_num$conserved[2]=length(h2[h2$distance_class=="0-100nt" |h2$distance_class=="100-200nt" ,]$ID)
gene_num$conserved[3]=length(h2[h2$distance_class=="0-100nt"|h2$distance_class=="100-200nt" |h2$distance_class=="200-300nt" ,]$ID)
gene_num$conserved[4]=length(h2[h2$distance_class=="0-100nt"|h2$distance_class=="100-200nt" |h2$distance_class=="200-300nt" |h2$distance_class=="300-400nt",]$ID)
gene_num$conserved[5]=length(h2[h2$distance_class=="0-100nt"|h2$distance_class=="100-200nt" |h2$distance_class=="200-300nt" |h2$distance_class=="300-400nt"|h2$distance_class=="400-500nt",]$ID)
gene_num$conserved[6]=length(h2$ID)

ggplot(data = gene_num, mapping = aes(x = region, y = (gene_num$gain_loss+gene_num$loss_gain)/gene_num$conserved,fill=region)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("compensation/conserved")+
  theme_classic() #compensation2conserved.pdf

####2.2 compensation when only considering position of anc G/L ####
#uORF_gain_loss_trlist_all<-readRDS("uORF_gain_loss_trlist_all.rds",refhook = NULL)
#uORF_gain_loss_trlist_all$gain_loss

uORF_matrix_canonical=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/8.new_27_mfa/PacBioSim/position_cds/uORF_matrix_dm6_PacBioSim_27species_update_noCDS_postoCDS_canonical.csv",sep=',',header=T)
h=uORF_matrix_canonical[,c(1:4,7,11,31:37)]
h=h[h$dm6+h$PacBioSim+h$droYak3+h$droAna3>0,]

gene_num=matrix(0,ncol=6,nrow=6) %>% as.data.table
colnames(gene_num)=c("mel_loss","mel_gain","anc_loss","anc_gain","gain_loss","loss_gain")

for (i in c(1:5)){
  s=as.character(i*100)
  x=matrix(0,ncol=2,nrow=4) %>% as.data.table
  colnames(x)=c("uORF","gene")
  x$type=c("mel_loss","mel_gain","anc_loss","anc_gain")
  
  h2=h[h$pos_to_cds<=i*100 &h$pos_to_cds>=(i-1)*100,]
  #class4
  h[h$dm6==0 & h$PacBioSim==1 & h$droYak3==1,]->z 
  x[x$type=="mel_loss",]$uORF=length(z$ID)
  gene_num$"mel_loss"[i]<-x[x$type=="mel_loss",]$gene<-length(unique(z$geneID))
  class4_list=unique(z$geneID)
  class4_tr_list=unique(z$transcriptID)
  #class3
  h[h$dm6==1 & h$PacBioSim==0 & h$droYak3==0,]->z
  x[x$type=="mel_gain",]$uORF=length(z$ID)
  gene_num$"mel_gain"[i]<-x[x$type=="mel_gain",]$gene<-length(unique(z$geneID))
  class3_list=unique(z$geneID)
  class3_tr_list=unique(z$transcriptID)
  #class2
  h2[h2$dm6==0 & h2$PacBioSim==0 & h2$droYak3==1 & h2$droAna3==1,]->z
  x[x$type=="anc_loss",]$uORF=length(z$ID)
  gene_num$"anc_loss"[i]<-x[x$type=="anc_loss",]$gene<-length(unique(z$geneID))
  class2_list=unique(z$geneID)
  class2_tr_list=unique(z$transcriptID)
  #class1
  h2[h2$dm6==1 & h2$PacBioSim==1 & h2$droYak3==0 & h2$droAna3==0,]->z
  x[x$type=="anc_gain",]$uORF=length(z$ID)
  gene_num$"anc_gain"[i]<-x[x$type=="anc_gain",]$gene<-length(unique(z$geneID)) #713
  class1_list=unique(z$geneID)
  class1_tr_list=unique(z$transcriptID)
  
  #fwrite(x,paste0("compensation_stat_",s,".txt"),sep='\t')
  #gene
  gain_loss=intersect(class1_list,class4_list)
  gene_num$"gain_loss"[i]=length(gain_loss)
  loss_gain=intersect(class2_list,class3_list)
  gene_num$"loss_gain"[i]=length(loss_gain)
  anc_gain=class1_list
  anc_loss=class2_list
  mel_gain=class3_list
  mel_loss=class4_list
  
  #transcript
  gain_loss_tr=intersect(class1_tr_list,class4_tr_list) 
  loss_gain_tr=intersect(class2_tr_list,class3_tr_list) 
  anc_gain_tr=class1_tr_list
  anc_loss_tr=class2_tr_list
  mel_gain_tr=class3_tr_list
  mel_loss_tr=class4_tr_list
  compensation_tr_list=list(anc_gain=anc_gain_tr,anc_loss=anc_loss_tr,mel_gain=mel_gain_tr,mel_loss=mel_loss_tr,gain_loss=gain_loss_tr,loss_gain=loss_gain_tr)
  #filename=paste0("uORF_gain_loss_trlist_",s,".rds")
  #saveRDS(compensation_tr_list, file = filename)
}

###all
h2=h
i=6
#class4
h2[h2$dm6==0 & h2$PacBioSim==1 & h2$droYak3==1,]->z 
x[x$type=="mel_loss",]$uORF=length(z$ID)
gene_num$"mel_loss"[i]<-x[x$type=="mel_loss",]$gene<-length(unique(z$geneID))
class4_list=unique(z$geneID)
class4_tr_list=unique(z$transcriptID)
#class3
h2[h2$dm6==1 & h2$PacBioSim==0 & h2$droYak3==0,]->z
x[x$type=="mel_gain",]$uORF=length(z$ID)
gene_num$"mel_gain"[i]<-x[x$type=="mel_gain",]$gene<-length(unique(z$geneID))
class3_list=unique(z$geneID)
class3_tr_list=unique(z$transcriptID)
#class2
h2[h2$dm6==0 & h2$PacBioSim==0 & h2$droYak3==1 & h2$droAna3==1,]->z
x[x$type=="anc_loss",]$uORF=length(z$ID)
gene_num$"anc_loss"[i]<-x[x$type=="anc_loss",]$gene<-length(unique(z$geneID))
class2_list=unique(z$geneID)
class2_tr_list=unique(z$transcriptID)
#class1
h2[h2$dm6==1 & h2$PacBioSim==1 & h2$droYak3==0 & h2$droAna3==0,]->z
x[x$type=="anc_gain",]$uORF=length(z$ID)
gene_num$"anc_gain"[i]<-x[x$type=="anc_gain",]$gene<-length(unique(z$geneID)) #713
class1_list=unique(z$geneID)
class1_tr_list=unique(z$transcriptID)

#fwrite(x,"compensation_stat_all.txt",sep='\t')
#gene
gain_loss=intersect(class1_list,class4_list)
gene_num$"gain_loss"[i]=length(gain_loss)
loss_gain=intersect(class2_list,class3_list)
gene_num$"loss_gain"[i]=length(loss_gain)
anc_gain=class1_list
anc_loss=class2_list
mel_gain=class3_list
mel_loss=class4_list

#transcript
gain_loss_tr=intersect(class1_tr_list,class4_tr_list) #231
loss_gain_tr=intersect(class2_tr_list,class3_tr_list) #82
anc_gain_tr=class1_tr_list
anc_loss_tr=class2_tr_list
mel_gain_tr=class3_tr_list
mel_loss_tr=class4_tr_list
compensation_tr_list=list(anc_gain=anc_gain_tr,anc_loss=anc_loss_tr,mel_gain=mel_gain_tr,mel_loss=mel_loss_tr,gain_loss=gain_loss_tr,loss_gain=loss_gain_tr)


gene_num$region=c("0-100nt","100-200nt","200-300nt","300-400nt","400-500nt","all")
gene_num$gain_loss2anc_gain=gene_num$gain_loss/gene_num$anc_gain
gene_num$loss_gain2anc_loss=gene_num$loss_gain/gene_num$anc_loss
gene_num$gain_loss2mel_loss=gene_num$gain_loss/gene_num$mel_loss
gene_num$loss_gain2mel_gain=gene_num$loss_gain/gene_num$mel_gain

gene_num$region=factor(gene_num$region,levels=c(gene_num$region))

ggplot(data = gene_num, mapping = aes(x = region, y = gain_loss2anc_gain,fill=region)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("gain_loss/anc_gain")+
  xlab("region of anc gain")+
  theme_classic() #gain_loss2anc_gain_by_region_of_anc_gain.pdf

fisher.test(matrix(c(gene_num[gene_num$region=="0-100nt",]$anc_gain,gene_num[gene_num$region=="0-100nt",]$gain_loss,
                     gene_num[gene_num$region=="100-200nt",]$anc_gain,gene_num[gene_num$region=="100-200nt",]$gain_loss),nrow=2))

ggplot(data = gene_num, mapping = aes(x = region, y = loss_gain2anc_loss,fill=region)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("loss_gain/anc_loss")+
  xlab("region of anc loss")+
  theme_classic() #loss_gain2anc_loss_by_region_of_anc_loss.pdf
fisher.test(matrix(c(gene_num[gene_num$region=="0-100nt",]$anc_loss,gene_num[gene_num$region=="0-100nt",]$loss_gain,
                     gene_num[gene_num$region=="200-300nt",]$anc_loss,gene_num[gene_num$region=="200-300nt",]$loss_gain),nrow=2))
fisher.test(matrix(c(gene_num[gene_num$region=="300-400nt",]$anc_loss,gene_num[gene_num$region=="300-400nt",]$loss_gain,
                     gene_num[gene_num$region=="200-300nt",]$anc_loss,gene_num[gene_num$region=="200-300nt",]$loss_gain),nrow=2))




####2.3 compensation when only considering position of mel G/L ####
uORF_matrix_canonical=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/8.new_27_mfa/PacBioSim/position_cds/uORF_matrix_dm6_PacBioSim_27species_update_noCDS_postoCDS_canonical.csv",sep=',',header=T)
h=uORF_matrix_canonical[,c(1:4,7,11,31:37)]
h=h[h$dm6+h$PacBioSim+h$droYak3+h$droAna3>0,]

gene_num=matrix(0,ncol=6,nrow=6) %>% as.data.table
colnames(gene_num)=c("mel_loss","mel_gain","anc_loss","anc_gain","gain_loss","loss_gain")

for (i in c(1:5)){
  s=as.character(i*100)
  x=matrix(0,ncol=2,nrow=4) %>% as.data.table
  colnames(x)=c("uORF","gene")
  x$type=c("mel_loss","mel_gain","anc_loss","anc_gain")
  
  h2=h[h$pos_to_cds<=i*100 &h$pos_to_cds>=(i-1)*100,]
  #class4
  h2[h2$dm6==0 & h2$PacBioSim==1 & h2$droYak3==1,]->z 
  x[x$type=="mel_loss",]$uORF=length(z$ID)
  gene_num$"mel_loss"[i]<-x[x$type=="mel_loss",]$gene<-length(unique(z$geneID))
  class4_list=unique(z$geneID)
  class4_tr_list=unique(z$transcriptID)
  #class3
  h2[h2$dm6==1 & h2$PacBioSim==0 & h2$droYak3==0,]->z
  x[x$type=="mel_gain",]$uORF=length(z$ID)
  gene_num$"mel_gain"[i]<-x[x$type=="mel_gain",]$gene<-length(unique(z$geneID))
  class3_list=unique(z$geneID)
  class3_tr_list=unique(z$transcriptID)
  #class2
  h[h$dm6==0 & h$PacBioSim==0 & h$droYak3==1 & h$droAna3==1,]->z
  x[x$type=="anc_loss",]$uORF=length(z$ID)
  gene_num$"anc_loss"[i]<-x[x$type=="anc_loss",]$gene<-length(unique(z$geneID))
  class2_list=unique(z$geneID)
  class2_tr_list=unique(z$transcriptID)
  #class1
  h[h$dm6==1 & h$PacBioSim==1 & h$droYak3==0 & h$droAna3==0,]->z
  x[x$type=="anc_gain",]$uORF=length(z$ID)
  gene_num$"anc_gain"[i]<-x[x$type=="anc_gain",]$gene<-length(unique(z$geneID)) #713
  class1_list=unique(z$geneID)
  class1_tr_list=unique(z$transcriptID)
  
  #fwrite(x,paste0("compensation_stat_",s,".txt"),sep='\t')
  #gene
  gain_loss=intersect(class1_list,class4_list)
  gene_num$"gain_loss"[i]=length(gain_loss)
  loss_gain=intersect(class2_list,class3_list)
  gene_num$"loss_gain"[i]=length(loss_gain)
  anc_gain=class1_list
  anc_loss=class2_list
  mel_gain=class3_list
  mel_loss=class4_list
  
  #transcript
  gain_loss_tr=intersect(class1_tr_list,class4_tr_list) 
  loss_gain_tr=intersect(class2_tr_list,class3_tr_list) 
  anc_gain_tr=class1_tr_list
  anc_loss_tr=class2_tr_list
  mel_gain_tr=class3_tr_list
  mel_loss_tr=class4_tr_list
  compensation_tr_list=list(anc_gain=anc_gain_tr,anc_loss=anc_loss_tr,mel_gain=mel_gain_tr,mel_loss=mel_loss_tr,gain_loss=gain_loss_tr,loss_gain=loss_gain_tr)
  #filename=paste0("uORF_gain_loss_trlist_",s,".rds")
  #saveRDS(compensation_tr_list, file = filename)
}

###all
h2=h
i=6
#class4
h2[h2$dm6==0 & h2$PacBioSim==1 & h2$droYak3==1,]->z 
x[x$type=="mel_loss",]$uORF=length(z$ID)
gene_num$"mel_loss"[i]<-x[x$type=="mel_loss",]$gene<-length(unique(z$geneID))
class4_list=unique(z$geneID)
class4_tr_list=unique(z$transcriptID)
#class3
h2[h2$dm6==1 & h2$PacBioSim==0 & h2$droYak3==0,]->z
x[x$type=="mel_gain",]$uORF=length(z$ID)
gene_num$"mel_gain"[i]<-x[x$type=="mel_gain",]$gene<-length(unique(z$geneID))
class3_list=unique(z$geneID)
class3_tr_list=unique(z$transcriptID)
#class2
h2[h2$dm6==0 & h2$PacBioSim==0 & h2$droYak3==1 & h2$droAna3==1,]->z
x[x$type=="anc_loss",]$uORF=length(z$ID)
gene_num$"anc_loss"[i]<-x[x$type=="anc_loss",]$gene<-length(unique(z$geneID))
class2_list=unique(z$geneID)
class2_tr_list=unique(z$transcriptID)
#class1
h2[h2$dm6==1 & h2$PacBioSim==1 & h2$droYak3==0 & h2$droAna3==0,]->z
x[x$type=="anc_gain",]$uORF=length(z$ID)
gene_num$"anc_gain"[i]<-x[x$type=="anc_gain",]$gene<-length(unique(z$geneID)) #713
class1_list=unique(z$geneID)
class1_tr_list=unique(z$transcriptID)

#fwrite(x,"compensation_stat_all.txt",sep='\t')
#gene
gain_loss=intersect(class1_list,class4_list)
gene_num$"gain_loss"[i]=length(gain_loss)
loss_gain=intersect(class2_list,class3_list)
gene_num$"loss_gain"[i]=length(loss_gain)
anc_gain=class1_list
anc_loss=class2_list
mel_gain=class3_list
mel_loss=class4_list

#transcript
gain_loss_tr=intersect(class1_tr_list,class4_tr_list) #231
loss_gain_tr=intersect(class2_tr_list,class3_tr_list) #82
anc_gain_tr=class1_tr_list
anc_loss_tr=class2_tr_list
mel_gain_tr=class3_tr_list
mel_loss_tr=class4_tr_list
compensation_tr_list=list(anc_gain=anc_gain_tr,anc_loss=anc_loss_tr,mel_gain=mel_gain_tr,mel_loss=mel_loss_tr,gain_loss=gain_loss_tr,loss_gain=loss_gain_tr)


gene_num$region=c("0-100nt","100-200nt","200-300nt","300-400nt","400-500nt","all")
gene_num$gain_loss2anc_gain=gene_num$gain_loss/gene_num$anc_gain
gene_num$loss_gain2anc_loss=gene_num$loss_gain/gene_num$anc_loss
gene_num$gain_loss2mel_loss=gene_num$gain_loss/gene_num$mel_loss
gene_num$loss_gain2mel_gain=gene_num$loss_gain/gene_num$mel_gain

gene_num$region=factor(gene_num$region,levels=c(gene_num$region))

ggplot(data = gene_num, mapping = aes(x = region, y = gain_loss2mel_loss,fill=region)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("gain_loss/mel_loss")+
  xlab("region of mel loss")+
  theme_classic() #gain_loss2mel_loss_by_region_of_mel_loss.pdf

fisher.test(matrix(c(gene_num[gene_num$region=="0-100nt",]$mel_loss,gene_num[gene_num$region=="0-100nt",]$gain_loss,
                     gene_num[gene_num$region=="100-200nt",]$mel_loss,gene_num[gene_num$region=="100-200nt",]$gain_loss),nrow=2))

ggplot(data = gene_num, mapping = aes(x = region, y = loss_gain2mel_gain,fill=region)) + 
  geom_bar(color="black",stat = 'identity')+
  scale_fill_brewer(palette="BuPu")+
  ylab("loss_gain/mel_gain")+
  xlab("region of mel gain")+
  theme_classic() #loss_gain2anc_loss_by_region_of_mel_gain.pdf
fisher.test(matrix(c(gene_num[gene_num$region=="0-100nt",]$mel_gain,gene_num[gene_num$region=="0-100nt",]$loss_gain,
                     gene_num[gene_num$region=="200-300nt",]$mel_gain,gene_num[gene_num$region=="200-300nt",]$loss_gain),nrow=2))
fisher.test(matrix(c(gene_num[gene_num$region=="300-400nt",]$mel_gain,gene_num[gene_num$region=="300-400nt",]$loss_gain,
                     gene_num[gene_num$region=="200-300nt",]$mel_gain,gene_num[gene_num$region=="200-300nt",]$loss_gain),nrow=2))



####2.4 distribution of mel_gain/loss uATG that participate in compensation####  
#?????

library(data.table)
library(ggplot2)
library(dplyr)
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
hh<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS_postoCDS.csv")
hh_0<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv")
hh$ID=factor(hh$ID,levels=hh_0$V1)
hh=hh[order(hh$ID),]

hh[,c(2,3,30:35,37)]->hh
hh2 <- tidyr::extract(hh, col = 'ID', into = c('tr', 'tr_pos_Gap','tr_pos_NoGap'),
                      regex = '(FBtr\\d+)_(\\d+)[(](\\d+)[)]',remove=F) %>% as.data.table()
hh2[,geno_pos:=paste(seqname,gene_pos,sep = "_")]
hh2[,c("ID","dm6","tr","tr_pos_Gap","tr_pos_NoGap","geneID","geno_pos",
       "genesymbol","pos_to_cds")]->hh2
hh2[,POS:=1:.N]
hh2[order(geno_pos,pos_to_cds,POS)]->hh2 #total 435425
hh2[!duplicated(geno_pos)]->hh3 #total 224879

## add uATG gain/loss events with detailed information
jj<-fread("/data/uORF_matrix_triCas2_ATG_GLOOME_result/gainLossMP.1.PerPosPerBranch.txt")
names(jj)[1]<-"GL"
merge(jj,hh2,by="POS")->jj2
merge(jj,hh3,by="POS")->jj3

h2=jj2
h3=jj3

#canonical
canonical=fread("/data/dmel-all-r6.04.canonical.transcripts.txt",header=F,sep='\t')

#bins: 100

h2$distance_class=0 #"0-100nt"
h2[h2$pos_to_cds>100,]$distance_class=100 #"100-200nt"
h2[h2$pos_to_cds>200,]$distance_class=200 #"200-300nt"
h2[h2$pos_to_cds>300,]$distance_class=300 #"300-400nt"
h2[h2$pos_to_cds>400,]$distance_class=400 #"400-500nt"
h2[h2$pos_to_cds>500,]$distance_class=500 #"500-600nt"
h2[h2$pos_to_cds>600,]$distance_class=600 #"600-700nt"
h2[h2$pos_to_cds>700,]$distance_class=700 #">700nt"

h2[,.(num=.N),by=.(geneID,tr,GL,distance_class)]->GL_num_gene
dcast(GL_num_gene,geneID+tr+distance_class~GL,value.var ="num")->GL_num_gene2
GL_num_gene2[is.na(GL_num_gene2)]<-0
#GL_num_gene2=GL_num_gene2[GL_num_gene2$transcriptID %in% canonical$V3,]
GL_num_gene2=GL_num_gene2[GL_num_gene2$tr %in% canonical$V3,]
m=data.table(bin=c("0-100nt","100-200nt","200-300nt","300-400nt","400-500nt","500-600nt","600-700nt",">700nt"),cor=0)

m[m$bin=="0-100nt",]$cor=cor.test(GL_num_gene2[GL_num_gene2$distance_class==0,]$gain,GL_num_gene2[GL_num_gene2$distance_class==0,]$loss,method = "spearman")$estimate
m[m$bin=="100-200nt",]$cor=cor.test(GL_num_gene2[GL_num_gene2$distance_class==100,]$gain,GL_num_gene2[GL_num_gene2$distance_class==100,]$loss,method = "spearman")$estimate
m[m$bin=="200-300nt",]$cor=cor.test(GL_num_gene2[GL_num_gene2$distance_class==200,]$gain,GL_num_gene2[GL_num_gene2$distance_class==200,]$loss,method = "spearman")$estimate
m[m$bin=="300-400nt",]$cor=cor.test(GL_num_gene2[GL_num_gene2$distance_class==300,]$gain,GL_num_gene2[GL_num_gene2$distance_class==300,]$loss,method = "spearman")$estimate
m[m$bin=="400-500nt",]$cor=cor.test(GL_num_gene2[GL_num_gene2$distance_class==400,]$gain,GL_num_gene2[GL_num_gene2$distance_class==400,]$loss,method = "spearman")$estimate
m[m$bin=="500-600nt",]$cor=cor.test(GL_num_gene2[GL_num_gene2$distance_class==500,]$gain,GL_num_gene2[GL_num_gene2$distance_class==500,]$loss,method = "spearman")$estimate
m[m$bin=="600-700nt",]$cor=cor.test(GL_num_gene2[GL_num_gene2$distance_class==600,]$gain,GL_num_gene2[GL_num_gene2$distance_class==600,]$loss,method = "spearman")$estimate
m[m$bin==">700nt",]$cor=cor.test(GL_num_gene2[GL_num_gene2$distance_class==700,]$gain,GL_num_gene2[GL_num_gene2$distance_class==700,]$loss,method = "spearman")$estimate

m$bin=factor(m$bin,levels=c("0-100nt","100-200nt","200-300nt","300-400nt","400-500nt","500-600nt","600-700nt",">700nt"))
m$bin=factor(m$bin,levels=c(">700nt","600-700nt","500-600nt","400-500nt","300-400nt","200-300nt","100-200nt","0-100nt"))

#figS3

#figS5
pdf("/results/figS5A.pdf",width=4,height=4)
ggplot(data = m, mapping = aes(x = bin, y = cor,fill=bin)) + 
  geom_bar(color="black", stat = 'identity',fill="#4772B9",width=0.8)+
  #scale_fill_brewer(palette="BuPu")+
  ylab("rho")+
  coord_flip()+
  theme_classic()  #figS3A.pdf
dev.off()

pdf("/results/figS5A.pdf",width=4,height=4)
ggplot(GL_num_gene2[GL_num_gene2$distance_class==0,], aes(x=log10(gain+1), y=log10(loss+1))) +
  geom_point(size=3,shape=16,fill="black",alpha=0.05)+
  coord_cartesian(xlim = c(0, 2),ylim = c(0, 2))+ 
  labs(x = 'log10(gain occurrence)', y = 'log10(loss occurrence)', color = NULL)+
  geom_abline(intercept = 0, slope = 1,color="red",linetype="dashed", size=0.8)+
  theme_classic() #figS3B.pdf
cor.test(GL_num_gene2[GL_num_gene2$distance_class==0,]$gain,GL_num_gene2[GL_num_gene2$distance_class==0,]$loss,method = "spearman")
dev.off()
#	Spearman's rank correlation rho

#data:  GL_num_gene2[GL_num_gene2$distance_class == 0, ]$gain and GL_num_gene2[GL_num_gene2$distance_class == 0, ]$loss
#S = 109910000000, p-value < 0.00000000000000022
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.3852892 

######GO list######
#bins: 0-300, only canonical
h2=h2[h2$tr %in% canonical$V3,]
h2$distance_class=0
h2[h2$pos_to_cds>300,]$distance_class=300

h2[,.(num=.N),by=.(geneID,GL,distance_class)]->GL_num_gene
dcast(GL_num_gene,geneID+distance_class~GL,value.var ="num")->GL_num_gene2
#GL_num_gene2[!is.na(geneID)]->GL_num_gene2
GL_num_gene2[is.na(GL_num_gene2)]<-0
tmp=GL_num_gene2[GL_num_gene2$distance_class==0,]
tmp2=tmp[order(-(tmp$gain+tmp$loss)),]
fwrite(tmp2[1:500,1],"/results/top500gainloss_genelist_canonical_300nt.txt") #for metascape

tmp2=tmp[order(-(tmp$gain-tmp$loss)),]
fwrite(tmp2[1:500,1],"/results/top500netgain_genelist_canonical_300nt.txt")

#tmp2=tmp[order(-(tmp$gain)),]
#fwrite(tmp2[1:500,1],"/results/top500gain_genelist_canonical_300nt.txt")

#tmp2=tmp[order(-(tmp$loss)),]
#fwrite(tmp2[1:500,1],"/results/top500loss_genelist_canonical_300nt.txt")

#all (merged 5UTR)
h3[,.(num=.N),by=.(geneID,GL)]->GL_num_gene_merged
dcast(GL_num_gene_merged,geneID~GL,value.var ="num")->GL_num_gene_merged2
GL_num_gene_merged2[is.na(GL_num_gene_merged2)]<-0
tmp=GL_num_gene_merged2

tmp2=tmp[order(-(tmp$gain+tmp$loss)),]
fwrite(tmp2[1:500,1],"/results/top500gainloss_genelist_all.txt")

tmp2=tmp[order(-(tmp$gain-tmp$loss)),]
fwrite(tmp2[1:500,1],"/results/top500netgain_genelist_all.txt")

tmp2=tmp[order(-(tmp$gain)),]
fwrite(tmp2[1:500,1],"/results/top500gain_genelist_all.txt")

tmp2=tmp[order(-(tmp$loss)),]
fwrite(tmp2[1:500,1],"/results/top500loss_genelist_all.txt")



#normalized by length

GL_num_gene_merged2$gainandloss=GL_num_gene_merged2$gain+GL_num_gene_merged2$loss
GL_num_gene_merged2$netgain=GL_num_gene_merged2$gain-GL_num_gene_merged2$loss

#merged length
merged_utr_len<-fread("/data/dmel-all-r6.04-gene-merged-5utr.length",header = T)
colnames(merged_utr_len)=c("geneID", "mean", "median", "longest_isoform","merged")
GL_num_gene_merged3=merge(GL_num_gene_merged2,merged_utr_len,by="geneID")
GL_num_gene_merged3$gainandloss_by_mergedlen=GL_num_gene_merged3$gainandloss/GL_num_gene_merged3$merged
GL_num_gene_merged3$netgain_by_mergedlen=GL_num_gene_merged3$gainandloss/GL_num_gene_merged3$merged
GL_num_gene_merged3$gain_by_mergedlen=GL_num_gene_merged3$gain/GL_num_gene_merged3$merged
GL_num_gene_merged3$loss_by_mergedlen=GL_num_gene_merged3$loss/GL_num_gene_merged3$merged

tmp1=GL_num_gene_merged3[order(-(GL_num_gene_merged3$gainandloss_by_mergedlen)),]
tmp2=GL_num_gene_merged3[order(-(GL_num_gene_merged3$netgain_by_mergedlen)),]
#tmp3=GL_num_gene_merged3[order(-(GL_num_gene_merged3$gain_by_mergedlen)),]
#tmp4=GL_num_gene_merged3[order(-(GL_num_gene_merged3$loss_by_mergedlen)),]

fwrite(tmp1[1:500,1],"/results/top500gainloss_normalized_genelist_all.txt")
fwrite(tmp2[1:500,1],"/results/top500netgain_normalized_genelist_all.txt")
#fwrite(tmp3[1:500,1],"/results/top500gain_normalized_genelist_all.txt")
#fwrite(tmp4[1:500,1],"/results/top500loss_normalized_genelist_all.txt")



###bg
hh<-fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS_postoCDS.csv")
fwrite(unique(hh[,c("geneID")]),"/results/bg_gene_list")

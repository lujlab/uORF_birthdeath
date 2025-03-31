library(data.table)
library(ggplot2)
library(magrittr)
library(dplyr)
#library(lessR)
#setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/pop_var/uORF/2.vep_UTRannotator_new/plot")
#MAF
# MAF distribution
f1=fread("/data/1356_uORF_var_type_SNPonly_rmCDS_new.bed",header=F)
f1=f1[!duplicated(f1$V4),] #39101
f2=fread("/data/1356_uORF_var_type_INDELonly_rmCDS_new.bed",header=F)
f2=f2[!duplicated(f2$V4),] #17570
f1$type="SNP"
f2$type="indel"
f=rbind(f1,f2)
f=f[f$V14>0 &f$V14<1,]
length(unique(f$V5))#8321, not include all genescan not represent genes
#stats
f$MAF=f$V14
f[f$V14>0.5,]$MAF=1-f[f$V14>0.5,]$V14

table(f[f$type=="SNP",]$V1)
#2L   2R   3L   3R    4    X 
#7844 8105 7389 8334  233 7196

table(f[f$type=="SNP"&f$MAF>=0.01,]$V1)
#2L   2R   3L   3R    4    X 
#1351 1388 1285 1407   26 1006 

table(f[f$type=="SNP"&f$MAF>=0.05,]$V1)
#2L   2R   3L   3R    4    X 
#640 634 557 629   9 481 

table(f[f$type=="SNP"&f$MAF>=0.1,]$V1)
# 2L  2R  3L  3R   4   X 
#439 456 387 426   6 357 

table(f[f$type=="indel",]$V1)
#  2L   2R   3L   3R    4    X 
#3546 3618 3068 3725   99 3514

table(f[f$type=="indel"&f$MAF>=0.01,]$V1)
# 2L  2R  3L  3R   4   X 
#704 692 605 698  20 638
table(f[f$type=="indel"&f$MAF>=0.05,]$V1)
# 2L  2R  3L  3R   4   X 
#336 304 269 295   4 277
table(f[f$type=="indel"&f$MAF>=0.1,]$V1)
# 2L  2R  3L  3R   4   X 
#232 192 185 197   4 187


summary(f$V14)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0003687 0.0007396 0.0014771 0.0264716 0.0048255 0.9996300 
tmp=f

tmp[tmp$MAF>0.05,]$MAF<-0.051
tmp[tmp$V14>0.05,]$V14<-0.051
f_maf=tmp[,c("V4","V14","MAF")]
f_maf=tmp
length(unique(tmp$V4)) # 56671 variants
length(unique(f$V5)) #8321 genes

pdf("/results/fig6b_MAF_uORF_var.pdf",width=4,height=4)
hist(tmp$MAF,breaks = 50,col = "gray",xlim=c(0,0.051),
     ylim=c(0,30000),
     main = "",las=1,xlab = "Minor allele frequency of uORF-related variants",ylab="Variants Counts")
dev.off()


#hist(tmp$V14,breaks = 50,col = "gray",xlim=c(0,0.051),
#     ylim=c(0,30000),
#     main = "",las=1,xlab = "Alternative allele frequency of uORF-related variants",ylab="Variants Counts")
#4.18*4.42 fig6B_alternative_uORF_var.pdf

length(f[f$type=="SNP"&f$MAF<0.05,]$V1)/39101 #0.9245544
length(f[f$type=="indel"&f$MAF<0.05,]$V1)/17570 #0.9154809

# gene with variants number
g=fread("/data/var_gene_number",header=F)
tmp=g
tmp[tmp$V1>50,]$V1<-51

pdf("/results/figS25_number_of_uORFs_and_genes_Distribution.pdf",width=4,height=4)
hist(tmp$V1,breaks= seq(0,51,1),col = "gray",ylim=c(0,2500),main = "",las=1,xlab = "Number of uORF-related variants",ylab="Gene Counts")
dev.off()
#figS_number_of_uORFs_and_genes_Distribution.pdf 
length(tmp[tmp$V1==1,]$V1) #2218
length(tmp[tmp$V1!=1,]$V1) #6216

# effect distribution
vt=fread("/data/varid_type_sort",header=F)
colnames(vt)=c("var_id","type")

tmp=f[,c("V4","V14","MAF")] #wrong var type in f, use var type in vt
colnames(tmp)=c("var_id","AF","MAF")
x=merge(tmp,vt)
#PieChart(type, data = vt,main = NULL,values = "input",values_size=1,labels_cex = 1.2,values_digits=0)
#compound  uAUG_gained  uAUG_lost  uFrameShift  uSTOP_gained  uSTOP_lost       Total 
#Frequencies:       3562        18269      10819        11695          6364        5962       56671 
#Proportions:      0.063        0.322      0.191        0.206         0.112       0.105       1.000 


#type_maf1=data.table(type=c("uAUG_gained","uAUG_lost","uFrameShift","uSTOP_gained","uSTOP_lost","compound"),
#                     MAF1=c(0.8560416,0.8126964,0.8008355,0.8255814,0.8137469,0.8160179))

x$type=factor(x$type,levels=c("uAUG_gained","uSTOP_gained","compound","uSTOP_lost","uAUG_lost","uFrameShift"))

breaks <- c(0, 0.01, 0.02, 0.03, 0.04, 0.05, Inf)


#derived allele
anc=fread("/data/anc_stat_56671.txt",sep='\t')
anc=anc[,c("id","stat")]
vt=fread("/data/varid_type_sort",header=F)
colnames(vt)=c("id","type")

tmp=f[,c("V4","V14","type","MAF")] #wrong var type in f
colnames(tmp)=c("id","AF","type1","MAF")
x=merge(tmp,vt,by="id")
x2=merge(anc,x,by="id")
x2$type_new=x2$type
x2[(x2$stat=="der"|x2$stat=="der_len")&x2$type=="uAUG_gained",]$type_new="uAUG_lost"
x2[(x2$stat=="der"|x2$stat=="der_len")&x2$type=="uAUG_lost",]$type_new="uAUG_gained"
x2[(x2$stat=="der"|x2$stat=="der_len")&x2$type=="uSTOP_gained",]$type_new="uSTOP_lost"
x2[(x2$stat=="der"|x2$stat=="der_len")&x2$type=="uSTOP_lost",]$type_new="uSTOP_gained"
table(x2$type_new)
#    compound  uAUG_gained    uAUG_lost  uFrameShift uSTOP_gained   uSTOP_lost 
#3562        18459        10629        11695         6383         5943 
table(x2$type)
#length(x2[x2$type_new=="uAUG_gained"&x2$MAF<0.01,]$id)/56671

#fwrite(x2,"/results/varid_type_sort_MAF.txt",sep='\t')


x2$DAF=x2$AF
x2[x2$stat=="der",]$DAF=1-x2[x2$stat=="der",]$AF


breaks <- c(0, 0.01, 0.02, 0.03, 0.04, 0.05, Inf)

# 将数据分配到这些 bin 中,DAF
x3 <- x2 %>%
  mutate(bin = cut(DAF, breaks = breaks, right = FALSE, include.lowest = TRUE, labels = c("0-0.01", "0.01-0.02", "0.02-0.03", "0.03-0.04", "0.04-0.05", ">0.05")))
#ggplot(x3, aes(x = bin, fill = type_new)) +
#  geom_histogram(stat = "count", position = "dodge") +
#  labs(x = "Value Bin", y = "Count", title = "Grouped Histogram with Custom Bins") +
#  scale_fill_brewer(palette = "Set1") +
#  theme_minimal()

x3_summary <- x3 %>%
  group_by(type_new, bin) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(type_new) %>%
  mutate(proportion = count / sum(count))


#add neutral control, DAF
si=fread("/data/all_shortintron.DAF",header=F,sep='\t')
colnames(si)=c("chr","pos","ref","alt","DAF")
si2 <- si %>%
  mutate(bin = cut(DAF, breaks = breaks, right = FALSE, include.lowest = TRUE, labels = c("0-0.01", "0.01-0.02", "0.02-0.03", "0.03-0.04", "0.04-0.05", ">0.05")))
si2$type_new="neutral"
si2_summary <- si2 %>%
  group_by(type_new, bin) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(type_new) %>%
  mutate(proportion = count / sum(count))
x4_summary=rbind(x3_summary,si2_summary)
x4_summary$type_new=factor(x4_summary$type_new,levels=c("uAUG_gained","compound","uSTOP_gained","uAUG_lost","uSTOP_lost","uFrameShift","neutral"))

p=ggplot(x4_summary, aes(x = bin, y = proportion, fill = type_new)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Derived frequency", y = "Proportion") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() #fig6C_6class_derived_frequency
pdf("/results/fig6c_6class_derived_frequency.pdf",width=6,height=6)
print(p)
dev.off()

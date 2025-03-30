library(data.table)
library(ggplot2)
library(gridExtra)
library(cowplot)
coor=fread("/data/rpkm/AbdB-J/Abd-B.1nt.bed",header=F)
####mRNA
#reads.cnt <- fread('paste /data/rpkm/AbdB-J/mrna/*_coverage', header = F)
f1=fread('/data/rpkm/AbdB-J/mrna/dmel_em_0_2h_coverage', header = F)
f2=fread('/data/rpkm/AbdB-J/mrna/dmel_em_2_6h_coverage', header = F)
f3=fread('/data/rpkm/AbdB-J/mrna/dmel_em_6_12h_coverage', header = F)
f4=fread('/data/rpkm/AbdB-J/mrna/dmel_em_12_24h_coverage', header = F)
reads.cnt=cbind(f1,f2,f3,f4)
x <- seq(7, ncol(reads.cnt), 7)
pos=reads.cnt[, .SD, .SDcols = c(1:3)]
reads.cnt=reads.cnt[, .SD, .SDcols = c(x)]

colnames(pos) <- c('chr',"pos0","pos")

#files <- fread('ls /data/rpkm/AbdB-J/mrna/*_coverage', header = F)
colnames(reads.cnt) <- c("dmel_em_0_2h_coverage",
                         "dmel_em_2_6h_coverage",
                         "dmel_em_6_12h_coverage",
                         "dmel_em_12_24h_coverage")

reads.cnt[reads.cnt=="."]=0
reads.cnt[, names(reads.cnt) := lapply(.SD, as.numeric)]
mrna_reads=cbind(pos,reads.cnt)

mrna_reads=mrna_reads[,c('chr',"pos0","pos",
                         "dmel_em_0_2h_coverage",
                         "dmel_em_2_6h_coverage",
                         "dmel_em_6_12h_coverage",
                         "dmel_em_12_24h_coverage")]

mrna_reads$relative_pos=(45025:1)
p_em02=ggplot(data=mrna_reads,aes(x=relative_pos,y=dmel_em_0_2h_coverage))+
  geom_bar(stat = 'identity',fill="blue")+
  theme_classic()
p_em26=ggplot(data=mrna_reads,aes(x=relative_pos,y=dmel_em_2_6h_coverage))+
  geom_bar(stat = 'identity',fill="blue")+
  theme_classic()
p_em612=ggplot(data=mrna_reads,aes(x=relative_pos,y=dmel_em_6_12h_coverage))+
  geom_bar(stat = 'identity',fill="blue")+
  theme_classic()
p_em1214=ggplot(data=mrna_reads,aes(x=relative_pos,y=dmel_em_12_24h_coverage))+
  geom_bar(stat = 'identity',fill="blue")+
  theme_classic()

#grid.arrange(p_em02,p_em26,p_em612,p_em1214,ncol=1)
pdf("/results/figS10_AbdBJ_Abd-b_mrna_cov_96block.pdf",width=8,height=8)
plot_grid(p_em02,p_em26,p_em612,p_em1214,ncol=1,align="v")
dev.off()
 #Abd-b_mrna_cov_96block,8*20


mrna_reads_96exon=mrna_reads[mrna_reads$pos>=16963159 & mrna_reads$pos<=16968574]
p_em02=ggplot(data=mrna_reads_96exon,aes(x=relative_pos,y=dmel_em_0_2h_coverage))+
  geom_bar(stat = 'identity',fill="blue")+
  theme_classic()
p_em26=ggplot(data=mrna_reads_96exon,aes(x=relative_pos,y=dmel_em_2_6h_coverage))+
  geom_bar(stat = 'identity',fill="blue")+
  theme_classic()
p_em612=ggplot(data=mrna_reads_96exon,aes(x=relative_pos,y=dmel_em_6_12h_coverage))+
  geom_bar(stat = 'identity',fill="blue")+
  theme_classic()
p_em1224=ggplot(data=mrna_reads_96exon,aes(x=relative_pos,y=dmel_em_12_24h_coverage))+
  geom_bar(stat = 'identity',fill="blue")+
  theme_classic()

#plot_grid(p_em02,p_em26,p_em612,p_em1224,ncol=1,align="v") #8*15


#####Psite#####
f1=fread('/data/rpkm/AbdB-J/ribo/dmel_em_0_2h_ribo_coverage', header = F)
f2=fread('/data/rpkm/AbdB-J/ribo/dmel_em_2_6h_ribo_coverage', header = F)
f3=fread('/data/rpkm/AbdB-J/ribo/dmel_em_6_12h_ribo_coverage', header = F)
f4=fread('/data/rpkm/AbdB-J/ribo/dmel_em_12_24h_ribo_coverage', header = F)
reads.cnt=cbind(f1,f2,f3,f4)

#reads.cnt <- fread('paste /data/rpkm/AbdB-J/ribo/*_coverage', header = F)
x <- seq(7, ncol(reads.cnt), 7)

pos=reads.cnt[, .SD, .SDcols = c(1:3)]
reads.cnt=reads.cnt[, .SD, .SDcols = c(x)]

colnames(pos) <- c('chr',"pos0","pos")


colnames(reads.cnt) <- c("dmel_em_0_2h_ribo_coverage",
                         "dmel_em_2_6h_ribo_coverage",
                         "dmel_em_6_12h_ribo_coverage",
                         "dmel_em_12_24h_ribo_coverage")

reads.cnt[reads.cnt=="."]=0
reads.cnt[, names(reads.cnt) := lapply(.SD, as.numeric)]
Psite_reads=cbind(pos,reads.cnt)
sample0<-c("dmel_em_0_2h","dmel_em_2_6h","dmel_em_6_12h","dmel_em_12_24h")
Psite_reads=Psite_reads[,c('chr',"pos0","pos",
                         "dmel_em_0_2h_ribo_coverage",
                         "dmel_em_2_6h_ribo_coverage",
                         "dmel_em_6_12h_ribo_coverage",
                         "dmel_em_12_24h_ribo_coverage")]

Psite_reads$relative_pos=(45025:1)
#p_em02_P=ggplot(data=Psite_reads,aes(x=relative_pos,y=dmel_em_0_2h_ribo_coverage))+
#  geom_bar(stat = 'identity',fill="blue")+
#  theme_classic()
#p_em26_P=ggplot(data=Psite_reads,aes(x=relative_pos,y=dmel_em_2_6h_ribo_coverage))+
#  geom_bar(stat = 'identity',fill="blue")+
#  theme_classic()
#p_em612_P=ggplot(data=Psite_reads,aes(x=relative_pos,y=dmel_em_6_12h_ribo_coverage))+
#  geom_bar(stat = 'identity',fill="blue")+
#  theme_classic()
#p_em1214_P=ggplot(data=Psite_reads,aes(x=relative_pos,y=dmel_em_12_24h_ribo_coverage))+
#  geom_bar(stat = 'identity',fill="blue")+
#  theme_classic()


Psite_reads_96exon=Psite_reads[Psite_reads$pos>=16963159 & Psite_reads$pos<=16968574]
p_em02_P=ggplot(data=Psite_reads_96exon,aes(x=relative_pos,y=dmel_em_0_2h_ribo_coverage))+
  geom_bar(stat = 'identity',fill="red")+
  theme_classic()
p_em26_P=ggplot(data=Psite_reads_96exon,aes(x=relative_pos,y=dmel_em_2_6h_ribo_coverage))+
  geom_bar(stat = 'identity',fill="red")+
  theme_classic()
p_em612_P=ggplot(data=Psite_reads_96exon,aes(x=relative_pos,y=dmel_em_6_12h_ribo_coverage))+
  geom_bar(stat = 'identity',fill="red")+
  theme_classic()
p_em1224_P=ggplot(data=Psite_reads_96exon,aes(x=relative_pos,y=dmel_em_12_24h_ribo_coverage))+
  geom_bar(stat = 'identity',fill="red")+
  theme_classic()

pdf("/results/figS10_AbdBJ_Psite_mRNA_4stage_96exon.pdf",width=8,height=16)
plot_grid(p_em02_P,p_em02,p_em26_P,p_em26,p_em612_P,p_em612,p_em1224_P,p_em1224,ncol=1,align="v") #8*15,Psite_mRNA_4stage_96exon
dev.off()


####TE
#uORF_matrix
uORF_matrix <- fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",sep=',',header=T)
colnames(uORF_matrix)[1]="uorf_align"
uORF_matrix=uORF_matrix[uORF_matrix$dm6>0 |uORF_matrix$PacBioSim>0 , ]
uORF_matrix=uORF_matrix[,c("uorf_align","dm6","PacBioSim","droYak3","transcriptID","position_withoutgap",
                           "seqname","geneID","gene_pos","genesymbol")]
uORF_matrix[, c("tmp1", "tmp2") := tstrsplit(uorf_align, "(", fixed = TRUE)]
uORF_matrix[, c("tmp1", "tmp2") := tstrsplit(tmp1, "_", fixed = TRUE)]
uORF_matrix$tmp2=as.numeric(uORF_matrix$tmp2)+1
uORF_matrix$uorf_id=paste(uORF_matrix$tmp1,uORF_matrix$tmp2,sep="_")
uORF_matrix=uORF_matrix[,c("uorf_align","uorf_id","dm6","PacBioSim","droYak3","transcriptID","position_withoutgap",
                           "seqname","geneID","gene_pos","genesymbol")]
# uORF TE,with cds overlapped
uORF_TE <- fread("/data/rpkm/dmel_uORF_mRNA_ribo_merge.txt",sep='\t',header=T)
uORF_matrix_TE=merge(uORF_matrix,uORF_TE,by="uorf_id") #rm overlap with cds

uORF_matrix_TE_mel=uORF_matrix_TE[uORF_matrix_TE$dm6==1,]
#unique(uORF_matrix_TE[uORF_matrix_TE$dm6==1,c("seqname","gene_pos")])

#max rpkm
#uORF_matrix_TE_mel$max_rpkm=0
#colnames(uORF_matrix_TE_mel)[c(78:97,109)] #rpkm
#colnames(uORF_matrix_TE_mel)[c(142:161,173)] #TE
uORF_matrix_TE_mel[is.na(uORF_matrix_TE_mel)]=0

uORF_matrix_TE_mel[, max_rpkm := apply(.SD, 1, max), .SDcols = c(78:97,109)]
uORF_matrix_TE_mel[, max_TE := apply(.SD, 1, max), .SDcols = c(142:161,173)]


uORF_matrix_TE_mel[uORF_matrix_TE_mel$geneID=="FBgn0000015",]
x=uORF_matrix_TE_mel[uORF_matrix_TE_mel$transcriptID=="FBtr0415465",]
length(uORF_matrix_TE_mel[uORF_matrix_TE_mel$transcriptID=="FBtr0415465",]$uorf_id)
x=x[order(x$position_withoutgap),]
#x$ymin=0:95
x$ymin=c(0:9,0:9,0:9,0:9,0:9,0:9,0:9,0:9,0:9,0:5)
#x$ymax=1:96
x$ymax=c(1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:6)
x$color="non"
x[x$max_TE>0.1,]$color="TE0.1"
nrow(x[x$max_TE>0.1,])
x[x$max_TE>0.5,]$color="TE0.5"


pdf("/results/figS10_AbdBJ_uORF_TE.pdf",width=8,height=4)
ggplot()+
  geom_rect(aes(xmin=x$position_withoutgap,xmax=x$dmel_uORF_len+x$position_withoutgap-1,ymin=x$ymin,ymax=x$ymax,fill=x$color))+
  geom_rect(aes(xmin=1,xmax=5649,ymin=12,ymax=13,fill="grey"))+
  scale_fill_manual(values=c("#CECFD1","#00BFC4", "#C77CFF","#F8766D"))+
  theme(axis.line=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank())
dev.off()




x=uORF_matrix_TE_mel[uORF_matrix_TE_mel$transcriptID=="FBtr0415463",]
length(uORF_matrix_TE_mel[uORF_matrix_TE_mel$transcriptID=="FBtr0415463",]$uorf_id)
x=x[order(x$position_withoutgap),]
#x$ymin=0:95
x$ymin=c(0:3)
#x$ymax=1:96
x$ymax=c(1:4)
x$color="non"
x[x$max_TE>0.1,]$color="TE0.1"
nrow(x[x$max_TE>0.1,])
x[x$max_TE>0.5,]$color="TE0.5"

pdf("/results/figS11a.pdf",width=8,height=2)
ggplot()+
  geom_rect(aes(xmin=x$position_withoutgap,xmax=x$dmel_uORF_len+x$position_withoutgap-1,ymin=x$ymin,ymax=x$ymax,fill=x$color))+
  geom_rect(aes(xmin=1,xmax=1211,ymin=5,ymax=6,fill="grey"))+
  scale_fill_manual(values=c("#CECFD1","#52A4DF"))+
  theme(legend.position = "none",axis.line=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank())
dev.off()


x=uORF_matrix_TE_mel[uORF_matrix_TE_mel$transcriptID=="FBtr0346764",]
length(uORF_matrix_TE_mel[uORF_matrix_TE_mel$transcriptID=="FBtr0346764",]$uorf_id)
x=x[order(x$position_withoutgap),]
#x$ymin=0:95
x$ymin=c(0:5)
#x$ymax=1:96
x$ymax=c(1:6)
x$color="non"
x[x$max_TE>0.1,]$color="TE0.1"
nrow(x[x$max_TE>0.1,])
x[x$max_TE>0.5,]$color="TE0.5"

pdf("/results/figS11b.pdf",width=8,height=2)
ggplot()+
  geom_rect(aes(xmin=x$position_withoutgap,xmax=x$dmel_uORF_len+x$position_withoutgap-1,ymin=x$ymin,ymax=x$ymax,fill=x$color))+
  geom_rect(aes(xmin=1,xmax=725,ymin=7,ymax=8,fill="grey"))+
  scale_fill_manual(values=c("#CECFD1","#52A4DF","#52A4DF"))+
  theme(legend.position = "none",axis.line=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank())
dev.off()

x=uORF_matrix_TE_mel[uORF_matrix_TE_mel$transcriptID=="FBtr0415465",]
length(uORF_matrix_TE_mel[uORF_matrix_TE_mel$transcriptID=="FBtr0415465",]$uorf_id)
x=x[order(x$position_withoutgap),]
#x$ymin=0:95
x$ymin=c(0:9,0:9,0:9,0:9,0:9,0:9,0:9,0:9,0:9,0:5)
#x$ymax=1:96
x$ymax=c(1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:6)
x$color="non"
x[x$max_TE>0.1,]$color="TE0.1"
nrow(x[x$max_TE>0.1,])
x[x$max_TE>0.5,]$color="TE0.5"
pdf("/results/figS11c.pdf",width=8,height=2)
ggplot()+
  geom_rect(aes(xmin=x$position_withoutgap,xmax=x$dmel_uORF_len+x$position_withoutgap-1,ymin=x$ymin,ymax=x$ymax,fill=x$color))+
  geom_rect(aes(xmin=1,xmax=5649,ymin=12,ymax=13,fill="grey"))+
  scale_fill_manual(values=c("#CECFD1","#52A4DF","#52A4DF","#52A4DF"))+
  theme(legend.position = "none",axis.line=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
        panel.background = element_blank()) #8*2
dev.off()

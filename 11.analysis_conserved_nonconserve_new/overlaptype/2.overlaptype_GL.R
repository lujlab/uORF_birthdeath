library(data.table)
library(ggplot2)
library(RColorBrewer)

#class
pos_class=fread("/data/rpkm/overlaptype/pos_uORF_class.txt",header=T)
pos_class$class="complex"
pos_class[nonoverlap==1,]$class="class1"
pos_class[inframe_overlapCDS==1 & outframe_overlapCDS==0 & nonoverlap==0 & inframe_overlapuORF==0 & outframe_overlapuORF==0,]$class="class2"
pos_class[inframe_overlapCDS==0 & outframe_overlapCDS==1 & nonoverlap==0 & inframe_overlapuORF==0 & outframe_overlapuORF==0,]$class="class3"
pos_class[inframe_overlapCDS==0 & outframe_overlapCDS==0 & nonoverlap==0 & inframe_overlapuORF==1 & outframe_overlapuORF==0,]$class="class4"
pos_class[inframe_overlapCDS==0 & outframe_overlapCDS==0 & nonoverlap==0 & inframe_overlapuORF==0 & outframe_overlapuORF==1,]$class="class5"
uATG_class=pos_class

#uORF matrix
uORF_matrix=fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",sep=',',header=T)
#remove duplicated
uORF_matrix_uniq=uORF_matrix[!duplicated(uORF_matrix, by = c("seqname", "gene_pos"))]
uORF_matrix_uniq$geno_pos=paste0(uORF_matrix_uniq$seqname,"_",uORF_matrix_uniq$gene_pos)
#only in mel
uORF_matrix_uniq_mel=uORF_matrix_uniq[uORF_matrix_uniq$dm6==1,]

#merge
uORF_matrix_uniq_mel_class=merge(uORF_matrix_uniq_mel,uATG_class,by="geno_pos")

#stat GL number for each class
stat_GLnum<-function(table_uORF_matrix){
  mel_sim_conserved=nrow(table_uORF_matrix[table_uORF_matrix$dm6==1&table_uORF_matrix$PacBioSim==1,])
  lost_in_sim=nrow(table_uORF_matrix[table_uORF_matrix$dm6==1 & table_uORF_matrix$PacBioSim==0 & table_uORF_matrix$droYak3==1,])
  gain_in_mel=nrow(table_uORF_matrix[table_uORF_matrix$dm6==1 & table_uORF_matrix$PacBioSim==0 & table_uORF_matrix$droYak3==0,])
  anc_gain=nrow(table_uORF_matrix[table_uORF_matrix$dm6==1 & table_uORF_matrix$PacBioSim==1 & table_uORF_matrix$droYak3==0 & table_uORF_matrix$droAna3==0,])
  sp4_1111=nrow(table_uORF_matrix[table_uORF_matrix$dm6==1 & table_uORF_matrix$PacBioSim==1 & table_uORF_matrix$droYak3==1 & table_uORF_matrix$droAna3==1,])
  sp4_1011=nrow(table_uORF_matrix[table_uORF_matrix$dm6==1 & table_uORF_matrix$PacBioSim==1 & table_uORF_matrix$droYak3==0 & table_uORF_matrix$droAna3==1,])
  sp4_0111=nrow(table_uORF_matrix[table_uORF_matrix$dm6==1 & table_uORF_matrix$PacBioSim==1 & table_uORF_matrix$droYak3==1 & table_uORF_matrix$droAna3==0,])
  out=data.table(sp4_1111=sp4_1111,sp4_1011=sp4_1011,sp4_0111=sp4_0111,mel_sim_conserved=mel_sim_conserved,
                 lost_in_sim=lost_in_sim,gain_in_mel=gain_in_mel,anc_gain=anc_gain)
  return(out)
}

uORF_matrix_uniq_mel_class1_num=stat_GLnum(uORF_matrix_uniq_mel_class[class=="class1",])
uORF_matrix_uniq_mel_class1_num$class="class1"
uORF_matrix_uniq_mel_class2_num=stat_GLnum(uORF_matrix_uniq_mel_class[class=="class2",])
uORF_matrix_uniq_mel_class2_num$class="class2"
uORF_matrix_uniq_mel_class3_num=stat_GLnum(uORF_matrix_uniq_mel_class[class=="class3",])
uORF_matrix_uniq_mel_class3_num$class="class3"
uORF_matrix_uniq_mel_class4_num=stat_GLnum(uORF_matrix_uniq_mel_class[class=="class4",])
uORF_matrix_uniq_mel_class4_num$class="class4"
uORF_matrix_uniq_mel_class5_num=stat_GLnum(uORF_matrix_uniq_mel_class[class=="class5",])
uORF_matrix_uniq_mel_class5_num$class="class5"
uORF_matrix_uniq_mel_complex_num=stat_GLnum(uORF_matrix_uniq_mel_class[class=="complex",])
uORF_matrix_uniq_mel_complex_num$class="complex"
uORF_matrix_uniq_mel_total_num=stat_GLnum(uORF_matrix_uniq_mel_class)
uORF_matrix_uniq_mel_total_num$class="total"

uORF_class_GL=rbind(uORF_matrix_uniq_mel_class1_num,uORF_matrix_uniq_mel_class2_num,
                    uORF_matrix_uniq_mel_class3_num,uORF_matrix_uniq_mel_class4_num,
                    uORF_matrix_uniq_mel_class5_num,uORF_matrix_uniq_mel_complex_num,
                    uORF_matrix_uniq_mel_total_num)

fwrite(uORF_class_GL,"/results/uORF_class_GL_mel_number.txt",sep='\t')


uORF_class_GL$total=uORF_class_GL$mel_sim_conserved+uORF_class_GL$lost_in_sim+uORF_class_GL$gain_in_mel
uORF_class_GL_p=uORF_class_GL
uORF_class_GL_oddsratio=uORF_class_GL
uORF_class_GL_ratio=uORF_class_GL
#fisher test, gain_in_mel in class1 / class1 vs. gain_in_mel in total / total = 30 / 227 vs. 4515 / 36565

#                   class0        total-class0                 total
#gain_in_mel        30          4515-3038                    4515
#total-gain_in_mel  227-3038    36565-227-(4515-30)      36565-4515
#total              21498         36565-21498                  36565
#a=30
#b=227
#c=4515
#d=36565
#fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))
#odds ratio>1  ->  more gain_in_mel in class0
for (i in c(1:6)){
  a=uORF_class_GL$sp4_1111[i]
  b=uORF_class_GL$total[i]
  c=uORF_class_GL$sp4_1111[7]
  d=uORF_class_GL$total[7]
  uORF_class_GL_p$sp4_1111[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_class_GL_oddsratio$sp4_1111[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_class_GL_ratio$sp4_1111[i]=a/b/(c/d)
  
  a=uORF_class_GL$sp4_1011[i]
  b=uORF_class_GL$total[i]
  c=uORF_class_GL$sp4_1011[7]
  d=uORF_class_GL$total[7]
  uORF_class_GL_p$sp4_1011[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_class_GL_oddsratio$sp4_1011[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_class_GL_ratio$sp4_1011[i]=a/b/(c/d)
  
  a=uORF_class_GL$sp4_0111[i]
  b=uORF_class_GL$total[i]
  c=uORF_class_GL$sp4_0111[7]
  d=uORF_class_GL$total[7]
  uORF_class_GL_p$sp4_0111[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_class_GL_oddsratio$sp4_0111[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_class_GL_ratio$sp4_0111[i]=a/b/(c/d)
  
  a=uORF_class_GL$mel_sim_conserved[i]
  b=uORF_class_GL$total[i]
  c=uORF_class_GL$mel_sim_conserved[7]
  d=uORF_class_GL$total[7]
  uORF_class_GL_p$mel_sim_conserved[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_class_GL_oddsratio$mel_sim_conserved[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_class_GL_ratio$mel_sim_conserved[i]=a/b/(c/d)
  
  a=uORF_class_GL$lost_in_sim[i]
  b=uORF_class_GL$total[i]
  c=uORF_class_GL$lost_in_sim[7]
  d=uORF_class_GL$total[7]
  uORF_class_GL_p$lost_in_sim[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_class_GL_oddsratio$lost_in_sim[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_class_GL_ratio$lost_in_sim[i]=a/b/(c/d)
  
  a=uORF_class_GL$gain_in_mel[i]
  b=uORF_class_GL$total[i]
  c=uORF_class_GL$gain_in_mel[7]
  d=uORF_class_GL$total[7]
  uORF_class_GL_p$gain_in_mel[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_class_GL_oddsratio$gain_in_mel[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_class_GL_ratio$gain_in_mel[i]=a/b/(c/d)
  
  a=uORF_class_GL$anc_gain[i]
  b=uORF_class_GL$total[i]
  c=uORF_class_GL$anc_gain[7]
  d=uORF_class_GL$total[7]
  uORF_class_GL_p$anc_gain[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_class_GL_oddsratio$anc_gain[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_class_GL_ratio$anc_gain[i]=a/b/(c/d)
}
fwrite(uORF_class_GL_p[1:6,1:8],"/results/uORF_class_GL_mel_number_p.txt",sep='\t')
fwrite(uORF_class_GL_oddsratio[1:6,1:8],"/results/uORF_class_GL_mel_number_oddsratio.txt",sep='\t')
#fwrite(uORF_class_GL_ratio[1:6,1:8],"uORF_class_GL_mel_number_enrichment.txt",sep='\t')

#heatmap
uORF_class_GL_p=fread("/results/uORF_class_GL_mel_number_p.txt")
uORF_class_GL_p2=melt(uORF_class_GL_p,id.vars="class",measure.vars=c("sp4_1111","sp4_1011","sp4_0111","mel_sim_conserved","lost_in_sim","gain_in_mel","anc_gain"))
colnames(uORF_class_GL_p2)=c("class","type","p")

uORF_class_GL_oddsratio=fread("/results/uORF_class_GL_mel_number_oddsratio.txt")
uORF_class_GL_oddsratio2=melt(uORF_class_GL_oddsratio,id.vars="class",measure.vars=c("sp4_1111","sp4_1011","sp4_0111","mel_sim_conserved","lost_in_sim","gain_in_mel","anc_gain"))
colnames(uORF_class_GL_oddsratio2)=c("class","type","oddsratio")

uORF_class_GL_oddsratio_p=merge(uORF_class_GL_oddsratio2,uORF_class_GL_p2,by=c("class","type"))
uORF_class_GL_oddsratio_p$p_label=""
uORF_class_GL_oddsratio_p[p<0.05,]$p_label="*"
uORF_class_GL_oddsratio_p[p<0.01,]$p_label="**"
uORF_class_GL_oddsratio_p[p<0.001,]$p_label="***"
uORF_class_GL_oddsratio_p$log2oddsratio=log2(uORF_class_GL_oddsratio_p$oddsratio)

uORF_class_GL_oddsratio_p2=uORF_class_GL_oddsratio_p[type!="mel_sim_conserved",]
uORF_class_GL_oddsratio_p2$class=factor(uORF_class_GL_oddsratio_p2$class,levels=c("complex","class5","class4","class3","class2","class1"))
uORF_class_GL_oddsratio_p2$type=factor(uORF_class_GL_oddsratio_p2$type,levels=c("sp4_1111","sp4_1011","sp4_0111","anc_gain","gain_in_mel","lost_in_sim"))
p=ggplot(data=uORF_class_GL_oddsratio_p2,aes(x=type,y=class,fill=log2oddsratio))+
  geom_tile()+
  scale_fill_gradient2(low = brewer.pal(11,"RdBu")[11],
                       mid = brewer.pal(11,"RdBu")[6],
                       high = brewer.pal(11,"RdBu")[1],
                       midpoint = 0,
                       limits = c(-0.67, 0.67),
                       breaks = c(-0.67, 0, 0.67),
                       labels = c("-0.67", "0", "0.67")) +
  geom_text(aes(label=p_label))+
  theme_minimal() #class5_type_GL_gainloss_oddratio.pdf
pdf("/results/class5_type_GL_gainloss_oddratio.pdf",width=8,height=4)
print(p)
dev.off()


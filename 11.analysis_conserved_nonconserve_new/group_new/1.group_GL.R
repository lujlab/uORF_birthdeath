library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(Rmisc)
library(RColorBrewer)

#group
uATG_group=fread("/data/rpkm/group/uATG_genopos_TEgroup_combined.txt",sep='\t')[,c("geno_pos","sum","group")]

#uORF matrix
uORF_matrix=fread("/data/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",sep=',',header=T)
#remove duplicated
uORF_matrix_uniq=uORF_matrix[!duplicated(uORF_matrix, by = c("seqname", "gene_pos"))]
uORF_matrix_uniq$geno_pos=paste0(uORF_matrix_uniq$seqname,"_",uORF_matrix_uniq$gene_pos)
#only in mel
uORF_matrix_uniq_mel=uORF_matrix_uniq[uORF_matrix_uniq$dm6==1,]

#merge
uORF_matrix_uniq_mel_group=merge(uORF_matrix_uniq_mel,uATG_group,by="geno_pos")

#stat GL number for each group
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

uORF_matrix_uniq_mel_group0_num=stat_GLnum(uORF_matrix_uniq_mel_group[group=="group0",])
uORF_matrix_uniq_mel_group0_num$group="group0"
uORF_matrix_uniq_mel_group1_num=stat_GLnum(uORF_matrix_uniq_mel_group[group=="group1",])
uORF_matrix_uniq_mel_group1_num$group="group1"
uORF_matrix_uniq_mel_group2_num=stat_GLnum(uORF_matrix_uniq_mel_group[group=="group2",])
uORF_matrix_uniq_mel_group2_num$group="group2"
uORF_matrix_uniq_mel_group3_num=stat_GLnum(uORF_matrix_uniq_mel_group[group=="group3",])
uORF_matrix_uniq_mel_group3_num$group="group3"
uORF_matrix_uniq_mel_total_num=stat_GLnum(uORF_matrix_uniq_mel_group)
uORF_matrix_uniq_mel_total_num$group="total"

uORF_group_GL=rbind(uORF_matrix_uniq_mel_group0_num,uORF_matrix_uniq_mel_group1_num,
                    uORF_matrix_uniq_mel_group2_num,uORF_matrix_uniq_mel_group3_num,
                    uORF_matrix_uniq_mel_total_num)

fwrite(uORF_group_GL,"/results/figS19_uORF_groupcombined_GL_mel_number.txt",sep='\t')


uORF_group_GL$total=uORF_group_GL$mel_sim_conserved+uORF_group_GL$lost_in_sim+uORF_group_GL$gain_in_mel
uORF_group_GL_p=uORF_group_GL
uORF_group_GL_oddsratio=uORF_group_GL
uORF_group_GL_ratio=uORF_group_GL
#fisher test, gain_in_mel in group0 / group0 vs. gain_in_mel in total / total = 3038 / 21498 vs. 4515 / 36565

#                   group0        total-group0                 total
#gain_in_mel        3038          4515-3038                    4515
#total-gain_in_mel  21498-3038    36565-21498-(4515-3038)      36565-4515
#total              21498         36565-21498                  36565
#a=3038
#b=21498
#c=4515
#d=36565
#fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))
#odds ratio>1  ->  more gain_in_mel in group0
for (i in c(1:4)){
  a=uORF_group_GL$sp4_1111[i]
  b=uORF_group_GL$total[i]
  c=uORF_group_GL$sp4_1111[5]
  d=uORF_group_GL$total[5]
  uORF_group_GL_p$sp4_1111[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_group_GL_oddsratio$sp4_1111[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_group_GL_ratio$sp4_1111[i]=a/b/(c/d)
  
  a=uORF_group_GL$sp4_1011[i]
  b=uORF_group_GL$total[i]
  c=uORF_group_GL$sp4_1011[5]
  d=uORF_group_GL$total[5]
  uORF_group_GL_p$sp4_1011[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_group_GL_oddsratio$sp4_1011[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_group_GL_ratio$sp4_1011[i]=a/b/(c/d)
  
  a=uORF_group_GL$sp4_0111[i]
  b=uORF_group_GL$total[i]
  c=uORF_group_GL$sp4_0111[5]
  d=uORF_group_GL$total[5]
  uORF_group_GL_p$sp4_0111[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_group_GL_oddsratio$sp4_0111[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_group_GL_ratio$sp4_0111[i]=a/b/(c/d)
  
  a=uORF_group_GL$mel_sim_conserved[i]
  b=uORF_group_GL$total[i]
  c=uORF_group_GL$mel_sim_conserved[5]
  d=uORF_group_GL$total[5]
  uORF_group_GL_p$mel_sim_conserved[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_group_GL_oddsratio$mel_sim_conserved[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_group_GL_ratio$mel_sim_conserved[i]=a/b/(c/d)
  
  a=uORF_group_GL$lost_in_sim[i]
  b=uORF_group_GL$total[i]
  c=uORF_group_GL$lost_in_sim[5]
  d=uORF_group_GL$total[5]
  uORF_group_GL_p$lost_in_sim[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_group_GL_oddsratio$lost_in_sim[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_group_GL_ratio$lost_in_sim[i]=a/b/(c/d)
  
  a=uORF_group_GL$gain_in_mel[i]
  b=uORF_group_GL$total[i]
  c=uORF_group_GL$gain_in_mel[5]
  d=uORF_group_GL$total[5]
  uORF_group_GL_p$gain_in_mel[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_group_GL_oddsratio$gain_in_mel[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_group_GL_ratio$gain_in_mel[i]=a/b/(c/d)
  
  a=uORF_group_GL$anc_gain[i]
  b=uORF_group_GL$total[i]
  c=uORF_group_GL$anc_gain[5]
  d=uORF_group_GL$total[5]
  uORF_group_GL_p$anc_gain[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$p
  uORF_group_GL_oddsratio$anc_gain[i]=fisher.test(matrix(c(a,c-a,b-a,d-b-(c-a)),nrow=2))$estimate
  uORF_group_GL_ratio$anc_gain[i]=a/b/(c/d)
}
fwrite(uORF_group_GL_p[1:4,1:8],"/results/uORF_groupcombined_GL_mel_number_p.txt",sep='\t')
fwrite(uORF_group_GL_oddsratio[1:4,1:8],"/results/uORF_groupcombined_GL_mel_number_oddsratio.txt",sep='\t')
#fwrite(uORF_group_GL_ratio[1:4,1:8],"uORF_groupcombined_GL_mel_number_enrichment.txt",sep='\t')

#heatmap
uORF_group_GL_p=fread("/results/uORF_groupcombined_GL_mel_number_p.txt")
uORF_group_GL_p2=melt(uORF_group_GL_p,id.vars="group",measure.vars=c("sp4_1111","sp4_1011","sp4_0111","mel_sim_conserved","lost_in_sim","gain_in_mel","anc_gain"))
colnames(uORF_group_GL_p2)=c("group","type","p")

uORF_group_GL_oddsratio=fread("/results/uORF_groupcombined_GL_mel_number_oddsratio.txt")
uORF_group_GL_oddsratio2=melt(uORF_group_GL_oddsratio,id.vars="group",measure.vars=c("sp4_1111","sp4_1011","sp4_0111","mel_sim_conserved","lost_in_sim","gain_in_mel","anc_gain"))
colnames(uORF_group_GL_oddsratio2)=c("group","type","oddsratio")

uORF_group_GL_oddsratio_p=merge(uORF_group_GL_oddsratio2,uORF_group_GL_p2,by=c("group","type"))
uORF_group_GL_oddsratio_p$p_label=""
uORF_group_GL_oddsratio_p[p<0.05,]$p_label="*"
uORF_group_GL_oddsratio_p[p<0.01,]$p_label="**"
uORF_group_GL_oddsratio_p[p<0.001,]$p_label="***"
uORF_group_GL_oddsratio_p$log2oddsratio=log2(uORF_group_GL_oddsratio_p$oddsratio)

#ggplot(data=uORF_group_GL_oddsratio_p,aes(x=type,y=group,fill=log2ES))+
#  geom_tile()+
#  scale_fill_gradient2(low = brewer.pal(11,"RdBu")[11],
#                       mid = brewer.pal(11,"RdBu")[6],
#                       high = brewer.pal(11,"RdBu")[1],
#                       midpoint = 0,
#                       limits = c(-0.5, 0.5),
#                       breaks = c(-0.5, 0, 0.5),
#                       labels = c("-0.5", "0", "0.5")) +
#  geom_text(aes(label=p_label))+
#  theme_minimal()


uORF_group_GL_oddsratio_p2=uORF_group_GL_oddsratio_p[type!="mel_sim_conserved",]
uORF_group_GL_oddsratio_p2$group=factor(uORF_group_GL_oddsratio_p2$group,levels=c("group3","group2","group1","group0"))
uORF_group_GL_oddsratio_p2$type=factor(uORF_group_GL_oddsratio_p2$type,levels=c("sp4_1111","sp4_1011","sp4_0111","anc_gain","gain_in_mel","lost_in_sim"))

pdf("/results/figS19group_type_GL_gainloss_oddsratio.pdf",width=6.6,height=2.5)
ggplot(data=uORF_group_GL_oddsratio_p2,aes(x=type,y=group,fill=log2oddsratio))+
  geom_tile()+
  scale_fill_gradient2(low = brewer.pal(11,"RdBu")[11],
                       mid = brewer.pal(11,"RdBu")[6],
                       high = brewer.pal(11,"RdBu")[1],
                       midpoint = 0,
                       limits = c(-0.75, 0.75),
                       breaks = c(-0.75, 0, 0.75),
                       labels = c("-0.75", "0", "0.75")) +
  geom_text(aes(label=p_label))+
  theme_minimal() #figS19group_type_GL_gainloss_oddsratio.pdf
dev.off()

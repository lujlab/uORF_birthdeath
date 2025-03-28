library(data.table)
# uORF matrix
uORF_matrix_27species <- fread("/data/uORF_matrix_dm6_pacbiosim_27species_update_noCDS.csv",sep=',',header=T)
# canon
dm6_canon <- fread("/data/dmel-all-r6.04.canonical.transcripts.txt",data.table=F,header=F)
colnames(dm6_canon) <- c("geneName","geneID","FBtr","len")
# geneInfo
dm6_geneInfo <- fread("/data/dmel6.gtf.geneInfo",data.table=F,header=T)
dm6_geneInfo_coding <- dm6_geneInfo[dm6_geneInfo$cds_len>0,]
dm6_geneInfo_coding_canon <- dm6_geneInfo_coding[dm6_geneInfo_coding$tx_name %in% dm6_canon$FBtr,]
#
tmp <- uORF_matrix_27species
tmp <- tmp[tmp$transcriptID %in% dm6_canon$FBtr,]
tmp$gainloss <- "none"
tmp$gainloss[tmp$dm6==1 & tmp$PacBioSim==1 & tmp$droYak3==0 & tmp$droAna3==0] <- "anc_gain"
tmp$gainloss[tmp$dm6==0 & tmp$PacBioSim==0 & tmp$droYak3==1 & tmp$droAna3==1] <- "anc_loss"
tmp$gainloss[tmp$dm6==1 & tmp$PacBioSim==0 & tmp$droYak3==0] <- "dmel_gain"
tmp$gainloss[tmp$dm6==0 & tmp$PacBioSim==1 & tmp$droYak3==1] <- "dmel_loss"
tmp$gainloss[tmp$dm6==0 & tmp$PacBioSim==1 & tmp$droYak3==0] <- "dsim_gain"
tmp$gainloss[tmp$dm6==1 & tmp$PacBioSim==0 & tmp$droYak3==1] <- "dsim_loss"
# unique genes
anc_gain <- unique(tmp$transcriptID[tmp$gainloss=="anc_gain"]) 
anc_loss <- unique(tmp$transcriptID[tmp$gainloss=="anc_loss"]) 
dmel_gain <- unique(tmp$transcriptID[tmp$gainloss=="dmel_gain"]) 
dmel_loss <- unique(tmp$transcriptID[tmp$gainloss=="dmel_loss"]) 
dsim_gain <- unique(tmp$transcriptID[tmp$gainloss=="dsim_gain"]) 
dsim_loss <- unique(tmp$transcriptID[tmp$gainloss=="dsim_loss"]) 

## bootstrap
# Anc loss and Dmel gain
l1 <- 1:1000
l1 <- sapply(l1,function(x){
  FBtr1 <- 1:length(anc_loss)
  for(i in 1:length(FBtr1)){
    FBtr <- anc_loss[i]
    utr5_len <- dm6_geneInfo_coding_canon$utr5_len[dm6_geneInfo_coding_canon$tx_name==FBtr]
    FBtr1[i] <- sample(dm6_geneInfo_coding_canon$tx_name[dm6_geneInfo_coding_canon$utr5_len>=utr5_len-20 & dm6_geneInfo_coding_canon$utr5_len<=utr5_len+20],1)
  }
  length(intersect(FBtr1,dmel_gain))
})
sum(l1>82)/length(l1)
mean(l1) #52
quantile(l1,c(0.025,0.975)) # 41    62 

# Anc gain & Dmel loss
l2 <- 1:1000
l2 <- sapply(l2,function(x){
  FBtr1 <- 1:length(anc_gain)
  for(i in 1:length(FBtr1)){
    FBtr <- anc_gain[i]
    utr5_len <- dm6_geneInfo_coding_canon$utr5_len[dm6_geneInfo_coding_canon$tx_name==FBtr]
    FBtr1[i] <- sample(dm6_geneInfo_coding_canon$tx_name[dm6_geneInfo_coding_canon$utr5_len>=utr5_len-20 & dm6_geneInfo_coding_canon$utr5_len<=utr5_len+20],1)
  }
  length(intersect(FBtr1,dmel_loss))
})
sum(l2>231)/length(l2)
mean(l2) #133.43
quantile(l2,c(0.025,0.975)) # 115 152 

# plot
# Anc. loss and Dmel gain
pdf("/results/fig3c.pdf",width=4,height=4)
hist(l1,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. loss and D.mel gain",
     xlim=c(30,90),ylim=c(0,100),axes=F)
abline(v=82,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
text(74,45,"P < 0.001",cex=1.5,font=3,xpd=T)
text(65,80,"Expected",cex=1.5,xpd=T)
text(74,60,"Observed",cex=1.5,xpd=T)
dev.off()
# Anc. gain & Dmel loss
pdf("/results/fig3d.pdf",width=4,height=4)
hist(l2,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. gain and D.mel loss",
     xlim=c(100,240),ylim=c(0,120),axes=F)
abline(v=231,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
#axis(1,c(90,120),c("",""),lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
text(200,75,"P < 0.001",cex=1.5,font=3,xpd=T)
text(140,110,"Expected",cex=1.5,xpd=T)
text(200,90,"Observed",cex=1.5,xpd=T)
dev.off()
fwrite(data.table(l1),"/results/loss_gain_shuffle_out")
fwrite(data.table(l2),"/results/gain_loss_shuffle_out")

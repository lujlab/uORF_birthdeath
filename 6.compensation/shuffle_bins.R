library(data.table)
setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/12.pos_to_cds/2.compensation")
#canonical=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/from_dyg/dmel-all-r6.04.canonical.transcripts.txt",header=F,sep='\t')
uORF_matrix_canonical=fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/8.new_27_mfa/PacBioSim/position_cds/uORF_matrix_dm6_PacBioSim_27species_update_noCDS_postoCDS_canonical.csv",sep=',',header=T)

# canon
dm6_canon <- fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/from_dyg/dmel-all-r6.04.canonical.transcripts.txt",data.table=F,header=F)
colnames(dm6_canon) <- c("geneName","geneID","FBtr","len")

# geneInfo
dm6_geneInfo <- fread("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/from_dyg/dmel6.gtf.geneInfo",data.table=F,header=T)
dm6_geneInfo_coding <- dm6_geneInfo[dm6_geneInfo$cds_len>0,]
dm6_geneInfo_coding_canon <- dm6_geneInfo_coding[dm6_geneInfo_coding$tx_name %in% dm6_canon$FBtr,]



#### bootstrap for all transcript ####
tmp=uORF_matrix_canonical
tmp$gainloss <- "none"
tmp$gainloss[tmp$dm6==1 & tmp$PacBioSim==1 & tmp$droYak3==0 & tmp$droAna3==0] <- "anc_gain"
tmp$gainloss[tmp$dm6==0 & tmp$PacBioSim==0 & tmp$droYak3==1 & tmp$droAna3==1] <- "anc_loss"
tmp$gainloss[tmp$dm6==1 & tmp$PacBioSim==0 & tmp$droYak3==0] <- "dmel_gain"
tmp$gainloss[tmp$dm6==0 & tmp$PacBioSim==1 & tmp$droYak3==1] <- "dmel_loss"
tmp$gainloss[tmp$dm6==0 & tmp$PacBioSim==1 & tmp$droYak3==0] <- "dsim_gain"
tmp$gainloss[tmp$dm6==1 & tmp$PacBioSim==0 & tmp$droYak3==1] <- "dsim_loss"
# unique genes
anc_gain <- unique(tmp$transcriptID[tmp$gainloss=="anc_gain"]) # 1970
anc_loss <- unique(tmp$transcriptID[tmp$gainloss=="anc_loss"]) # 203
dmel_gain <- unique(tmp$transcriptID[tmp$gainloss=="dmel_gain"]) # 1604
dmel_loss <- unique(tmp$transcriptID[tmp$gainloss=="dmel_loss"]) # 553
dsim_gain <- unique(tmp$transcriptID[tmp$gainloss=="dsim_gain"]) # 1364
dsim_loss <- unique(tmp$transcriptID[tmp$gainloss=="dsim_loss"]) # 451

# bootstrap 
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
hist(l1,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. loss and D.mel gain",
     xlim=c(30,90),ylim=c(0,100),axes=F)
abline(v=82,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
text(74,45,"P < 0.001",cex=1.5,font=3,xpd=T)
text(65,80,"Expected",cex=1.5,xpd=T)
text(74,60,"Observed",cex=1.5,xpd=T)

# Anc. gain & Dmel loss
hist(l2,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. gain and D.mel loss",
     xlim=c(100,240),ylim=c(0,120),axes=F)
abline(v=231,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
#axis(1,c(90,120),c("",""),lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
text(200,75,"P < 0.001",cex=1.5,font=3,xpd=T)
text(140,110,"Expected",cex=1.5,xpd=T)
text(200,90,"Observed",cex=1.5,xpd=T)

#fwrite(data.table(l1),"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/9.compensation/loss_gain_shuffle_out")
#fwrite(data.table(l2),"/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/9.compensation/gain_loss_shuffle_out")


#### bootstrap for <100nt regions #### 
tmp=uORF_matrix_canonical[uORF_matrix_canonical$pos_to_cds<=100,]
tmp$gainloss <- "none"
tmp$gainloss[tmp$dm6==1 & tmp$PacBioSim==1 & tmp$droYak3==0 & tmp$droAna3==0] <- "anc_gain"
tmp$gainloss[tmp$dm6==0 & tmp$PacBioSim==0 & tmp$droYak3==1 & tmp$droAna3==1] <- "anc_loss"
tmp$gainloss[tmp$dm6==1 & tmp$PacBioSim==0 & tmp$droYak3==0] <- "dmel_gain"
tmp$gainloss[tmp$dm6==0 & tmp$PacBioSim==1 & tmp$droYak3==1] <- "dmel_loss"
tmp$gainloss[tmp$dm6==0 & tmp$PacBioSim==1 & tmp$droYak3==0] <- "dsim_gain"
tmp$gainloss[tmp$dm6==1 & tmp$PacBioSim==0 & tmp$droYak3==1] <- "dsim_loss"
# unique genes
anc_gain <- unique(tmp$transcriptID[tmp$gainloss=="anc_gain"]) # 713
anc_loss <- unique(tmp$transcriptID[tmp$gainloss=="anc_loss"]) # 59
dmel_gain <- unique(tmp$transcriptID[tmp$gainloss=="dmel_gain"]) # 536
dmel_loss <- unique(tmp$transcriptID[tmp$gainloss=="dmel_loss"]) # 183
dsim_gain <- unique(tmp$transcriptID[tmp$gainloss=="dsim_gain"]) # 508
dsim_loss <- unique(tmp$transcriptID[tmp$gainloss=="dsim_loss"]) # 127

# bootstrap 
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
sum(l1>5)/length(l1) #0.041
mean(l1) #2.549
quantile(l1,c(0.025,0.975)) # 0   6 

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
sum(l2>29)/length(l2) #0
mean(l2) #10.595
quantile(l2,c(0.025,0.975)) # 5    17 

# plot
# Anc. loss and Dmel gain
hist(l1,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. loss and D.mel gain",
     xlim=c(0,10),ylim=c(0,300),axes=F)
abline(v=5,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
text(6,150,"P = 0.041",cex=1.5,font=3,xpd=T)
text(2,280,"Expected",cex=1.5,xpd=T)
text(6,200,"Observed",cex=1.5,xpd=T)
# Anc. gain & Dmel loss
hist(l2,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. gain and D.mel loss",
     xlim=c(0,35),ylim=c(0,150),axes=F)
abline(v=29,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
#axis(1,c(90,120),c("",""),lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
#text(200,75,"P < 0.001",cex=1.5,font=3,xpd=T)
#text(140,10,"Expected",cex=1.5,xpd=T)
#text(200,90,"Observed",cex=1.5,xpd=T)

fwrite(data.table(l1),"./shuffle/loss_gain_shuffle_out_100nt")
fwrite(data.table(l2),"./shuffle/gain_loss_shuffle_out_100nt")



#### bootstrap for <200nt regions #### 
tmp=uORF_matrix_canonical[uORF_matrix_canonical$pos_to_cds<=200,]
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

# bootstrap 
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
sum(l1>14)/length(l1) #0.018
mean(l1) #8.391
quantile(l1,c(0.025,0.975)) # 4    14 

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
sum(l2>66)/length(l2) #0
mean(l2) #31.112
quantile(l2,c(0.025,0.975)) # 21.000 41.025 

# plot
# Anc. loss and Dmel gain
hist(l1,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. loss and D.mel gain",
     xlim=c(0,20),ylim=c(0,200),axes=F)
abline(v=14,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
text(15,150,"P = 0.018",cex=1.5,font=3,xpd=T)
#text(2,280,"Expected",cex=1.5,xpd=T)
#text(6,200,"Observed",cex=1.5,xpd=T)
# Anc. gain & Dmel loss
hist(l2,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. gain and D.mel loss",
     xlim=c(0,70),ylim=c(0,120),axes=F)
abline(v=66,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
#axis(1,c(90,120),c("",""),lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
#text(200,75,"P < 0.001",cex=1.5,font=3,xpd=T)
#text(140,10,"Expected",cex=1.5,xpd=T)
#text(200,90,"Observed",cex=1.5,xpd=T)

fwrite(data.table(l1),"./shuffle/loss_gain_shuffle_out_200nt")
fwrite(data.table(l2),"./shuffle/gain_loss_shuffle_out_200nt")



#### bootstrap for <300nt regions #### 
tmp=uORF_matrix_canonical[uORF_matrix_canonical$pos_to_cds<=300,]
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

# bootstrap 
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
sum(l1>24)/length(l1) #0.01
mean(l1) #15.439
quantile(l1,c(0.025,0.975)) # 9    23 

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
sum(l2>107)/length(l2) #0
mean(l2) #55.088
quantile(l2,c(0.025,0.975)) # 42.975 69.000

# plot
# Anc. loss and Dmel gain
hist(l1,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. loss and D.mel gain",
     xlim=c(0,30),ylim=c(0,200),axes=F)
abline(v=24,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
text(15,150,"P = 0.01",cex=1.5,font=3,xpd=T)
#text(2,280,"Expected",cex=1.5,xpd=T)
#text(6,200,"Observed",cex=1.5,xpd=T)
# Anc. gain & Dmel loss
hist(l2,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. gain and D.mel loss",
     xlim=c(20,120),ylim=c(0,150),axes=F)
abline(v=107,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
#axis(1,c(90,120),c("",""),lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
#text(200,75,"P < 0.001",cex=1.5,font=3,xpd=T)
#text(140,10,"Expected",cex=1.5,xpd=T)
#text(200,90,"Observed",cex=1.5,xpd=T)
#5.82 5.27

fwrite(data.table(l1),"./shuffle/loss_gain_shuffle_out_300nt")
fwrite(data.table(l2),"./shuffle/gain_loss_shuffle_out_300nt")





#### bootstrap for <500nt regions #### 
tmp=uORF_matrix_canonical[uORF_matrix_canonical$pos_to_cds<=500,]
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

# bootstrap 
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
sum(l1>49)/length(l1) #0.01
mean(l1) #15.439
quantile(l1,c(0.025,0.975)) # 9    23 

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
sum(l2>173)/length(l2) #0
mean(l2) #55.088
quantile(l2,c(0.025,0.975)) # 42.975 69.000

# plot
# Anc. loss and Dmel gain
hist(l1,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. loss and D.mel gain",
     xlim=c(0,55),ylim=c(0,150),axes=F)
abline(v=24,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
text(15,150,"P = 0.01",cex=1.5,font=3,xpd=T)
#text(2,280,"Expected",cex=1.5,xpd=T)
#text(6,200,"Observed",cex=1.5,xpd=T)
# Anc. gain & Dmel loss
hist(l2,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. gain and D.mel loss",
     xlim=c(50,200),ylim=c(0,150),axes=F)
abline(v=107,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
#axis(1,c(90,120),c("",""),lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
#text(200,75,"P < 0.001",cex=1.5,font=3,xpd=T)
#text(140,10,"Expected",cex=1.5,xpd=T)
#text(200,90,"Observed",cex=1.5,xpd=T)
#5.82 5.27

fwrite(data.table(l1),"./shuffle/loss_gain_shuffle_out_500nt")
fwrite(data.table(l2),"./shuffle/gain_loss_shuffle_out_500nt")

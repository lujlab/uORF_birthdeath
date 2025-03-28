library(data.table)
l1=fread("/data/loss_gain_shuffle_out_300nt")
l2=fread("/data/gain_loss_shuffle_out_300nt")

pdf("/results/figS21a.pdf",width=6,height=6)
hist(l1$l1,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. loss and D.mel gain",
     xlim=c(0,30),ylim=c(0,200),axes=F)
abline(v=24,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
text(15,150,"P = 0.01",cex=1.5,font=3,xpd=T)
#text(2,280,"Expected",cex=1.5,xpd=T)
#text(6,200,"Observed",cex=1.5,xpd=T)
dev.off()

pdf("/results/figS21b.pdf",width=6,height=6)
# Anc. gain & Dmel loss
hist(l2$l2,breaks=30,col="aquamarine",xlab="Number of genes",main="Genes with Anc. gain and D.mel loss",
     xlim=c(20,120),ylim=c(0,150),axes=F)
abline(v=107,col="red",lty=5,lwd=2)
axis(1,lwd=4,lwd.ticks=4,cex.axis=2)
#axis(1,c(90,120),c("",""),lwd=4,lwd.ticks=4,cex.axis=2)
axis(2,lwd=4,lwd.ticks=4,cex.axis=2)
text(100,75,"P < 0.001",cex=1.5,font=3,xpd=T)
dev.off()
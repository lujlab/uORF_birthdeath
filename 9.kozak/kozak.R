library(data.table)
library(parallel)
library(ggplot2)
library(pheatmap)
suppressPackageStartupMessages(library(Biostrings))
library(RColorBrewer)
args <- commandArgs(trailingOnly = TRUE)

cds.kozak <- fread(args[1],data.table=F,header=T)[,2]
cds.kozak <- cds.kozak[nchar(cds.kozak)==10]
cds.kozak <- cds.kozak[subseq(cds.kozak,7,9)=='ATG']
cds.pfm <- consensusMatrix(cds.kozak)
toPWM <- function(x){apply(x,2,function(x){x/sum(x)})}
cds.pwm <- toPWM(cds.pfm[c(1:4), ]) #for dm6
#cds.pwm <- toPWM(cds.pfm[c(1:3,5), ]) #for sim1
PwmPadding <- function(pwm){
    pwm <- rbind(pwm, N = rep(0.25, 10))
    return(pwm[, c(1:6, 10)])
}
cds.pwm <- PwmPadding(cds.pwm)

uorf.kozak <- fread(input="cat /dev/stdin",data.table=F,header=F)[,1]
uorf.kozak0 <- uorf.kozak
uorf.kozak <- xscat(subseq(uorf.kozak, s = 1, e = 6), subseq(uorf.kozak, s = 10, e = 10))
sc <- sapply(strsplit(as.character(uorf.kozak), '', fixed = TRUE), function(x){
        x[!x %in% c('A', 'C', 'G', 'T')] <- 'N'
        sum(diag(log2(cds.pwm[x, ]/0.25)))
})
uorf.kozak <- data.frame(seq=uorf.kozak0,kozak_score=sc)
fwrite(uorf.kozak,file="",col.names=F,row.names=F,sep="\t",quote=F)

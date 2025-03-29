setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/9.compensation/MK_test_1356/")
#setwd("/lustre/user/lulab/luoshq/syq/project/uorf_MK/asymptoticMK_wtw")
library(data.table)

## read polymorphic sites and incooperate ancestral base ####################################################
anc.base <- fread('./pop_AF/1356_snp_dm6_anc.txt', header = F)
colnames(anc.base) <- c('id', 'dm6', 'sim', 'type')
setkey(anc.base, id)
#add gene id 
gene_tr<-fread("dmel-all-geneID_transcriptID-r6.04.tr_sorted",header = F,stringsAsFactors = F)
names(gene_tr)<-c("isoform","gene2")

CorrectAnc <- function(dtt.raw){
  dtt <- dtt.raw[id %in% anc.base$id] # syq: filter SNPs that: 1) not correspoding site in Dsim; 2) which is lowercase in Dmel
  dtt[anc.base, sim := i.sim, on = 'id'] # syq: if ref != sim & alt != sec, it can't determine which allele is derived.  
  dtt <- dtt[ref == sim | alt == sim]
  dtt[, reftype := ifelse(ref == sim, 'anc', 'der')]
  dtt[, daf := ifelse(reftype == 'anc', af, 1 - af)]
  dtt[, anc_seq := ifelse(reftype == 'anc', ref_seq, alt_seq)]
  dtt[, der_seq := ifelse(reftype == 'anc', alt_seq, ref_seq)]
  return(dtt)
}
GetPolyAUGLoss <- function(x){
  poly.raw.region <- fread(x)
  poly.region <- CorrectAnc(poly.raw.region)
  poly.region.aug.loss <- poly.region[anc_seq == 'ATG'][!duplicated(id)]
  return(poly.region.aug.loss)
}

poly.utr5.aug.loss <- GetPolyAUGLoss('poly_diff_utr5.txt')
poly.intron.strict.aug.loss <- GetPolyAUGLoss('poly_diff_short_intron_strict_neutral.txt')


## read fixed AUGloss and remove polymorphic sites ####################################################
all.pol.sites <- fread('pop_AF/1356_snp_dm6_all_coods_nohead.txt',header=F)
#colnames(all.pol.sites)=c("CHROM","POS","ID","AF")
setnames(all.pol.sites, c('chrm', 'posn', 'id', 'af'))
all.pol.sites[, chrmpos := paste(chrm, posn, sep = '_')]
GetFixedTriplet <- function(x){
  dtt <- fread(x)
  dtt.final <- dtt[
    !paste(chrm, posn1, sep = '_') %in% all.pol.sites$chrmpos &
      !paste(chrm, posn2, sep = '_') %in% all.pol.sites$chrmpos &
      !paste(chrm, posn3, sep = '_') %in% all.pol.sites$chrmpos]
  return(dtt.final)
}

# /lustre/user/lulab/luoshq/syq/project/uorf_MK/asymptoticMK_wtw
tmp=fread('fixed_diff_utr5.txt')
tmphead=colnames(tmp)
tmp2=merge(tmp,gene_tr,by="isoform")
tmp2$gene=tmp2$gene2
tmp3=tmp2[,..tmphead]
fwrite(tmp3, 'fixed_diff_utr5_gene.txt', sep = '\t')

fix.utr5 <- fread('fixed_diff_utr5_gene.txt')
fix.utr5[, id := paste(chrm, posn1, posn2, posn3, sep = '_')]
fix.utr5.site.bed <- melt(fix.utr5, id.vars = c('chrm', 'id', 'strand'),measure.vars = c('posn1', 'posn2', 'posn3'))
fix.utr5.site.bed <- fix.utr5.site.bed[, .(chrm = sub('^chr', '', chrm), start = value - 1, end = value,name = id, score = substr(variable, 5, 5), strand)] #syq: there a lot duplication 
fwrite(fix.utr5.site.bed, 'fixed_diff_utr5_site.bed', col.names = FALSE, sep = '\t')
# bedtools intersect -a fixed_diff_utr5_site.bed -b dmel-all-r6.04.gtf.CDS.bed >fixed_diff_utr5_site_cds.ov.bed 
fix.utr5.cds.ov <- fread("./fixed_diff_utr5_site_cds.ov.bed")
fix.utr5.nocds <- fix.utr5[!id %in% fix.utr5.cds.ov$V4]
fwrite(fix.utr5.nocds, file = 'fixed_diff_utr5_nocds.txt', sep = '\t')
fix.utr5.aug <- fix.utr5.nocds[toupper(dm6_seq) == 'ATG'][
  dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0][
    !duplicated(id)]
fix.utr5.augloss <- fix.utr5.nocds[
  dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0]
fix.utr5.augloss[, yak_seq := gsub('#|-', '', tstrsplit(alignment, ';')[[3]])]
fix.utr5.augloss <- fix.utr5.augloss[dm6_seq != 'ATG' & sim_seq == 'ATG' & yak_seq == 'ATG']
fix.utr5.augloss <- fix.utr5.augloss[!duplicated(id)]
fwrite(fix.utr5.aug, 'fixed_aug_utr5.txt', sep = '\t')
fwrite(fix.utr5.augloss, 'fixed_augloss_utr5.txt')

fix.utr5.aug.loss <- GetFixedTriplet('./fixed_augloss_utr5.txt')

# zcat short_intron_strict_neutral_syq.maf.gz | python find_newly_fixed_triplet_dsim.py intron > fixed_diff_intron_strict_neutral_dsim.txt
#zcat short_intron_strict_neutral.maf.gz | python find_newly_fixed_triplet_lcl.py intron > fixed_diff_intron_strict_neutral.txt
fix.si.strict <- fread('fixed_diff_intron_strict_neutral.txt', sep = '\t', quote = '')
fix.si.strict.aug <- fix.si.strict[toupper(dm6_seq) == 'ATG'][dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0] # 如果diff_in_3sp！=0，则无法判断突变是否是在 dmel lineage 上发生的
fix.si.strict.augloss <- fix.si.strict[dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0]
fix.si.strict.augloss[, yak_seq := gsub('#|-', '', tstrsplit(alignment, ';')[[3]])]
fix.si.strict.augloss <- fix.si.strict.augloss[dm6_seq != 'ATG' & sim_seq == 'ATG' & yak_seq == 'ATG']
fwrite(fix.si.strict.aug, 'fixed_aug_intron_strict_neutral.txt', sep = '\t')
fwrite(fix.si.strict.augloss, 'fixed_augloss_intron_strict_neutral.txt', sep = '\t')
fix.intron.strict.aug.loss <- GetFixedTriplet('./fixed_augloss_intron_strict_neutral.txt')

## MKtest for AUG loss #####
source('./asymptoticMK_src/asymptoticMK_local_clean.R')

K80.fix <- function(ts, tv){
  pts = ts/(ts+tv) *0.05
  qtv = tv/(ts+tv) *0.05
  w1 = 1 - 2 * pts - qtv
  w2 = 1 - 2 * qtv
  return( (-0.5 * log(w1) - 0.25 *log(w2)) * (ts + tv)/0.05)
}

RunAsymtoticMK <- function(fixed.dt, fixed.dt0, pol.dt, pol.dt0, nbin = 20, resfmt = 'df',
                           daf_low = 0.05, daf_high = 0.95, ...){
  if('cnt_ts' %in% names(fixed.dt)){
    d <- K80.fix(sum(fixed.dt$cnt_ts), sum(fixed.dt$cnt_tv))
  }else{
    d <- K80.fix(sum(fixed.dt$mut_tag == 'S'), sum(fixed.dt$mut_tag == 'V'))
  }
  if('cnt_ts' %in% names(fixed.dt0)){
    d0 <- K80.fix(sum(fixed.dt0$cnt_ts), sum(fixed.dt0$cnt_tv))
  }else{
    d0 <- K80.fix(sum(fixed.dt0$mut_tag == 'S'), sum(fixed.dt0$mut_tag == 'V'))
  }
  
  d <- round(d)
  d0 <- round(d0)
  
  pol.dt <- copy(pol.dt)
  pol.dt0 <- copy(pol.dt0)
  
  # "polymorphic sites in neutral regions were grouped into bins of equal size based on increasing derived allele frequency"
  set.seed(123)
  x1 <- sort(pol.dt0[, daf])
  x2 <- x1[ seq(1, length(x1), by = length(x1) / nbin) ]
  x2[1] <- 0
  x2[length(x2)] <- 1
  daf_breaks <- unique(x2)
  df.breaks <- data.frame(low = daf_breaks[-length(daf_breaks)], high = daf_breaks[-1])
  
  dtf.raw <- rbind(test = pol.dt, neutral = pol.dt0, idcol = 'grp', fill = TRUE)
  dtf.raw[, freq := cut(daf, breaks = daf_breaks, labels = FALSE, include.lowest = TRUE)]
  dtf <- dtf.raw[
    , .(x = df.breaks$high[findInterval(median(daf), daf_breaks, all.inside = TRUE)],
        p = sum(grp == 'test'), p0 = sum(grp == 'neutral')),
    by = 'freq'][order(freq)]
  dtf[, freq := NULL]
  
  res <- asymptoticMK(d0, d, df = dtf, xlow = daf_low, xhigh = daf_high, ...)
  res$d <- d
  res$d0 <- d0
  res$p <- nrow(pol.dt[daf >= daf_low & daf <= daf_high])
  res$p0 <- nrow(pol.dt0[daf >= daf_low & daf <= daf_high])
  res$alpha_raw <- 1 - (d0/res$p0) / (d/res$p)
  if(resfmt == 'df'){
    print(res)
    invisible(dtf)
  }else{
    invisible(res)
  }
}

res<-RunAsymtoticMK(fix.utr5.aug.loss, fix.intron.strict.aug.loss,
                    poly.utr5.aug.loss, poly.intron.strict.aug.loss,
                    daf_low = 0.05, daf_high = 0.95, nbin = 20, output = 'none')
#model         a          b        c alpha_asymptotic    CI_low   CI_high alpha_original   d  d0   p  p0 alpha_raw
#1 exponential 0.6242544 -0.4494757 11.22187        0.6242484 0.5977541 0.6496715       0.539027 617 152 994 612 0.5998771

#5*5
#barplot(c(0.5831,0.6316),space = 0.8,las=1,ylim = c(0,0.8))

## MKtest for AUG loss in different gene set:compensatory uORF #####
gene_tr<-fread("dmel-all-geneID_transcriptID-r6.04.tr_sorted",header = F,stringsAsFactors = F)
names(gene_tr)<-c("tr","gene")
merge(poly.utr5.aug.loss,gene_tr,by.x="feature",by.y="tr")->poly.utr5.aug.loss2

#hh<-readRDS("/lustre/user/lulab/luoshq/syq/project/uorf_MK/asymptoticMK_wtw/FBtr_uORF_gain_loss_v2_Dsim.rds")
hh<-readRDS("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/9.compensation/uORF_gain_loss_trlist.rds",refhook = NULL)
for (i in 1:length(hh)) {
  x<-as.data.table(hh[[i]])
  names(x)<-"tr"
  merge(x,gene_tr,by="tr")->x
  assign(names(hh)[i],x)
}
as.data.table(intersect(anc_gain$gene,mel_loss$gene))->ancGain_melLoss
names(ancGain_melLoss)<-"gene"
mel_loss[!gene%in%anc_gain$gene]->mel_loss_323
as.data.table(intersect(anc_loss$gene,mel_gain$gene))->ancLoss_melGain
names(ancLoss_melGain)<-"gene"
mel_gain[!gene%in%anc_loss$gene]->mel_gain_1522

## total mel loss:DYG 的dmel_loss没有全部在 SYQ 的fix.utr5.aug.loss找到，因为1）fix.utr5.aug.loss过滤了多态性位点，dmel_loss没有；2）
res<-RunAsymtoticMK(fix.utr5.aug.loss[gene%in%mel_loss$gene], fix.intron.strict.aug.loss,
                    poly.utr5.aug.loss2[gene%in%mel_loss$gene], poly.intron.strict.aug.loss,
                    daf_low = 0.05, daf_high = 0.95, nbin = 50, output = 'none')
#model         a           b  c alpha_asymptotic    CI_low   CI_high alpha_original   d  d0   p
#1 linear 0.8746758 -0.03257846 NA        0.8420973 0.7594218 0.9247728      0.8622748 318 152 179
#p0 alpha_raw
#1 612 0.86019655

## mel_loss & anc_gain
res<-RunAsymtoticMK(fix.utr5.aug.loss[gene%in%ancGain_melLoss$gene], fix.intron.strict.aug.loss,
                    poly.utr5.aug.loss2[gene%in%ancGain_melLoss$gene], poly.intron.strict.aug.loss,
                    daf_low = 0.05, daf_high = 0.95, nbin = 50, output = 'none')
#model         a          b  c alpha_asymptotic    CI_low   CI_high alpha_original   d  d0   p
#1 linear 0.8049147 0.06555686 NA        0.8704716 0.7615535 0.9793897      0.8298689 159 152 108
#p0 alpha_raw
#1 612 0.8312986

## mel_loss without anc_gain
res<-RunAsymtoticMK(fix.utr5.aug.loss[gene%in%mel_loss_323$gene], fix.intron.strict.aug.loss,
                    poly.utr5.aug.loss2[gene%in%mel_loss_323$gene], poly.intron.strict.aug.loss,
                    daf_low = 0.05, daf_high = 0.95, nbin = 50, output = 'none')
#model         a          b  c alpha_asymptotic    CI_low  CI_high alpha_original   d  d0  p  p0
#1 linear 0.9447841 -0.1298968 NA        0.8148873 0.7406635 0.889111       0.895339 160 152 71 612
#alpha_raw
#1 0.8897876

## MKtest for all triplets loss(total 64): 5' UTR vs. short intron ###########################################
## fixed triplets loss
TRIPLETS <- do.call(paste0, do.call(expand.grid, list(c('A', 'C', 'G', 'T'))[rep(1, 3)]))

fix.si.strict<- GetFixedTriplet('./fixed_diff_intron_strict_neutral.txt')
fix.si.strict[,yak_seq := gsub('#|-', '', tstrsplit(alignment, ';')[[3]])] # original: sub ??
fix.si.strict.triplet <- lapply(TRIPLETS, function(x){
  res <- fix.si.strict[toupper(dm6_seq) != x & toupper(sim_seq) == x & toupper(yak_seq) == x][
    dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0]
  return(res)
})

fix.utr5.nocds <- GetFixedTriplet('./fixed_diff_utr5_nocds.txt')
fix.utr5.nocds[,yak_seq := gsub('#|-', '', tstrsplit(alignment, ';')[[3]])]
fix.utr5.triplet <- lapply(TRIPLETS, function(x){
  res <- fix.utr5.nocds[toupper(dm6_seq) != x & toupper(sim_seq) == x & toupper(yak_seq) == x][
    dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0][
      !duplicated(id)]
  return(res)
})

## polymorphic triplets loss
poly.utr5 <- CorrectAnc(fread('poly_diff_utr5.txt'))
poly.intron.strict <- CorrectAnc(fread('poly_diff_short_intron_strict_neutral.txt'))

poly.utr5.triplet <- lapply(TRIPLETS, function(x){
  poly.utr5[anc_seq == x][!duplicated(id)]
})

poly.intron.strict.triplet <- lapply(TRIPLETS, function(x){
  poly.intron.strict[anc_seq == x][!duplicated(id)]
})

## MKtest for all triplets loss
alpha.all.triplets.loss <- mapply(
  function(df, df0, dp, dp0){
    res <- RunAsymtoticMK(
      df, df0, dp, dp0,
      daf_low = 0.1, daf_high = 0.9, nbin = 10, output = 'table', resfmt = 'alpha')
    res$p <- sum(dp$daf > 0.1 & dp$daf < 0.9)
    res$p0 <- sum(dp0$daf > 0.1 & dp0$daf < 0.9)
    return(res)
  },
  fix.utr5.triplet, fix.si.strict.triplet,
  poly.utr5.triplet, poly.intron.strict.triplet, SIMPLIFY = FALSE)

names(alpha.all.triplets) <- TRIPLETS
alpha.all.triplets.df <- rbindlist(alpha.all.triplets, idcol = 'triplet')
alpha.all.triplets.df[, ratio.dp := d / p]
fwrite(alpha.all.triplets.df,"alpha_all_TripletsLoss.txt",col.names =T,row.names = F,sep="\t",quote = F)

## MKtest for all triplets gain: 5' UTR vs. short intron. syq: ???###########################################
TRIPLETS <- do.call(paste0, do.call(expand.grid, list(c('A', 'C', 'G', 'T'))[rep(1, 3)]))
## fixed triplets
fix.si.strict <- GetFixedTriplet('./fixed_diff_intron_strict_neutral.txt')
fix.utr5.nocds <- GetFixedTriplet('./fixed_diff_utr5_nocds.txt')

fix.si.strict.triplet <- lapply(TRIPLETS, function(x){
  res <- fix.si.strict[toupper(dm6_seq) == x][
    dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0]
})
fix.utr5.triplet <- lapply(TRIPLETS, function(x){
  res <- fix.utr5.nocds[toupper(dm6_seq) == x][
    dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0][
      !duplicated(id)]
})

## polymorphic triplets
poly.utr5 <- CorrectAnc(fread('poly_diff_utr5.txt'))
poly.intron.strict <- CorrectAnc(fread('poly_diff_short_intron_strict_neutral.txt'))

poly.utr5.triplet <- lapply(TRIPLETS, function(x){
  poly.utr5[der_seq == x][!duplicated(id)]
})
poly.intron.strict.triplet <- lapply(TRIPLETS, function(x){
  poly.intron.strict[der_seq == x][!duplicated(id)]
})


RunAsymtoticMK2 <- function(fixed.dt, fixed.dt0, pol.dt, pol.dt0, nbin = 20, resfmt = 'df', ...){
  # different with RunAsymtoticMK in bins divide
  if('cnt_ts' %in% names(fixed.dt)){
    d <- K80.fix(sum(fixed.dt$cnt_ts), sum(fixed.dt$cnt_tv))
  }else{
    d <- K80.fix(sum(fixed.dt$mut_tag == 'S'), sum(fixed.dt$mut_tag == 'V'))
  }
  if('cnt_ts' %in% names(fixed.dt0)){
    d0 <- K80.fix(sum(fixed.dt0$cnt_ts), sum(fixed.dt0$cnt_tv))
  }else{
    d0 <- K80.fix(sum(fixed.dt0$mut_tag == 'S'), sum(fixed.dt0$mut_tag == 'V'))
  }
  
  d <- round(d)
  d0 <- round(d0)
  
  pol.dt <- copy(pol.dt)
  pol.dt0 <- copy(pol.dt0)
  
  pol.dt[, freq := cut(daf, breaks = seq(0, nbin) / nbin, labels = FALSE)]
  tmp1 <- pol.dt[, .(cnt = .N, daf0 = median(daf)), by = 'freq'][order(freq)]
  tmp1[, daf := (freq - 0.5) / nbin]
  
  pol.dt0[, freq := cut(daf, breaks = seq(0, nbin) / nbin, labels = FALSE)]
  tmp2 <- pol.dt0[, .(cnt = .N, daf0 = median(daf)), by = 'freq'][order(freq)]
  tmp2[, daf := (freq - 0.5) / nbin]
  
  dtf <- merge(tmp1, tmp2, by = c('freq', 'daf'))
  dtf <- dtf[, .(x = daf, p = cnt.x, p0 = cnt.y)]
  
  res <- asymptoticMK(d0, d, df = dtf, ...)
  res$d <- d
  res$d0 <- d0
  res$p <- nrow(pol.dt[daf  >= 0.1 & daf <= 0.9])
  res$p0 <- nrow(pol.dt0[daf  >= 0.1 & daf <= 0.9])
  res$alpha_raw <- 1 - (d0/res$p0) / (d/res$p)
  if(resfmt == 'df'){
    print(res)
    invisible(dtf)
  }else{
    invisible(res)
  }
}
# MKtest for all triplets
alpha.all.triplets2 <- mapply(
  function(df, df0, dp, dp0){
    res <- RunAsymtoticMK2(
      df, df0, dp, dp0,
      xlow = 0.1, xhigh = 0.9, nbin = 10, output = 'table', resfmt = 'alpha')
    res$p <- sum(dp$daf > 0.1 & dp$daf < 0.9)
    res$p0 <- sum(dp0$daf > 0.1 & dp0$daf < 0.9)
    return(res)
  },
  fix.utr5.triplet, fix.si.strict.triplet,
  poly.utr5.triplet, poly.intron.strict.triplet, SIMPLIFY = FALSE)

names(alpha.all.triplets2) <- TRIPLETS
alpha.all.triplets.df2 <- rbindlist(alpha.all.triplets2, idcol = 'triplet')
alpha.all.triplets.df2[, ratio.dp := d / p]

## MKtest for AUG gain (based on above MKtest for all triplets gain)###################################
#TRIPLETS <- do.call(paste0, do.call(expand.grid, list(c('A', 'C', 'G', 'T'))[rep(1, 3)]))
TRIPLETS<-"ATG"
## fixed triplets
fix.si.strict <- GetFixedTriplet('./fixed_diff_intron_strict_neutral.txt')
fix.utr5.nocds <- GetFixedTriplet('./fixed_diff_utr5_nocds.txt')

fix.si.strict.triplet <- lapply(TRIPLETS, function(x){
  res <- fix.si.strict[toupper(dm6_seq) == x][
    dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0]
})
fix.utr5.triplet <- lapply(TRIPLETS, function(x){
  res <- fix.utr5.nocds[toupper(dm6_seq) == x][
    dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0][
      !duplicated(id)]
})

## polymorphic triplets
poly.utr5 <- CorrectAnc(fread('poly_diff_utr5.txt'))
poly.intron.strict <- CorrectAnc(fread('poly_diff_short_intron_strict_neutral.txt'))

poly.utr5.triplet <- lapply(TRIPLETS, function(x){
  poly.utr5[der_seq == x][!duplicated(id)]
})
poly.intron.strict.triplet <- lapply(TRIPLETS, function(x){
  poly.intron.strict[der_seq == x][!duplicated(id)]
})


RunAsymtoticMK2 <- function(fixed.dt, fixed.dt0, pol.dt, pol.dt0, nbin = 20, resfmt = 'df', ...){
  # different with RunAsymtoticMK in bins divide
  if('cnt_ts' %in% names(fixed.dt)){
    d <- K80.fix(sum(fixed.dt$cnt_ts), sum(fixed.dt$cnt_tv))
  }else{
    d <- K80.fix(sum(fixed.dt$mut_tag == 'S'), sum(fixed.dt$mut_tag == 'V'))
  }
  if('cnt_ts' %in% names(fixed.dt0)){
    d0 <- K80.fix(sum(fixed.dt0$cnt_ts), sum(fixed.dt0$cnt_tv))
  }else{
    d0 <- K80.fix(sum(fixed.dt0$mut_tag == 'S'), sum(fixed.dt0$mut_tag == 'V'))
  }
  
  d <- round(d)
  d0 <- round(d0)
  
  pol.dt <- copy(pol.dt)
  pol.dt0 <- copy(pol.dt0)
  
  pol.dt[, freq := cut(daf, breaks = seq(0, nbin) / nbin, labels = FALSE)]
  tmp1 <- pol.dt[, .(cnt = .N, daf0 = median(daf)), by = 'freq'][order(freq)]
  tmp1[, daf := (freq - 0.5) / nbin]
  
  pol.dt0[, freq := cut(daf, breaks = seq(0, nbin) / nbin, labels = FALSE)]
  tmp2 <- pol.dt0[, .(cnt = .N, daf0 = median(daf)), by = 'freq'][order(freq)]
  tmp2[, daf := (freq - 0.5) / nbin]
  
  dtf <- merge(tmp1, tmp2, by = c('freq', 'daf'))
  dtf <- dtf[, .(x = daf, p = cnt.x, p0 = cnt.y)]
  
  res <- asymptoticMK(d0, d, df = dtf, ...)
  res$d <- d
  res$d0 <- d0
  res$p <- nrow(pol.dt[daf  >= 0.05 & daf <= 0.95])
  res$p0 <- nrow(pol.dt0[daf  >= 0.05 & daf <= 0.95])
  res$alpha_raw <- 1 - (d0/res$p0) / (d/res$p)
  if(resfmt == 'df'){
    print(res)
    invisible(dtf)
  }else{
    invisible(res)
  }
}
# MKtest for all triplets
alpha.all.triplets2 <- mapply(
  function(df, df0, dp, dp0){
    res <- RunAsymtoticMK2(
      df, df0, dp, dp0,
      xlow = 0.05, xhigh = 0.95, nbin = 10, output = 'table', resfmt = 'alpha')
    res$p <- sum(dp$daf > 0.05 & dp$daf < 0.95)
    res$p0 <- sum(dp0$daf > 0.05 & dp0$daf < 0.95)
    return(res)
  },
  fix.utr5.triplet, fix.si.strict.triplet,
  poly.utr5.triplet, poly.intron.strict.triplet, SIMPLIFY = FALSE)

names(alpha.all.triplets2) <- TRIPLETS
alpha.all.triplets.df2 <- rbindlist(alpha.all.triplets2, idcol = 'triplet')
alpha.all.triplets.df2[, ratio.dp := d / p]
## MKtest for AUG gain in different gene set:compensatory uORF ####
df<-fix.utr5.triplet[[1]]
df0<-fix.si.strict.triplet[[1]]
dp<-poly.utr5.triplet[[1]]
dp0<-poly.intron.strict.triplet[[1]]
merge(dp,gene_tr,by.x="feature",by.y="tr")->dp

## total mel gain:
res <- RunAsymtoticMK2(
  df[gene%in%mel_gain$gene], df0, dp[gene%in%mel_gain$gene], dp0,
  xlow = 0.05, xhigh = 0.95, nbin = 50, output = 'table', resfmt = 'alpha')
#model         a          b  c alpha_asymptotic    CI_low   CI_high alpha_original   d  d0   p
#1 linear 0.4971024 0.09099365 NA         0.588096 0.4792017 0.6969904      0.5754426 712 334 565
#p0 alpha_raw
#1 628 0.5779584

## mel_gain & anc_loss
res1 <- RunAsymtoticMK2(
  df[gene%in%ancLoss_melGain$gene], df0, dp[gene%in%ancLoss_melGain$gene], dp0,
  xlow = 0.05, xhigh = 0.95, nbin = 10, output = 'table', resfmt = 'alpha')
#model         a         b  c alpha_asymptotic     CI_low   CI_high alpha_original
#1 linear 0.4762567 -0.301705 NA        0.1745516 -0.3810899 0.7301932      0.3949275
#d  d0  p  p0 alpha_raw
#1 46 334 49 628 0.4334672

## mel_gain without anc_loss
res2 <- RunAsymtoticMK2(
  df[gene%in%mel_gain_1522$gene], df0, dp[gene%in%mel_gain_1522$gene], dp0,
  xlow = 0.05, xhigh = 0.95, nbin = 50, output = 'table', resfmt = 'alpha')
#model         a         b  c alpha_asymptotic    CI_low  CI_high alpha_original   d  d0   p  p0
#1 linear 0.4971063 0.1109354 NA        0.6080417 0.4903495 0.725734      0.5825781 666 334 516 628
#alpha_raw
#1 0.5879383

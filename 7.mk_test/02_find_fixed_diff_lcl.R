#setwd("~/syq/project/uorf_MK/asymptoticMK_wtw/")
setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/9.compensation/MK_test_1356/")
library(data.table)

# short intron triplets ====================================================================================
# strict
fix.si.strict <- fread('fixed_diff_intron_strict_neutral.txt', sep = '\t', quote = '')

fix.si.strict.aug <- fix.si.strict[toupper(dm6_seq) == 'ATG'][
  dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0] # 如果diff_in_3sp！=0，则无法判断突变是否是在 dmel lineage 上发生的

fix.si.strict.augloss <- fix.si.strict[
  dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0]
fix.si.strict.augloss[, yak_seq := gsub('#|-', '', tstrsplit(alignment, ';')[[3]])]
fix.si.strict.augloss <- fix.si.strict.augloss[dm6_seq != 'ATG' & sim_seq == 'ATG' & yak_seq == 'ATG']

fwrite(fix.si.strict.aug, 'fixed_aug_intron_strict_neutral.txt', sep = '\t')
fwrite(fix.si.strict.augloss, 'fixed_augloss_intron_strict_neutral.txt', sep = '\t')

# looser
fix.si.looser <- fread('fixed_diff_intron_looser_neutral.txt', sep = '\t', quote = '')

fix.si.looser.aug <- fix.si.looser[toupper(dm6_seq) == 'ATG'][
  dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0]

fwrite(fix.si.looser.aug, 'fixed_aug_intron_looser_neutral.txt', sep = '\t')

# short intron sites ==============================================================================
## strict
fix.site.si.strict <- fread('fixed_site_intron_strict_neutral.txt')

# remove diff_in_3sp, repeat, gap and non dm6 lineage difference
fix.site.si.strict <- fix.site.si.strict[diff_in_3sp == 0 & flag_repeat == 0 & flag_gap == 0 & mut_tag != 'O'] # 如果diff_in_3sp！=0 (三个物种序列都不一样)，则无法判断突变是否是在 dmel lineage 上发生的
fix.site.si.strict[, site_id := paste(chrm, posn, sep = '_')]
fix.site.si.strict <- fix.site.si.strict[!duplicated(site_id)]

# divide into AUG and non-AUG differences
raw.si.strict.aug <- fix.si.strict[toupper(dm6_seq) == 'ATG']
# syq: associated the mutation site with ATG triplet
x <- fix.site.si.strict$site_id %in% raw.si.strict.aug[, paste(c(chrm, chrm, chrm), c(posn1, posn2, posn3), sep = '_')]
fix.site.si.strict.aug <- fix.site.si.strict[x]
fix.site.si.strict.nonaug <- fix.site.si.strict[!x]

fwrite(fix.site.si.strict, 'fixed_site_intron_strict_neutral_clean.txt', sep = '\t')
fwrite(fix.site.si.strict.nonaug, 'fixed_site_intron_strict_neutral_clean_nonaug.txt', sep = '\t')

## looser
fix.site.si.looser <- fread('fixed_site_intron_looser_neutral.txt')

# remove diff_in_3sp, repeat, gap and non dm6 lineage difference
fix.site.si.looser <- fix.site.si.looser[diff_in_3sp == 0 & flag_repeat == 0 & flag_gap == 0 & mut_tag != 'O']
fix.site.si.looser[, site_id := paste(chrm, posn, sep = '_')]
fix.site.si.looser <- fix.site.si.looser[!duplicated(site_id)]

# divide into AUG and non-AUG differences
raw.si.looser.aug <- fix.si.looser[toupper(dm6_seq) == 'ATG']
x <- fix.site.si.looser$site_id %in% raw.si.looser.aug[, paste(c(chrm, chrm, chrm), c(posn1, posn2, posn3), sep = '_')]
fix.site.si.looser.aug <- fix.site.si.looser[x]
fix.site.si.looser.nonaug <- fix.site.si.looser[!x]

fwrite(fix.site.si.looser, 'fixed_site_intron_looser_neutral_clean.txt', sep = '\t')
fwrite(fix.site.si.looser.nonaug, 'fixed_site_intron_looser_neutral_clean_nonaug.txt', sep = '\t')

# ## lazy way of get all the fixed sites
fix.si.strict.all <- fix.si.strict[dm6_mut > 0, .(chrm, posn1, posn2, posn3, mut_tag)]

fix.si.strict.all[, paste0('tag', 1:3) := tstrsplit(mut_tag, '')]
fix.si.strict.all[, mut_tag := NULL]

fix.si.strict.all <- melt(fix.si.strict.all, id.vars = 'chrm', measure.vars = patterns('^posn', '^tag'))
fix.si.strict.all[, variable := NULL]

setnames(fix.si.strict.all, c('chrm', 'posn', 'mutype'))
fix.si.strict.all <- fix.si.strict.all[mutype != 'O']
fix.si.strict.all[, chrmpos := paste(chrm, posn, sep = '_')]

fix.si.strict.all <- fix.si.strict.all[!duplicated(chrmpos)]
# fwrite(fix.si.strict.all, 'fixed_diff_intron_strict_neutral.txt.allsite', sep = '\t')

# 5' UTR ============================================================================================
## triplet
fix.utr5 <- fread('fixed_diff_utr5.txt')
fix.utr5[, id := paste(chrm, posn1, posn2, posn3, sep = '_')]

# make a site bed
fix.utr5.site.bed <- melt(
  fix.utr5, id.vars = c('chrm', 'id', 'strand'),
  measure.vars = c('posn1', 'posn2', 'posn3'))
fix.utr5.site.bed <- fix.utr5.site.bed[, .(
  chrm = sub('^chr', '', chrm), start = value - 1, end = value,
  name = id, score = substr(variable, 5, 5), strand)] #syq: there a lot duplication 
fwrite(fix.utr5.site.bed, 'fixed_diff_utr5_site.bed', col.names = FALSE, sep = '\t')

# find those overlapping with CDS
# syq: bedtools intersect -a fixed_diff_utr5_site.bed -b /gpfs2/zhangh/backup/database/flybase/dmel_r6.04_FB2015_01/dmel-all-r6.04.gtf.CDS.bed >fixed_diff_utr5_site_cds.ov.bed 
fix.utr5.cds.ov <- fread("./fixed_diff_utr5_site_cds.ov.bed")

# remove triplets containing CDS overlapping sites
fix.utr5.nocds <- fix.utr5[!id %in% fix.utr5.cds.ov$V4]
fwrite(fix.utr5.nocds, file = 'fixed_diff_utr5_nocds.txt', sep = '\t')

# extract fixed AUGs
fix.utr5.aug <- fix.utr5.nocds[toupper(dm6_seq) == 'ATG'][
  dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0][
    !duplicated(id)]
fwrite(fix.utr5.aug, 'fixed_aug_utr5.txt', sep = '\t')

fix.utr5.augloss <- fix.utr5.nocds[
  dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0]
fix.utr5.augloss[, yak_seq := gsub('#|-', '', tstrsplit(alignment, ';')[[3]])]
fix.utr5.augloss <- fix.utr5.augloss[dm6_seq != 'ATG' & sim_seq == 'ATG' & yak_seq == 'ATG']
fix.utr5.augloss <- fix.utr5.augloss[!duplicated(id)]
fwrite(fix.utr5.augloss, 'fixed_augloss_utr5.txt')

## site
fix.site.utr5 <- fread('fixed_site_utr5.txt')

# remove diff_in_3sp, repeat, gap and non dm6 lineage difference
fix.site.utr5 <- fix.site.utr5[diff_in_3sp == 0 & flag_repeat == 0 & flag_gap == 0 & mut_tag != 'O']
fix.site.utr5[, site_id := paste(chrm, posn, sep = '_')]
fix.site.utr5 <- fix.site.utr5[!duplicated(site_id)]

# remove CDS overlapping sites
fix.site.utr5.bed <- fix.site.utr5[, .(
  chrm = sub('chr', '', chrm), start = posn - 1, end = posn, name = site_id, score = 0, strand
)]
fwrite(fix.site.utr5.bed, 'Rtmp_fix_site_utr5.bed', col.names = FALSE, sep = '\t')


# bedtools intersect -a Rtmp_fix_site_utr5.bed -b /gpfs2/zhangh/backup/database/flybase/dmel_r6.04_FB2015_01/dmel-all-r6.04.gtf.CDS.bed  > fix.site.utr5.cds.ov.bed
fix.site.utr5.cds.ov <- fread("./fix.site.utr5.cds.ov.bed")

fix.site.utr5.final <- fix.site.utr5[!site_id %in% fix.site.utr5.cds.ov$V4]

# get non-AUG sites
raw.utr5.aug <- fix.utr5[toupper(dm6_seq) == 'ATG']
x <- fix.site.utr5.final$site_id %in% raw.utr5.aug[, paste(c(chrm, chrm, chrm), c(posn1, posn2, posn3), sep = '_')]

fix.site.utr5.final.nonaug <- fix.site.utr5.final[!x]

fwrite(fix.site.utr5.final, 'fixed_site_utr5_clean.txt', sep = '\t')
fwrite(fix.site.utr5.final.nonaug, 'fixed_site_utr5_clean_nonaug.txt', sep = '\t')

# 3' UTR ============================================================================================
## triplet
fix.utr3 <- fread('fixed_diff_utr3.txt')
fix.utr3[, id := paste(chrm, posn1, posn2, posn3, sep = '_')]

# make a site bed
fix.utr3.site.bed <- melt(
  fix.utr3, id.vars = c('chrm', 'id', 'strand'),
  measure.vars = c('posn1', 'posn2', 'posn3'))
fix.utr3.site.bed <- fix.utr3.site.bed[, .(
  chrm = sub('^chr', '', chrm), start = value - 1, end = value,
  name = id, score = substr(variable, 5, 5), strand)]
fwrite(fix.utr3.site.bed, 'fixed_diff_utr3_site.bed', col.names = FALSE, sep = '\t')

# find those overlapping with CDS
# bedtools intersect -a fixed_diff_utr3_site.bed -b /gpfs2/zhangh/backup/database/flybase/dmel_r6.04_FB2015_01/dmel-all-r6.04.gtf.CDS.bed  > fix.site.utr3.cds.ov.bed
fix.utr3.cds.ov <- fread("./fix.site.utr3.cds.ov.bed")

# remove triplets containing CDS overlapping sites
fix.utr3.nocds <- fix.utr3[!id %in% fix.utr3.cds.ov$V4]
fwrite(fix.utr3.nocds, file = 'fixed_diff_utr3_nocds.txt', sep = '\t')

# extract fixed AUGs
fix.utr3.aug <- fix.utr3.nocds[toupper(dm6_seq) == 'ATG'][
  dm6_mut > 0 & diff_in_3sp == 0 & repeat_break == 0 & ngaps == 0 & other_indel == 0][
    !duplicated(id)]
fwrite(fix.utr3.aug, 'fixed_aug_utr3.txt', sep = '\t')

## site
fix.site.utr3 <- fread('fixed_site_utr3.txt')

# remove diff_in_3sp, repeat, gap and non dm6 lineage difference
fix.site.utr3 <- fix.site.utr3[diff_in_3sp == 0 & flag_repeat == 0 & flag_gap == 0 & mut_tag != 'O']
fix.site.utr3[, site_id := paste(chrm, posn, sep = '_')]
fix.site.utr3 <- fix.site.utr3[!duplicated(site_id)]

# remove CDS overlapping sites
fix.site.utr3.bed <- fix.site.utr3[, .(
  chrm = sub('chr', '', chrm), start = posn - 1, end = posn, name = site_id, score = 0, strand
)]
fwrite(fix.site.utr3.bed, 'Rtmp_fix_site_utr3.bed', col.names = FALSE, sep = '\t')

# bedtools intersect -a Rtmp_fix_site_utr3.bed -b /gpfs2/zhangh/backup/database/flybase/dmel_r6.04_FB2015_01/dmel-all-r6.04.gtf.CDS.bed  > fix.site.utr3.cds.ov.bed
fix.site.utr3.cds.ov <- fread("./fix.site.utr3.cds.ov.bed")


fix.site.utr3.final <- fix.site.utr3[!site_id %in% fix.site.utr3.cds.ov$V4]

# get non-AUG sites
raw.utr3.aug <- fix.utr3[toupper(dm6_seq) == 'ATG']
x <- fix.site.utr3.final$site_id %in% raw.utr3.aug[, paste(c(chrm, chrm, chrm), c(posn1, posn2, posn3), sep = '_')]

fix.site.utr3.final.nonaug <- fix.site.utr3.final[!x]

fwrite(fix.site.utr3.final, 'fixed_site_utr3_clean.txt', sep = '\t')
fwrite(fix.site.utr3.final.nonaug, 'fixed_site_utr3_clean_nonaug.txt', sep = '\t')

##not done: fixed sites in CDS ###########################################################################
fix.site.cds <- fread('fixed_site_cds.txt')

# remove diff_in_3sp, repeat, gap and non dm6 lineage difference
fix.site.cds <- fix.site.cds[diff_in_3sp == 0 & flag_repeat == 0 & flag_gap == 0 & mut_tag != 'O']
fix.site.cds[, site_id := paste(chrm, posn, sep = '_')]
fix.site.cds <- fix.site.cds[!duplicated(site_id)]

# make a vcf of fixed sites
# #CHROM POS     ID        REF    ALT     QUAL FILTER INFO
fix.site.cds.vcf <- fix.site.cds[, .(
  '#CHROM' = chrm, POS = posn, ID = site_id, REF = sim_seq,
  ALT = dm6_seq, QUAL = '100.0', FILTER = 'PASS', INFO = 'AF=1.00'
)]
fwrite(fix.site.cds.vcf, 'fixed_site_cds_clean.vcf', sep = '\t')
# cat fixed_site_cds_clean.vcf | java -jar ~/local/src/snpEff/snpEff.jar -nodownload -canon -onlyProtein -no-downstream -no-upstream -no-intron -no-utr -no-intergenic -v BDGP6.82 | java -jar ~/local/src/snpEff/SnpSift.jar extractFields -s "," -e "." -  CHROM POS ID REF ALT "ANN[*].EFFECT" "ANN[*].HGVS_P" "ANN[*].GENEID" 
# >fixed_site_cds_clean.snpeff
#lcl: cat fixed_site_cds_clean.vcf | java -jar ~/lcl/software/snpEff/snpEff.jar -nodownload -canon -onlyProtein -no-downstream -no-upstream -no-intron -no-utr -no-intergenic -v BDGP6.28.99 | java -jar ~/lcl/software/snpEff/SnpSift.jar extractFields -s "," -e "." -  CHROM POS ID REF ALT "ANN[*].EFFECT" "ANN[*].HGVS_P" "ANN[*].GENEID"   >fixed_site_cds_clean.snpeff
fix.site.cds.snpeff <- fread('fixed_site_cds_clean.snpeff')
setnames(fix.site.cds.snpeff, c('#CHROM', 'ANN[*].EFFECT', 'ANN[*].HGVS_P', 'ANN[*].GENEID'),
         c('CHROM', 'EFFECT', 'HGVS_P', 'GENEID'))
fix.site.cds.snpeff[, type := ifelse(
  grepl('missense_variant', EFFECT), 'Nonsyn',
  ifelse(grepl('synonymous_variant', EFFECT), 'Syn', 'Other')
)]
fix.site.cds.snpeff <- fix.site.cds.snpeff[type != 'Other']

fix.site.cds.final <- merge(fix.site.cds, fix.site.cds.snpeff[, .(site_id = ID, type)], by = 'site_id')
fwrite(fix.site.cds.final, 'fixed_site_cds_clean.txt', sep = '\t')

## fixed site in CDS containing repeats (for CCWU) ###########################################################
if(FALSE) local({
  fix.site.cds <- fread('fixed_site_cds.txt')
  
  # remove diff_in_3sp, repeat, gap and non dm6 lineage difference
  fix.site.cds <- fix.site.cds[diff_in_3sp == 0 & flag_gap == 0 & mut_tag != 'O']
  fix.site.cds[, site_id := paste(chrm, posn, sep = '_')]
  fix.site.cds <- fix.site.cds[!duplicated(site_id)]
  
  # make a vcf of fixed sites
  # #CHROM POS     ID        REF    ALT     QUAL FILTER INFO
  fix.site.cds.vcf <- fix.site.cds[, .(
    '#CHROM' = chrm, POS = posn, ID = site_id, REF = sim_seq,
    ALT = dm6_seq, QUAL = '100.0', FILTER = 'PASS', INFO = 'AF=1.00'
  )]
  fwrite(fix.site.cds.vcf, 'fixed_site_cds_wrepeat.vcf', sep = '\t')
  
  fix.site.cds.snpeff <- fread('fixed_site_cds_wrepeat.snpeff')
  setnames(fix.site.cds.snpeff, c('#CHROM', 'ANN[*].EFFECT', 'ANN[*].HGVS_P', 'ANN[*].GENEID'),
           c('CHROM', 'EFFECT', 'HGVS_P', 'GENEID'))
  fix.site.cds.snpeff[, type := ifelse(
    grepl('missense_variant', EFFECT), 'Nonsyn',
    ifelse(grepl('synonymous_variant', EFFECT), 'Syn', 'Other')
  )]
  fix.site.cds.snpeff <- fix.site.cds.snpeff[type != 'Other']
  
  fix.site.cds.final <- merge(fix.site.cds, fix.site.cds.snpeff[, .(site_id = ID, type)], by = 'site_id')
  fwrite(fix.site.cds.final, 'fixed_site_cds_wrepeat.txt', sep = '\t')
})
## get conserved 5' UTR AUGs ##########################################################################################

#zcat r6.04_utr5_12col.maf.gz | python3 find_conserved_aug_syq.py utr5 >conserved_aug_utr5_syq.txt
conserved.aug.utr5 <- fread('conserved_aug_utr5.txt')
conserved.aug.utr5[, gid := paste(chrm, posn1, sep = '_')]

dmel.cds <- fread('./dmel-all-r6.04.gtf.CDS.bed')
dmel.cds[, V1 := paste0('chr', V1)]

tmp1 <- dmel.cds[conserved.aug.utr5, on = .(V1 == chrm, V2 < posn1, V3 >= posn1), nomatch = 0]
tmp2 <- dmel.cds[conserved.aug.utr5, on = .(V1 == chrm, V2 < posn3, V3 >= posn3), nomatch = 0]

conserved.aug.utr5.clean <- conserved.aug.utr5[!gid %in% tmp1$gid & !gid %in% tmp2$gid]
conserved.aug.utr5.clean <- conserved.aug.utr5.clean[dm6_mut == 0 & diff_in_3sp == 0 & repeat_break == 0]

conserved.aug.utr5.clean[, yak_seq := tstrsplit(gsub('#|-', '', alignment), ';', keep = 3)]
conserved.aug.utr5.clean <- conserved.aug.utr5.clean[!(ngaps > 0 & sim_seq != 'ATG' & yak_seq != 'ATG')]
# syq: conserved.aug.utr5.clean3 <- conserved.aug.utr5.clean[ ngaps == 0 | sim_seq == 'ATG' | yak_seq == 'ATG'] 和上面等价

fwrite(conserved.aug.utr5.clean, file = 'conserved_aug_utr5_clean.txt', sep = '\t')


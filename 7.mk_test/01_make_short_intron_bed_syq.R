#setwd("~/syq/project/uorf_MK/asymptoticMK_wtw/")
setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/9.compensation/MK_test/")
library(data.table)

dmel.gtf <- fread('/lustre/user/lulab/luoshq/syq/reference/dmel-all-r6.04.gtf')
dmel.gtf <- tidyr::extract(
  dmel.gtf, col = 'V9', into = c('gene', 'isoform'),
  regex = 'gene_id "(.*?)".*?transcript_id "(.*?)"')

setDT(dmel.gtf)
dmel.gtf.exon <- dmel.gtf[V3 == 'exon']
dmel.gtf.exon[, total.exons := .N, by = 'isoform']

# get 1-based closed interval of introns
dmel.gtf.intron <- dmel.gtf.exon[total.exons > 1][order(isoform, V4)][
  , .(chrm = V1[1], start = V5[-.N] + 1, end = V4[-1] - 1, strand = V7[1], gene = gene[1]), by = 'isoform']
dmel.gtf.intron[, intron_name := paste(isoform, seq_len(.N), sep = '_'), by = 'isoform']

# output in bed format
dmel.gtf.intron.bed <- dmel.gtf.intron[, .(
  chrm, start = start - 1, end, intron_name, score = 0, strand)][order(chrm, start)]
fwrite(dmel.gtf.intron.bed, file = 'dmel_intron_r6.04.bed', col.names = FALSE, sep = '\t')

# remove intron from all exons (both use r6.04)
# bedtools subtract -a dmel_intron_r6.04.bed -b /gpfs2/zhangh/backup/database/flybase/dmel_r6.04_FB2015_01/dmel-all-r6.04.gtf.exon.bed >dmel_intron_exon_substracted.r6.04.bed
# bedtools sort -i dmel_intron_exon_substracted.r6.04.bed | bedtools merge -i stdin -s -c 6 -o distinct >dmel_intron_exon_substracted.merge.r6.04.bed

# readIN subtracted-exon introns and compare the introns before and after exon subtraction
dmel.intron.exon.subtracted <- fread('dmel_intron_exon_substracted.r6.04.bed')
dmel.intron.exon.subtracted[, score := V3 - V2]
dmel.intron.only.len <- dmel.intron.exon.subtracted[, .(score = sum(score)), by = 'V4']

short.intron.bed <- dmel.gtf.intron.bed[end - start <= 65]
short.intron.bed[, score := end - start]

short.intron.bed[dmel.intron.only.len, intron.only := i.score, on = .(intron_name = V4)]
short.intron.uniq <- short.intron.bed[!duplicated(paste(chrm, start, end))]
short.intron.uniq[, .(diff = score - intron.only)][,.N, by = diff][order(-N)][1:5]
short.intron.uniq[, ov.exon := is.na(intron.only) | intron.only < score]
#    diff-base     N
# 1:    0 18487 / 20158 == 92%
# 2:   NA   962 / 20158 == 4.8%
# 3:    3   171
# 4:    4    66
# 5:    6    60

# how many short intron are completely included by other intron intervals?
# bedtools sort -i dmel_intron_r6.04.bed | bedtools merge -i stdin -s -c 6 -o distinct >dmel_intron.merge.r6.04.bed
dmel.intron.merge <- fread('dmel_intron.merge.r6.04.bed')
dmel.intron.merge[, handle := paste(V1, V2, V3, V4, sep = '|')]
short.intron.uniq[, handle := paste(chrm, start, end, strand, sep = '|')]

short.intron.uniq[, consistent := handle %in% dmel.intron.merge$handle]

short.intron.uniq[, table(ov.exon, consistent)]
# ov.exon FALSE  TRUE
#   FALSE  1318 17169 (85%)
#   TRUE    348  1323

# strict definations of short introns
# consistent between different isoforms and do not overlap with exons #syq: actually do not overlap with exons and do not included in other longer introns(redundancy removal)
short.intron.strict <- short.intron.uniq[ov.exon == FALSE & consistent]
short.intron.strict.neutral <- short.intron.strict[
  , .(chrm,
      start = ifelse(strand == '+', start + 7,  end - 30),
      end   = ifelse(strand == '+', start + 30, end - 7),
      intron_name, score = 0, strand
  )]
short.intron.strict.neutral[, intron_name := paste(chrm, start, end, sep = '_')]

# looser definations of short introns
# 8-30nt of all short introns, subtract exonic region #syq : just substract the overlapped region, not drop the whole intron 
short.intron.looser <- copy(short.intron.uniq)

short.intron.looser.neutral <- short.intron.looser[
  , .(chrm,
      start = ifelse(strand == '+', start + 7,  end - 30),
      end   = ifelse(strand == '+', start + 30, end - 7),
      intron_name, score = 0, strand
  )][order(chrm, start)]
fwrite(short.intron.looser.neutral, 'Rtmp_silnbf.bed', col.names = FALSE, sep = '\t')

#syq: bedtools subtract -a Rtmp_silnbf.bed -b /gpfs2/zhangh/backup/database/flybase/dmel_r6.04_FB2015_01/dmel-all-r6.04.gtf.exon.bed|bedtools merge -i stdin -s -c 6 -o distinct > short.intron.looser.neutral.tmp
short.intron.looser.neutral.final <- fread("short.intron.looser.neutral.tmp")
short.intron.looser.neutral.final <- short.intron.looser.neutral.final[
  , .(chrm = V1, start = V2, end = V3, intron_name = paste(V1, V2, V3, sep = '_'), score = 0, strand = V4)]
file.remove('Rtmp_silnbf.bed','short.intron.looser.neutral.tmp')

## save results
fwrite(short.intron.strict.neutral, file = 'short_intron_strict_neutral.bed', col.names = FALSE, sep = '\t')
fwrite(short.intron.looser.neutral.final, file = 'short_intron_looser_neutral.bed', col.names = FALSE, sep = '\t')

# with chr prefix
short.intron.strict.neutral.chr <- copy(short.intron.strict.neutral)
short.intron.strict.neutral.chr[, `:=`(chrm = paste0('chr', chrm))]

short.intron.looser.neutral.final.chr <- copy(short.intron.looser.neutral.final)
short.intron.looser.neutral.final.chr[, `:=`(chrm = paste0('chr', chrm))]

fwrite(short.intron.strict.neutral.chr, file = 'short_intron_strict_neutral.chr.bed', col.names = FALSE, sep = '\t')
fwrite(short.intron.looser.neutral.final.chr, file = 'short_intron_looser_neutral.chr.bed', col.names = FALSE, sep = '\t')

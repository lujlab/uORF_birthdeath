#setwd("/lustre/user/lulab/luoshq/syq/project/uorf_MK/asymptoticMK_wtw")
setwd("/gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/9.compensation/MK_test/")
#### Get sequence of short intron neutral regions:####
# bedtools getfasta -fi /lustre/user/lulab/luoshq/syq/reference/dmel-all-chromosome-r6.04.fasta -bed short_intron_strict_neutral.bed -name -s -fo short_intron_strict_neutral.fa
# bedtools getfasta -fi /lustre/user/lulab/luoshq/syq/reference/dmel-all-chromosome-r6.04.fasta -bed short_intron_looser_neutral.bed -name -s -fo short_intron_looser_neutral.fa

#### Get alignment of 8-30nt of short introns ####
# pip install bx-python
# modified ZH's MAFtools3.py: maketrans -> str.maketrans https://www.runoob.com/python3/python3-string-maketrans.html
### python deprecated_ZH/MAFtools3_lcl.py bed /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/7.new_27_maf dm6 short_intron_strict_neutral.chr.bed > short_intron_strict_neutral.maf
# python deprecated_ZH/MAFtools3_syq.py bed /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/7.new_27_maf/chr dm6 short_intron_strict_neutral.chr.bed > short_intron_strict_neutral.maf
# python deprecated_ZH/MAFtools3_syq.py bed /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/7.new_27_maf/chr dm6 short_intron_looser_neutral.chr.bed > short_intron_looser_neutral.maf


# gzip short_intron_strict_neutral.maf
# gzip short_intron_looser_neutral.maf

#### Extract newly fixed differences in different regions ####
# newly fixed triplets
# zcat short_intron_strict_neutral.maf.gz | python3 find_newly_fixed_triplet_lcl.py intron > fixed_diff_intron_strict_neutral.txt
# zcat short_intron_looser_neutral.maf.gz | python3 find_newly_fixed_triplet_lcl.py intron > fixed_diff_intron_looser_neutral.txt

# awk '{print "chr"$0}' dmel-all-r6.04.gtf.3UTR.bed > dmel-all-r6.04.gtf.3UTR.chr.bed
# awk '{print "chr"$0}' dmel-all-r6.04.gtf.5UTR.bed > dmel-all-r6.04.gtf.5UTR.chr.bed
# awk '{print "chr"$0}' dmel-all-r6.04.gtf.CDS.bed > dmel-all-r6.04.gtf.CDS.chr.bed
# python deprecated_ZH/MAFtools3_lcl.py bed /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/7.new_27_maf/chr dm6 dmel-all-r6.04.gtf.3UTR.chr.bed > r6.04_utr3_12col.maf
# python deprecated_ZH/MAFtools3_lcl.py bed /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/7.new_27_maf/chr dm6 dmel-all-r6.04.gtf.5UTR.chr.bed > r6.04_utr5_12col.maf
# python deprecated_ZH/MAFtools3_lcl.py bed /gpfs2/liucl/lcl/analysis/uORF_gain_loss/species_27_maf/dm6_other_sim_genome/7.new_27_maf/chr dm6 dmel-all-r6.04.gtf.CDS.chr.bed > r6.04_cds_12col.maf
# gzip r6.04_utr5_12col.maf
# gzip r6.04_utr3_12col.maf
# gzip r6.04_cds_12col.maf
# zcat r6.04_utr5_12col.maf.gz | python3 find_newly_fixed_triplet_lcl.py utr5 >fixed_diff_utr5.txt
# zcat r6.04_utr3_12col.maf.gz | python3 find_newly_fixed_triplet_lcl.py utr3 >fixed_diff_utr3.txt


# newly fixed sites
# zcat short_intron_strict_neutral.maf.gz | python3 find_newly_fixed_site_lcl.py intron >fixed_site_intron_strict_neutral.txt
# zcat short_intron_looser_neutral.maf.gz | python3 find_newly_fixed_site_lcl.py intron >fixed_site_intron_looser_neutral.txt
# zcat r6.04_utr5_12col.maf.gz | python3 find_newly_fixed_site_lcl.py utr5 >fixed_site_utr5.txt
# zcat r6.04_utr3_12col.maf.gz | python3 find_newly_fixed_site_lcl.py utr3 >fixed_site_utr3.txt


# Refine fixed differences in R: `02_find_fixed_diff.R`


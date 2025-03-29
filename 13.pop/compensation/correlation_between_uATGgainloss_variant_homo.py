#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#by LIUchenlu 2024/11/19
#usage correlation_between_uATGgainloss_variant_homo.py  gain.vcf loss.vcf samplelist_file out.table

#pip install pyvcf
import sys
import vcf
import csv
import scipy.stats as stats 
uAUG_gained_vcf = sys.argv[1] 
uAUG_lost_vcf = sys.argv[2]
samplelist_file = sys.argv[3] 
outfile = sys.argv[4]

#import vcf files
vcf_gain = vcf.Reader(filename=uAUG_gained_vcf)
vcf_loss = vcf.Reader(filename=uAUG_lost_vcf)
#import sample file
def read_sample(file_path):
  with open(file_path, 'r') as file:
    reader = csv.reader(file)
    first_column = [row[0] for row in reader]
  return first_column
sample_list = read_sample(samplelist_file)


for record_gain in vcf_gain:
    allele1_chr = record_gain.CHROM
    allele1_pos = record_gain.POS
    allele1_ref = record_gain.REF
    allele1_alt = str(record_gain.ALT[0])
    vcf_loss = vcf.Reader(filename=uAUG_lost_vcf)
    for record_loss in vcf_loss:
        allele2_chr = record_loss.CHROM
        if allele1_chr != allele2_chr:
           continue
        allele2_pos = record_loss.POS
        allele2_ref = record_loss.REF
        allele2_alt = str(record_loss.ALT[0])
        type_AB = 0
        type_Ab = 0
        type_aB = 0
        type_ab = 0
        type_unknown = 0
        for s in sample_list:
            if record_gain.genotype(s)['GT'] in ["0|0","0/0"] and record_loss.genotype(s)['GT'] in ["0|0","0/0"]:
                type_AB = type_AB + 1
            elif record_gain.genotype(s)['GT'] in ["1|1","1/1"] and record_loss.genotype(s)['GT'] in ["0|0","0/0"]:
                type_aB = type_aB + 1
            elif record_gain.genotype(s)['GT'] in ["0|0","0/0"] and record_loss.genotype(s)['GT'] in ["1|1","1/1"]:
                type_Ab = type_Ab + 1
            elif record_gain.genotype(s)['GT'] in ["1|1","1/1"] and record_loss.genotype(s)['GT'] in ["1|1","1/1"]:
                type_ab = type_ab + 1
            else:
                type_unknown = type_unknown + 1
        #fisher test, odd_ratio=ab*AB/(aB*Ab)
        odd_ratio, p_value = stats.fisher_exact([[type_AB, type_Ab], [type_aB, type_ab]] )
        out=[allele1_chr,str(allele1_pos),allele1_ref,allele1_alt,allele2_chr,str(allele2_pos),allele2_ref,allele2_alt,
             str(type_AB),str(type_Ab),str(type_aB),str(type_ab),str(type_unknown),str(odd_ratio),str(p_value)]
        out2 = '\t'.join(out)
        print(out2, file=open(outfile, "a"))



# -*- coding: utf-8 -*-
"""
filter mRNA reads for calculating library size;
-------------------------------------------------------------------------------
@author: zh (mt1022)
@date: Sun Feb 14 16:20:16 2016
env: python3.4 or higher
"""
import pysam

sam_stream = pysam.Samfile('-', 'r')
min_len = 18
max_len = 1000

lib_size = 0
for seg in sam_stream:
    i = seg.query_alignment_length
    if seg.mapping_quality >= 10 and min_len <= i <= max_len:
        lib_size = lib_size + 1

print(lib_size)

"""
=================================
author: mt1022 (zh)
date: 2017-02-20 09:59
"""
import sys


for line in open(sys.argv[1], 'r'):
    ary = line.strip().split('\t')
    feature_name = ary[3] + '_' + ary[4]
    start = int(ary[1])
    end = int(ary[2])
    for i in range(end - start):
        out = [ary[0], str(start + i), str(start + i + 1), feature_name, str(i), ary[5]]
        print('\t'.join(out))


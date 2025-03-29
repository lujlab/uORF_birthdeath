"""
make meta-gene profile files for CDS and uORF
===================================================
author: mt1022 (zh)
date: 2017-02-21 19:50
"""
import sys
import argparse
import statistics as stat


def meta_profile(cov_path, out_prefix):
    """ """
    # prepare
    cds_path = open('cds_interval.txt', 'r')
    uorf_path = open('uorf_interval.txt', 'r')
    out_cds = open(out_prefix + '_cds.txt', 'w')
    out_uorf = open(out_prefix + '_uorf.txt', 'w')

    # get cds start and end on transcripts
    cds_interval = dict()
    for line in cds_path:
        ary = line.rstrip().split('\t')
        cds_interval[ary[0]] = tuple(int(i) for i in ary[1:])

    # get transcript coverage
    trans_cov = dict()
    for line in cov_path:
        ary = line.rstrip().split('\t')
        if ary[3] in cds_interval:
            trans_cov[ary[3]] = [float(i) for i in ary[6].split(',')]

    # get cds mean and cds median
    cds_stat = dict()
    for i in cds_interval:
        i_start, i_end = cds_interval[i]
        i_cov = trans_cov[i]
        i_codon = [sum(i_cov[j:(j + 3)]) for j in range(i_start - 1, i_end, 3)]
        cds_stat[i] = [stat.mean(i_codon), stat.median(i_codon)]

    # make profile of each CDS
    for i in cds_interval:
        i_start, i_end = cds_interval[i]
        i_left = i_start - 1 - 10 * 3
        i_right = i_start + 2 + 40 * 3
        i_cov = trans_cov[i]
        i_meta = []
        for j in range(i_left, i_right, 3):
            if j < 0:
                i_meta.append('NA')
            elif j + 3 > len(i_cov):
                i_meta.append('NA')
            else:
                i_meta.append('%.5f' % sum(i_cov[j:(j + 3)]))
        out = [i] + ['%.5f' % j for j in cds_stat[i]] + i_meta
        print('\t'.join(out), file=out_cds)
    out_cds.close()

    # make profile of each uORF
    for line in uorf_path:
        name, i, i_start = line.strip().split('\t')
        i_start = int(i_start)
        i_left = i_start - 1 - 5 * 3
        i_right = i_start + 2 + 10 * 3
        i_cov = trans_cov[i]
        i_meta = []
        for j in range(i_left, i_right, 3):
            if j < 0:
                i_meta.append('NA')
            elif j + 3 > len(i_cov):
                i_meta.append('NA')
            else:
                i_meta.append('%.5f' % sum(i_cov[j:(j + 3)]))
        out = [name] + ['%.5f' % j for j in cds_stat[i]] + i_meta
        print('\t'.join(out), file=out_uorf)
    out_uorf.close()
    return

if __name__ == '__main__':
    args_cov = open(sys.argv[1], 'r')
    args_prefix = sys.argv[2]
    meta_profile(args_cov, args_prefix)


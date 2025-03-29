"""
sum up reads mapped to a transcript interval based on transcript perbase cov;
=================================
author: mt1022 (zh)
date: 2016-10-04 19:51
"""
import sys


def sum_up_reads(pnt_cov_path=None, intv_path=None):
    assert pnt_cov_path is not None
    assert intv_path is not None
    if intv_path == '-':
        intv = sys.stdin
    else:
        intv = open(intv_path, 'r')

    pnt_cov = dict()
    for line in open(pnt_cov_path, 'r'):
        ary = line.rstrip().split()
        pnt_cov[ary[3]] = [float(i) for i in ary[6].split(',')]

    for line in intv:
        ary = line.rstrip().split()
        name = '{0}_{1}'.format(ary[0], ary[1])
        try:
            reads = sum(pnt_cov[ary[0]][(int(ary[1]) - 1):int(ary[2])])
            reads = str(reads)
        except KeyError:
            reads = 'NA'
        print(name + '\t' + reads)
    return


if __name__ == '__main__':
    sum_up_reads(*sys.argv[1:])

"""
group coverage
=================================
author: mt1022 (zh)
date: 2017-02-20 11:09
"""
import sys
import argparse
import itertools


def group_pntcov(in_stream, out_stream):
    """ """
    for k, g in itertools.groupby(in_stream, lambda line: line.split('\t')[3]):
        lines = [line.strip().split('\t') for line in g]
        start = lines[0][1]
        end = lines[-1][2]
        cov = ','.join([i[6] for i in lines])
        out = [lines[0][0], start, end, k, '.', lines[0][5], cov]
        print('\t'.join(out), file=out_stream)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='group per nt coverage')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--output', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

    args = parser.parse_args()
    group_pntcov(args.input, args.output)


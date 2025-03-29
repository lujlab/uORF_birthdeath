"""
filter reads and scale weight by read length
================================================
author: mt1022 (zh)
date: 2017-02-19 20:22
"""
import sys
import argparse


def scale_readlen(in_stream, out_stream):
    """ """
    for line in in_stream:
        ary = line.rstrip().split('\t')
        if int(ary[4]) < 10:
            continue
        rlen = sum([int(i) for i in ary[10].split(',')])
        if rlen < 3 or rlen > 10:
            continue
        ary[4] = '%.5f' % (1 / rlen)
        print('\t'.join(ary), file=out_stream)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='filter reads / scale read length')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--output', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

    args = parser.parse_args()
    scale_readlen(args.input, args.output)


"""
filter reads and scale weight by read length
    default params are for ribo-seq of 27-34nt after removing 12nt from each end;
================================================
author: mt1022 (zh)
date: 2017-02-19 20:22
"""
import sys
import argparse


def scale_readlen(in_stream, out_stream, min_len, max_len):
    """ """
    for line in in_stream:
        ary = line.rstrip().split('\t')
        if int(ary[4]) < 10:
            continue
        rlen = sum([int(i) for i in ary[10].split(',')])
        if rlen < min_len or rlen > max_len:
            continue
        ary[4] = '%.5f' % (1 / rlen)
        print('\t'.join(ary), file=out_stream)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='filter reads / scale read length')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--output', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-m', '--min_len', nargs='?', type=int, default=3)
    parser.add_argument('-M', '--max_len', nargs='?', type=int, default=10)

    args = parser.parse_args()
    scale_readlen(args.input, args.output, args.min_len, args.max_len)


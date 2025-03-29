"""
concatenate perbase coverage of all exons to transcripts;
Notice!!!
the input all exons from the same transcript should be grouped together and be ordered on 5th column;

input1: exon pnt coverage
input2: sorted reference bed file
    (if the feature is on minus strand, first exon generally is the right most exon)
=================================
author: mt1022 (zh)
date: 2017-02-20 14:28
"""
import sys


exon_cov = dict()
for line in open(sys.argv[1], 'r'):
    ary = line.rstrip().split()
    exon_cov[ary[3]] = ary[6]

pbcov = dict()
for line in open(sys.argv[2], 'r'):
    ary = line.split()
    name = ary[3]
    id = '{0}_{1}'.format(ary[3], ary[4])
    try:
        cov = exon_cov[id].split(',')
        cov = [cov[i] if cov[i] != '.' else '0' for i in range(len(cov))]
        if len(cov) != (int(ary[2]) - int(ary[1])):
            print('error!', file=sys.stderr)
        if ary[5] == '-':
            cov = cov[::-1]
    except KeyError:
        cov = ['0'] * (int(ary[2]) - int(ary[1]))

    if name in pbcov:
        pbcov[name]['left'] = min(int(ary[1]), pbcov[name]['left'])
        pbcov[name]['right'] = max(int(ary[2]), pbcov[name]['right'])
        pbcov[name]['len'] += len(cov)
        pbcov[name]['cov'] = pbcov[name]['cov'] + cov
    else:
        pbcov[name] = dict(chrm=ary[0], name=name, left=int(ary[1]), right=int(ary[2]),
                           cov=cov, len=len(cov), strand=ary[5])

for name in sorted(pbcov.keys()):
    s = [str(pbcov[name][i]) for i in ['chrm', 'left', 'right', 'name', 'len', 'strand']]
    cov = pbcov[name]['cov']
    s.append(','.join(cov))
    print('\t'.join(s))


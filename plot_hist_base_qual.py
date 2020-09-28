#!/usr/bin/env python

"""
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.

usage: python %prog mapped_to_HDV.bam ref_HDV.fa
python3.7
"""


import os,sys,glob
import pysam
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


bam=sys.argv[1]
fa=sys.argv[2]

nt_pos=666  # 0-based

seq=[]
with open(fa) as infile:
    for line in infile:
        if not '>' in line:
            seq.append(line.replace('U', 'T').strip())
seq=''.join(seq)


infile=pysam.AlignmentFile(bam, 'rb')
quals=[]
rects=[]
mut_quals=[]
mut_rects=[]
for line in infile:
    if not line.is_qcfail is True:
        if line.reference_start <= nt_pos < line.reference_end:
            poss=line.get_reference_positions()
            rseq=seq[line.reference_start:line.reference_end]
            qseq=line.query_sequence
            rpos=[]
            qpos=[]
            total=poss[0]
            qp=0
            for b,l in line.cigar:
                if b == 0 or b == 7 or b == 8:
                    for _ in range(l):
                        rpos.append(total)
                        qpos.append(qp)
                        total += 1
                        qp += 1
                elif b == 2 or b == 3:
                    for _ in range(l):
                        total += 1
                elif b == 4:
                    for _ in range(l):
                        qp += 1
            mut_read=False
            for rp,qp in zip(rpos,qpos):
                if not seq[rp] == qseq[qp]:
                    if seq[rp] == 'A' and qseq[qp] == 'G':
                        if rp == nt_pos:
                            mut_read=True
            for pos,n in zip(poss, range(len(poss))):
                if pos == nt_pos:
                    pos_target=n
                    break
            qual_target=line.query_alignment_qualities[pos_target]
            if mut_read is False:
                quals.append(qual_target)
            else:
                mut_quals.append(qual_target)


# plot hist of qual
plt.figure(figsize=(4, 3))  # (x, y)
gs=gridspec.GridSpec(2, 1)  # (y, x)
gs.update(wspace=0.05)

ax=plt.subplot(gs[0])
ax.hist(mut_quals, range=(0, max(quals + mut_quals)), bins=(max(quals + mut_quals)), color='tab:red')
ax.xaxis.set_ticks([])

ax=plt.subplot(gs[1])
ax.hist(quals, range=(0, max(quals + mut_quals)), bins=(max(quals + mut_quals)), color='tab:blue')

ax.set_xlabel('Base quality, pos=%d' % (nt_pos + 1))  # 1-based
ax.set_ylabel('Count')

plt.suptitle('histogram of base quality')
plt.savefig('plot_out_qual_hist.pdf')


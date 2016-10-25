#!/usr/bin/python
# ## Usage
# - Put this script to folder where you store your scripts (e.g.: user/bin), and make it executable `chmod +X path/to/Concatenator.CELseq,py`
# - Requires homefolder/var/CELSeq1_96BC.py, homefolder/var/CELSeq2_384BC.py to be in place
# - Put read1 and read2 fastq files (concatenated by lane) in a folder
# - The files minimum name is referred as root.
# - Test it by running in bash: `python path/to/Concatenator.CELseq.py`
# - Run it with your data: python `python path/to/Concatenator.CELseq.py path/to/FastqFiles 1` # Put 1 or 2 if you used CELseq2 primers
# - Original author: Anna Alemany, van Oudenaarden group, 25-10-2016
# - Modified  by Abel Vertesy, van Oudenaarden group, 25-10-2016
# - To be improved: allowing 1 mismatch in CBCs will yield 3% more mappability (P. Lijnzaad)

import sys, os
from itertools import izip # to iterate over two files in parallel

sys.path.append('/home/hub_oudenaarden/avertesy/var/')
from CELSeq1_96BC import bc2sample as bccelseq1     # cell barcodes for CELseq 1
from CELSeq2_384BC import bc2sample as bccelseq2   # cell barcodes for CELseq 2

try:
    fqr= sys.argv[1]
    protocol = sys.argv[2]
except:
    print "Give root to fastq files and 1 or 2 according to celseq protocol"
    sys.exit()

if protocol not in ['1', '2']:
    print "Only cellseq 1 or 2 accepted so far. Please change input"
    sys.exit()

if protocol == '1':
    bclen = 8
    umilen = 4
elif protocol == '2':
    bclen = 8
    umilen = 6

# find fastq files. It assumes you concatenated

fq1 = fqr + '_R1_cat.fastq'
fq2 = fqr + '_R2_cat.fastq'

if not os.path.isfile(fq1) or not os.path.isfile(fq2):
    print 'fastq files not found'
    sys.exit()

fout = open(fqr + '_cbc.fastq', 'w+')

NrReadsIn = NrReadsOut = badCBCs = badUMIs = 0
with open(fq1) as f1, open(fq2) as f2:
    for l1, l2 in izip(f1, f2):
        l1, l2 = l1.rstrip(), l2.rstrip()
        if l1[0] == '@' and l2[0] == '@':
            l1 = l1.rsplit(' ')
            l2 = l2.rsplit(' ')
            if l1[0] != l2[0]:
                print 'check'
                sys.exit()
            q = l1[0]

            NrReadsIn += 1
            # Read the sequences:
            s1, s2 = f1.next(), f2.next()           # read next line
            s1, s2 = s1.rstrip(), s2.rstrip()       # remove \n and spaces
            bcseq = s1[:bclen+umilen]
            if protocol == "1":
                # for celseq 1 (R1: celbc + umi + polyA)
                CBCseq = bcseq[:bclen]
                if CBCseq not in bccelseq1:        # is valid? Discard read if not.
                    badCBCs += 1
                    continue
                CBC_ID = bccelseq1[CBCseq]
                UMIseq = bcseq[bclen:]
            elif protocol == "2":
                # for celseq 2 (R1: umi + celbc + polyA)
                CBCseq = bcseq[umilen:]
                if CBCseq not in bccelseq2:        # is valid?  Discard read if not.
                    badCBCs += 1
                    continue
                CBC_ID = bccelseq2[CBCseq]
                UMIseq = bcseq[:umilen]
            if "N" in UMIseq:                      # is the UMI determined? Discard read if not.
            	badUMIs += 1
                continue
            # the next line should have a plus sign:
            p1, p2 = f1.next(), f2.next()
            p1, p2 = p1.rstrip(), p2.rstrip()
            if p1[0] != '+' and p2[0] != '+':
                print 'check'
                sys.exit()
            # phreds:
            q1, q2 = f1.next(), f2.next();
            q1, q2 = q1.rstrip(), q2.rstrip()
            NrReadsOut += 1

            # Write out read1.fastq with barcode in the name
            # print >> fout, '\n'.join([':'.join([q, CBCseq, UMIseq]), s2, p2, q2])
            print >> fout, '\n'.join([':'.join([q, CBCseq, UMIseq, CBC_ID]), s2, p2, q2]) # write out with cell index number appended

print "The number of input reads\t", NrReadsIn
print "The number of output reads\t", NrReadsOut
print "The number of incorrect cell barcodes\t", badCBCs
print "The number of UMIs with an N\t", badUMIs

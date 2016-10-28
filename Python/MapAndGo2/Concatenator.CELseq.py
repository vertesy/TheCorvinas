#!/usr/bin/python
# ## Input
# - Put this script to folder where you store your scripts (e.g.: user/bin), and make it executable `chmod +X path/to/Concatenator.CELseq.py`
# - Requires homefolder/var/CELSeq1_96BC.py, homefolder/var/CELSeq2_384BC.py to be in place
# - Put read1 and read2 fastq files (concatenated by lane) in a folder
# - Provide 3 arguments: (1) the full path to R1 fastq file; (2) 1 or 2 according to celseq protocol; (3) Hamming distance (integer) for matching cell barcodes.
# - The files minimum name is referred as root.
# - Allowing 1 mismatch in CBCs will yield 3%+ more mappability
# ## Usage
# - Test it by running in bash: `python path/to/Concatenator.CELseq.py`
# - Run it with your data: python `python path/to/Concatenator.CELseq.py path/to/FastqFiles 1` # Put 1 or 2 if you used CELseq2 primers
# - Original author: Anna Alemany, van Oudenaarden group, 25-10-2016
# - Modified  by Abel Vertesy, van Oudenaarden group, 27-10-2016

print "Concatenator.CELseq.py Started, eyy!"
import sys, os
from itertools import izip # to iterate over two files in parallel

def_bar = "/home/hub_oudenaarden/avertesy/var/"
sys.path.append(def_bar)
from CELSeq1_96BC import bc2sample as bccelseq1     # cell barcodes for CELseq 1
from CELSeq2_384BC import bc2sample as bccelseq2   # cell barcodes for CELseq 2

# Function and object definitions -------------------------------------------------------


#### Indetify the single closest barcode within MaxHammingDist away####
def findClosestBarcode(string, BarcodeList, MaxHammingDist):
    if (string not in BarcodeList) and (MaxHammingDist > 0): # if NOT in the set of barcodes, and you allow >0 Hamming distance
        k = []
        for DefinedCBC in BarcodeList:  # check for all defined BC-s
            hd = 0
            for i in range(len(DefinedCBC)): # check for all bases 1-by-1, and add 1 to hamming distance
                if DefinedCBC[i] != string[i]:
                    hd += 1
                if hd > MaxHammingDist: # Stop matching against this BC if MaxHammingDist is exceeded
                    break
            if hd <= MaxHammingDist:
                k.append(DefinedCBC)
        if len(k) == 1:
            string = k[0]  # overwrite it again if you did not find a direct match, but there is a single "hd <= MaxHammingDist" match.
    if (string not in BarcodeList): # overwrite BC it if you do not find it
        string = "NoCBCAssigned"
    return string

# Check command line inputs -------------------------------------------------------

try:
    fq1= sys.argv[1]
    protocol = sys.argv[2]
    MaximumHammingDist = int(sys.argv[3])
except:
    print "Provide 3 arguments: (1) the full path to R1 fastq file; (2) 1 or 2 according to celseq protocol; (3) Hamming distance (integer) for matching cell barcodes."
    sys.exit()

if protocol not in ['1', '2']:
    print "Either 1 or 2 is accepted for cellseq 1 or 2 as the 2nd argument."
    sys.exit()

if protocol == '1':
    bclen = 8
    umilen = 4
elif protocol == '2':
    bclen = 8
    umilen = 6

if (not MaximumHammingDist < int(bclen)):
	print "MaximumHammingDist:", MaximumHammingDist, "bclen", bclen
	print "MaximumHammingDist (in the 3rd argument) has to be an integer number smaller than the lenght of the cell barcode."
	sys.exit()

# Find fastq files. It assumes you concatenated the lanes by MapAndGo.py, thus have the name "*_cat.fastq"
fqr = fq1.split("R1_cat.fastq")[0]
fq2 = fqr + 'R2_cat.fastq'

if not os.path.isfile(fq1):
	print fq1
	print '...Fastq1 files not found'
	sys.exit()
if not os.path.isfile(fq2):
	print fq2
	print '...Fastq2 files not found'
	sys.exit()

# Run -------------------------------------------------------

fout = open(fqr + 'R2_cbc.fastq', 'w+')
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

            if protocol == "1":					# for celseq 1 (R1: celbc + umi + polyA)
                CBCseq = bcseq[:bclen]			# if bclen =4 it takes bcseq[0,1,2,3]
                UMIseq = bcseq[bclen:]
                CBCseq = findClosestBarcode(CBCseq, bccelseq1, MaximumHammingDist)        # is valid? Return closest CBC or "NoCBCAssigned"
                if CBCseq == "NoCBCAssigned":
                    badCBCs += 1
                    continue
                CBC_ID = bccelseq1[CBCseq]
            elif protocol == "2":                # for celseq 2 (R1: umi + celbc + polyA)
                CBCseq = bcseq[umilen:]
                UMIseq = bcseq[:umilen]
                CBCseq = findClosestBarcode(CBCseq, bccelseq2, MaximumHammingDist)        # is valid? Return closest CBC or "NoCBCAssigned"
                if CBCseq == "NoCBCAssigned":
                    badCBCs += 1
                    continue
                CBC_ID = bccelseq2[CBCseq]
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
            print >> fout, '\n'.join([':'.join([q, CBCseq, UMIseq, CBC_ID]), s2, p2, q2]) # write out with cell index number appended

print "		The number of input reads\t", NrReadsIn
print "		The number of output reads\t", NrReadsOut
print "		\t", round(100*NrReadsOut/NrReadsIn), '"%" of the reads have a valid CBC and UMI.'
print "		The number of incorrect cell barcodes\t", badCBCs
print "		The number of UMIs with an N\t", badUMIs


print "Concatenator.CELseq.py Finsihed"

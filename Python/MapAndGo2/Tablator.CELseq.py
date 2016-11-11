#!/usr/bin/python
# ## Usage
# - Put this script to folder where you store your scripts (e.g.: user/bin), and make it executable `chmod +X path/to/script`
# - Provide 2 arguments separated by a whitespace: InputSamFile, CelSeq protocol [1,2]
# - Test it by running in bash: `python path/to/Tablator.CELseq.py`
# - Run it with your data: python `python path/to/Tablator.CELseq.py path/to/InputSam 1` # Put 1 (or 2 if you used CELseq2 primers)
# - Original author: Anna Alemany, van Oudenaarden group, 25-10-2016
# - Modified  by Abel Vertesy, van Oudenaarden group, 25-10-2016
# ## Loaction: /home/hub_oudenaarden/avertesy/bin/Mapping/Tablator.CELseq.py

print "Tablator.CELseq.py Started"
import sys, os
import numpy as np

# def_bar = "/home/hub_oudenaarden/avertesy/var/"
def_bar = "/hpc/hub_oudenaarden/MapAndGo2/var/"
sys.path.append(def_bar)
from CELSeq1_96BC import bc2sample as bccelseq1		 # cell barcodes for CELseq 1
from CELSeq2_384BC import bc2sample as bccelseq2 	 # cell barcodes for CELseq 2

# Function and object definitions -------------------------------------------------------

class samSingleRead(object):

	def __init__(self, read):
		self.qname = read[0]
		self.flag = read[1]
		self.rname = read[2]
		self.pos = read[3]
		self.mapq = int(read[4])
		self.cigar = read[5]
		self.rnext = read[6]
		self.pnext = read[7]
		self.tlen = read[8]
		self.seq = read[9]
		self.qual = read[10]
		self.opt = {}
		if len(read) > 11:
			for i in range(11,len(read)):
				l = read[i].rsplit(':')
				if l[1] == 'i':
					self.opt[l[0]] = int(l[2])
				elif l[1] == 'f':
					self.opt[l[0]] = float(l[2])
				else:
					self.opt[l[0]] = l[2]
		return

	def isMapped(self):
		k = bin(int(self.flag))[2:].zfill(12)
		return k[-3] == '0'

	def isUniquelyMapped(self): # defined for bwa MEM but
		return self.mapq >= 20 and 'XA' not in self.opt	and 'SA' not in self.opt # If the mapping quality of a read alignment is Q, the probability that the alignment is wrong can be calculated with: probability = 10 ** -(Q/10)

	def isReverseStrand(self):
		k = bin(int(self.flag))[2:].zfill(12)
		return k[-5] == '1'


# Initialize -------------------------------------------------------

#  The readname ends in *:BC:UMI:CELL_ID , therefore the positions from the end are (0-indexing!):
bc_position = -3
umi_position = -2
cell_position = -1

try:
	InputSamFile = sys.argv[1]
	protocol = sys.argv[2]
	if not os.path.isfile(InputSamFile):
		print "samfile not found"
		sys.exit()

except:
	print "Provide 2 arguments separated by a whitespace: InputSamFile, CelSeq protocol [1,2]"
	sys.exit()

if protocol not in ['1', '2']:
	print  "Either 1 or 2 is accepted for cellseq 1 or 2 as the 2nd argument."
	sys.exit()

if protocol == '1': #### CEL Seq 1 parmeters ####
	bclen = 8
	umilen = 4
	bc2sample = bccelseq1

elif protocol == '2': #### CEL Seq 2 parmeters ####
	bclen = 8
	umilen = 6
	bc2sample = bccelseq2

K = 4**umilen

# Count -------------------------------------------------------
print "...Counting"

readscnt = {}
umicnt = {}
transcnt = {}
with open(InputSamFile) as f:
	for line in f:
		if line[0] == '@':
			continue
		r = line.rstrip()
		r = r.rsplit('\t')
		read = samSingleRead(r)
		q2 = read.qname
		q = q2.split(':')

		if read.isMapped() and read.isUniquelyMapped() and not read.isReverseStrand():
			gene = read.rname
			bc = q[len(q)-bc_position] 		# -3 for the beforelast element (-1  coz it's 0-indexed)
			if not len(bc) == bclen:
				print 'len(bc) != bclen'
				sys.exit()
			umi = q[len(q)-umi_position] 		# -2 for the beforelast element (-1  coz it's 0-indexed)
			if not len(umi) == umilen:
				print 'len(umi) != umilen'
				sys.exit()
			cell = q[len(q)-cell_position] 		# -1 for the last element (-1  coz it's 0-indexed)
			# cell = int(q[9]) # not a number any more

			try:
				umicnt[gene][cell].append(umi)
				readscnt[gene][cell] += 1
			except:
				umicnt[gene] = {c: [] for c in bc2sample.values()}
				umicnt[gene][cell].append(umi)
				readscnt[gene] = {c: 0 for c in bc2sample.values()}
				readscnt[gene][cell] = 1

# Writing out count tables -------------------------------------------------------
print "...Writing out count tables"

fc = open(InputSamFile[:-4] + '.ReadCounts.tsv', 'w+')
fb = open(InputSamFile[:-4] + '.BarcodeCounts.tsv', 'w+')
ft = open(InputSamFile[:-4] + '.TranscriptCounts.tsv', 'w+')

print >> fc, 'GENEID', '\t'.join( str(c) for c in sorted(bc2sample.values()))
print >> fb, 'GENEID', '\t'.join( str(c) for c in sorted(bc2sample.values()))
print >> ft, 'GENEID', '\t'.join( str(c) for c in sorted(bc2sample.values()))

for gene in sorted(umicnt):
	print >> fc, '\t'.join( [gene, '\t'.join(str(readscnt[gene][cell]) for cell in sorted(bc2sample.values()) ) ] )
	print >> fb, gene, '\t', '\t'.join( str(len(set(	 umicnt[gene][cell]))) for cell in sorted(bc2sample.values()) )
	t = []
	for cell in sorted(bc2sample.values()):
		x = 1.0 * len(set(umicnt[gene][cell]))
		if x > 0 and x < K:
			t.append( round(np.log(1.-x/K)/np.log(1.-1./K), 2) )
		elif x == K:
			t.append( round(np.log(1.-(K-1e-3)/K)/np.log(1.-1./K), 2) )
		elif x > K:
			print gene,
			print umicnt[gene][cell]
			print set(umicnt[gene][cell])
			print len(set(umicnt[gene][cell]))
		else:
			t.append(0)
	print >> ft, '\t', gene, '\t'.join(str(ti) for ti in t)

fc.close()
fb.close()
ft.close()

print "CountingTablator.CELseq.py Finished"

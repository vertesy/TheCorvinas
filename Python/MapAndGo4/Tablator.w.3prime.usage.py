#!/usr/bin/python
# ## Usage
# - Put this script to folder where you store your scripts (e.g.: user/bin), and make it executable `chmod +X path/to/script`
# - Provide 2 arguments separated by a whitespace: InputSamFile, CelSeq protocol [1,2]
# - Test it by running in bash: `python path/to/Tablator.CELseq.py`
# - Run it with your data: python `python path/to/Tablator.CELseq.py path/to/InputSam 1` # Put 1 (or 2 if you used CELseq2 primers)
# - Original author: Anna Alemany, van Oudenaarden group, 25-10-2016
# - Modified and extended to 3' analysis by Abel Vertesy, van Oudenaarden group, 11-06-2018
# ## Location (at HPC): /hpc/hub_oudenaarden/MapAndGo2
# ## Loaction (Local): ~/Github_repos/TheCorvinas/Python/MapAndGo4/Tablator.w.3prime.usage.py

print "Tablator.w.3prime.usage.py Started"
import sys, os
import numpy as np

def_bar = "/Users/abelvertesy/Github_repos/TheCorvinas/Python/MapAndGo2/var"
# def_bar = "/hpc/hub_oudenaarden/MapAndGo2/var/"

sys.path.append(def_bar)
from CELSeq1_96BC import bc2sample as bccelseq1		 # cell barcodes for CELseq 1
from CELSeq2_384BC import bc2sample as bccelseq2 	 # cell barcodes for CELseq 2

print "If you get an error with CELseq1: X had to be replaced in CELSeq1_96BC.py"

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
bc_position = 3
umi_position = 2
cell_position = 1

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

x=0
readscnt = {}
umicnt = {}
MappingPosList = {}
transcnt = {}
with open(InputSamFile) as f:
	for line in f:
 # These are for counting progress
		x=x+1
		i = x%100000
		if i==0:
			print "Lines processed: ",x

# Process
		if line[0] == '@':
			continue
		r = line.rstrip()
		r = r.rsplit('\t')
		read = samSingleRead(r)
		q2 = read.qname
		q = q2.split(':')

		if read.isMapped() and read.isUniquelyMapped() and not read.isReverseStrand():
			gene = read.rname
			bc = q[len(q) - bc_position] 		# -3 for the before last element (-1  coz it's 0-indexed)
			if not len(bc) == bclen:
				print 'len(bc) != bclen'
				sys.exit()
			umi = q[len(q) - umi_position] 		# -2 for the before last element (-1  coz it's 0-indexed)
			if not len(umi) == umilen:
				print 'len(umi) != umilen'
				sys.exit()
			cell = q[len(q) - cell_position] 		# -1 for the last element (-1  coz it's 0-indexed)

# Mapping Positions
			# print read.pos
			try:
				if(umi not in umicnt[gene][cell]):
					# "Only add to it if UMI is not in umicnt[gene][cell] (aka UMI collapse right when counting)"
					MappingPosList[gene][cell].append(read.pos)
			except:
				MappingPosList[gene] = {c: [] for c in bc2sample.values()}
				MappingPosList[gene][cell].append(read.pos)

# Count tables
			try:
				umicnt[gene][cell].append(umi) # "This is not really UMI count, its a list of UMI-s"
				readscnt[gene][cell] += 1
			except: 												# Exception, if you did not encounter this gene yet.
				umicnt[gene] = {c: [] for c in bc2sample.values()} # Barcode counts: Create an empty list for that gene.
				umicnt[gene][cell].append(umi)
				readscnt[gene] = {c: 0 for c in bc2sample.values()} # Read counts
				readscnt[gene][cell] = 1

cellIDs_sorted = sorted(bc2sample.values())

# Writing out count tables -------------------------------------------------------
print "...Writing out count tables"

ftrunk = InputSamFile[:-4]
fc_n = ftrunk + '.ReadCounts.tsv'
fb_n = ftrunk + '.BarcodeCounts.tsv'
ft_n = ftrunk + '.TranscriptCounts.tsv'
fmp_n = ftrunk + '.MappingPositions.tsv'


fc = open( fc_n, 'w+')
fb = open( fb_n, 'w+')
ft = open( ft_n, 'w+')
fmp = open( fmp_n, 'w+')

print >> fc, '\t'.join(['GENEID', '\t'.join( str(c) for c in cellIDs_sorted)])
print >> fb, '\t'.join(['GENEID', '\t'.join( str(c) for c in cellIDs_sorted)])
print >> ft, '\t'.join(['GENEID', '\t'.join( str(c) for c in cellIDs_sorted)])
print >> fmp, '\t'.join(['GENEID', '\t'.join( str(c) for c in cellIDs_sorted)])

for gene in sorted(umicnt):
	print gene
	print >> fc, '\t'.join([gene, '\t'.join([str(readscnt[gene][cell]) for cell in cellIDs_sorted ]) ] )
	print >> fb, '\t'.join([gene, '\t'.join( str(len(set(umicnt[gene][cell]))) for cell in cellIDs_sorted )])
	print >> fmp, '\t'.join([gene, '\t'.join( str(", ".join(MappingPosList[gene][cell])) for cell in cellIDs_sorted )])

	# Transcript count conversion ------------------------------------------------------------
	t = []
	for cell in cellIDs_sorted:
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
	# print >> ft, gene,'\t','\t'.join(str(ti) for ti in t)
	print >> ft, '\t'.join( [gene, '\t'.join(str(ti) for ti in t)] )

fc.close()
fb.close()
ft.close()
fmp.close()

print "Tablator.w.3prime.usage.py Finished"

print gene
print MappingPosList[gene]
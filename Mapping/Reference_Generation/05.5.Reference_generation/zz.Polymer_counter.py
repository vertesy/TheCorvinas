#!/usr/bin/env python
'''
### zz.Polymer_counter.py
- count poly A/C/G/T stretches of lenght N in each line of any file

'''

import os

# Setup ------------------------------------------------
N=5
Nr_of_lines=1e7
random_expectation=(0.25**N)*Nr_of_lines

# FastqDir1 = "/hpc/hub_oudenaarden/Abel/X_inact/FASTQ/";
FastqDir1 = "/hpc/hub_oudenaarden/Abel/SNP_mouse/FASTQ/";

os.chdir(FastqDir1)

# fname_Fastq = 'd5_24_af.fastq'
# fname_Fastq = 'd3_23_gm.fastq'
# fname_Fastq = 'R1/AVO55_CastB6_Lib2_R1.fastq'
fname_Fastq = 'R2/AVO55_CastB6_Lib2_R2.fastq'


bases = ['A', 'C', 'G', 'T']


# open files ------------------------------------------------
print 'Search for homo-polymers of length', N, 'in the first', Nr_of_lines,' lines of ', fname_Fastq
for base in bases:	
	print base
	polyX=base*N
	polymer_count=0

	i=0
	with open(fname_Fastq) as f:
		for line in f:
			i += 1
			# print i
			if line[0].isupper():
				if line.upper().count (polyX) > N:
					polymer_count=polymer_count+1
				if i > Nr_of_lines:
					enrichment = polymer_count / random_expectation
					print polymer_count, 'homo-',base,'-polymers are found in the first', i, 'lines, that is', enrichment, ' * enriched over random occurrence'
					# print line
					break
print "Done", i





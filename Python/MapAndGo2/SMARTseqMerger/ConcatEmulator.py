#!/usr/bin/python

import glob, os


# Setup ------------------------------------------------------
InputDir = "/Users/abelvertesy/Downloads/fasta"
inGenomeName = "hg19.fa"
OutDir = "/Users/abelvertesy/Downloads/fasta_cbc"

if not os.path.isdir(OutDir): #Check if output folders exists
		os.mkdir(OutDir)


# Init ------------------------------------------------------
os.chdir(InputDir)
filez = glob.glob('*.fastq')
print "\t",len(filez), " .fastq files are found in", InputDir
sq= range(0,len(filez))

# Read the fastq ------------------------------------------------------
ID=0
ccc = OutDir+"/"+'Concatenated.CBC.fastq'
cbc_fastq = open(ccc, "w+") 							# write to new file that will have the CBC:UMI:cellID tag


idx_file = open(OutDir+"/file_index.tsv", "w+") 		# write to new file that will have the CBC:UMI:cellID tag

for file in filez:
	print(file)
	ID+=1
	CELseqTags = "XXXXYYYY:XXXX:"+str(ID).zfill(3)
	print >> idx_file, file,"\t",CELseqTags
	print CELseqTags
	with open(file) as f:
		for line in f:
			if line.startswith('@') and '_embryo1_' in line:						# read name
				print >> cbc_fastq, line.strip()+':'+CELseqTags
			else:
				print >> cbc_fastq, line.strip()

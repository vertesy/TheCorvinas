#!/usr/bin/python
### This script counts bases in a reference file
# ~/Github_repos/TheCorvinas/Mapping/Reference_Stats/BaseFrequencies/BaseCounter.py

print "BaseCounter.py started, eyy!"
import glob
import sys, os

# Default parameters
# inputdir = "~/Github_repos/TheCorvinas/Mapping/Reference_Stats/BaseFrequencies/"
# inputfile = "head.fa"
inputdir = "~/Downloads/"
inputfile = "C.elegans.Aggregate_1003_genes_sorted_oriented_ERCC92.fa"
outputfile = "BaseFrequencies.C.elegans.tsv"

input_file_fullpath = inputdir+inputfile
print input_file_fullpath

# Function and object definitions -------------------------------------------------------
def LetterFreq(String, Letter):
	Freq = str(round(float(String.upper().count(Letter))/String.__len__(),	2))
	return Freq

# Go -------------------------------------------------------
FrequencyTable = open(inputdir + outputfile, 'w+')
print >> FrequencyTable, "Name\tA\tC\tG\tT\tN"

with open(input_file_fullpath) as f:
	for line in f:
		if line[0] == '>':
			TrName = line[1:].rstrip()
		else:
			TranscriptLength=  line.__len__()
			Az = LetterFreq(line, "A")
			Cz = LetterFreq(line, "C")
			Gz = LetterFreq(line, "G")
			Tz = LetterFreq(line, "T")
			Nz = LetterFreq(line, "N")
			AllFreq= TrName+"\t" + Az+"\t"+Cz+"\t"+Gz+"\t"+Tz+"\t"+Nz
			print >> FrequencyTable, AllFreq

# Go -------------------------------------------------------

print "Finsihed"


#!/usr/bin/python
'''
### GC_content_in_each_fasta.py
- GC content in any .reformat.fa file
	- "reformat" means, each sequence is in one line.
'''

import re

def main():
	file = open("/Users/abelvertesy/zz_local_backups/servers/UMC/bin/x_reactivation/DevCell_analysis/05.5.Reference_generation/ERCC/ERCC92.reformat.fa","r")
	# file = open("/Users/abelvertesy/zz_local_backups/servers/UMC/bin/x_reactivation/DevCell_analysis/05.5.Reference_generation/ERCC/ERCC92.reformat.fa.head","r")
	gcCount = 0
	totalBaseCount = 0
	for line in file:
		line = line.strip("\n")
		# if line.startswith(">"):
			# print line
		if not line.startswith(">"):
			gcCount = len(re.findall("[GCgc]", line))
			# print gcCount
			totalBaseCount = len(re.findall("[GCTAgcta]", line))
			# print totalBaseCount
			gcFraction = float(gcCount) / totalBaseCount
			# print(gcFraction * 100)
			print(gcFraction)
print "done"


"script is not correct"

if __name__ == '__main__':
	main()

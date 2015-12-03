#!/usr/bin/python
'''
### GC_content_in_each_fasta.py
- GC content in any .reformat.fa file
	- "reformat" means, each sequence is in one line.
source: http://saml.rilspace.org/calculating-gc-content-in-python-and-d-how-to-get-10x-speedup-in-d
'''


import re
# import string

def main():
	file = open("/Users/abelvertesy/zz_local_backups/servers/UMC/bin/x_reactivation/DevCell_analysis/05.5.Reference_generation/ERCC/ERCC92.reformat.fa","r")
	gcCount = 0
	totalBaseCount = 0
	for line in file:
		line = line.strip("\n")
		if not line.startswith(">"):
			gcCount += len(re.findall("[GCgc]", line))
			totalBaseCount += len(re.findall("[GCTAgcta]", line))
			gcFraction = float(gcCount) / totalBaseCount
			print(gcFraction * 100)
print "done"

if __name__ == '__main__':
	main()

#!/usr/bin/python
'''
### md.tsv_to_md_parser.py
- read a tab delimeted file and print the output in [github] .md format
- Usage: script tsv_filepath
'''
import csv, sys

if len(sys.argv) > 0:			# First poz File location
	tsv=sys.argv[1]
	if len(tsv) == 0:
		tsv="/Users/abelvertesy/bin/tmp.tsv"

if len(sys.argv) > 1:			# second poz False or True
	write_file=sys.argv[2]
	if len(write_file) == 0:
		write_file = False


Outfile= tsv+'.md'
print Outfile

# Execute --------------------------------------------------------
Markdown_file=open(Outfile,'w'); Markdown_file.close()				# Overwrite
Markdown_file=open(Outfile,'a')

i=0
with open(tsv, 'rU') as textfile:									# U for Universal lien ending
	for row in (list(csv.reader(textfile, delimiter='\t'))):
		i +=1
		if i ==1:
			nr_cols = len(row)

			line= '| '+" \t| ".join(row)+' |'
			if write_file: 	Markdown_file.write(line+"\n")
			print line

			line = ('|---'*nr_cols)+'|'
			if write_file: 	Markdown_file.write(line+"\n")
			print line

		else:
			line= '| '+" \t| ".join(row)+' |'
			if write_file: 	Markdown_file.write(line+"\n")
			print line


print "\n\nDONE"
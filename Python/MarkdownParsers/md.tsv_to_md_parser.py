#!/usr/bin/python
'''
### md.tsv_to_md_parser.py
- read a tab delimeted file and print the output in [github] .md format
- Usage: script tsv_filepath
'''
import csv, sys

# Defaults
tsv="/Users/abelvertesy/bin/tmp.tsv"
write_file = False

# Use Default, unless passed on
if len(sys.argv) > 1:			# First poz File location
	tsv=sys.argv[1]

if len(sys.argv) > 2:			# second poz False or True
	write_file=sys.argv[2]
	if write_file == "F":
		write_file = False

# print "write_file:", write_file
# print "tsv", tsv

Outfile= tsv+'.md'
print "\n",Outfile,"\n\n"

# Execute --------------------------------------------------------
Markdown_file=open(Outfile,'w'); Markdown_file.close()				# Overwrite
Markdown_file=open(Outfile,'a')

i=0
with open(tsv, 'rU') as textfile:									# U for Universal line ending
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
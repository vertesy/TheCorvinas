#!/usr/bin/python
'''
### md.csv_to_md_parser.py
- read a semicolon delimeted file and print the output in [github] .md format
- Usage: script tsv_filepath
'''
import csv, sys

write_file = True

tsv=sys.argv[1]
Outfile= tsv+'.md'
print Outfile

# Execute --------------------------------------------------------
Markdown_file=open(Outfile,'w'); Markdown_file.close()				# Overwrite
Markdown_file=open(Outfile,'a')

i=0
with open(tsv, 'rU') as textfile:									# U for Universal line ending
	for row in (list(csv.reader(textfile, delimiter=';'))):
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
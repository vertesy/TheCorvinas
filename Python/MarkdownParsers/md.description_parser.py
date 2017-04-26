#!/usr/bin/python
'''
### md.description_parser.py
	- Parse a github wiki page from the script headers
'''
import os, glob, re

# setup ------------------------------------------------
Outfile ="~/x_reactivation.wiki/Analysis_Methods/Pipeline-Diff-Mapping.md"
ScriptDir = "~/x_reactivation/analysis/DevCell_analysis/"
GitRepo = "https://github.com/bow/x_reactivation/"


# Prepare ------------------------------------------------
GitCode = GitRepo+"blob/master/analysis/DevCell_analysis/"
GitWikiPage = GitRepo+"wiki/"+re.split( '[/\.]', Outfile)[-2]

os.chdir(ScriptDir)
ls_Scripts = glob.glob("*.py")
ls_Scripts.extend(glob.glob("*.sh"))
ls_Scripts.extend(glob.glob("*.R"))

ls_Scripts= sorted(ls_Scripts); print ls_Scripts

Header=["# Differential Mapping Pipeline",
	"## Table of Contents"
	]

Markdown_file=open(Outfile,'w')					# Create Header
for l in Header:
	Markdown_file.writelines(l)
	Markdown_file.write("\n\n")

# Create Table of Contents ------------------------------------------------------

Markdown_file=open(Outfile,'a')
in_header = False
for s in ls_Scripts:
	with open(s) as Script:
		for line in Script:
			if line.startswith('\'\'\''):
				in_header = not in_header
			if in_header:
				if line.startswith('###'):
					script_fname = (line.split('### ')[1]).rstrip()
					entry= '- ['+script_fname+"\n"+']('+GitWikiPage+'#-'+(script_fname.replace('.','')).lower()+")\n"
					print entry
					Markdown_file.write(entry)
Markdown_file.write( "\n## Detailed Script Explanation \n")

# Creat Body ------------------------------------------------------
Markdown_file=open(Outfile,'a')
in_header = False
for s in ls_Scripts:
	with open(s) as Script:
		for line in Script:
			if line.startswith('\'\'\''):
				in_header = not in_header
			if in_header:
				if line.startswith('\'\'\''):
					Markdown_file.write("\n\n")
				else:
					print line
					if line.startswith('###'):
						entry = line.split('### ')[1].rstrip()
						line = '### [ '+entry+']('+GitCode+entry+")\n" 				# header with link
					Markdown_file.write(line)

Markdown_file.close()
print 'DONE'


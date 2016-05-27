#!/usr/bin/python -u
'''
### 00b.invite_Script_03_in_a_loop.py
- usage: replace
	- GenomeDir [modified] genomes
	- .GTF file's location
	- "RefSeq" tag in line 22
'''


import os, glob
from subprocess import call

GTF_fnp = "/Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/Data_analysis/UCSC_GTF_tables/hg19/RefSeq/03.hg19.Refseq_genes.gtf"
GenomeDir = '/Users/abelvertesy/Downloads/Genomes_Donors_Xinact'
os.chdir(GenomeDir)

Genomes = glob.glob("*.reformat.fa"); print Genomes
for Genome_name in Genomes:
	# print Genome_name
	Out_Ref_name = (Genome_name.split('reformat.fa')[0])+'RefSeq.fa'
	print Out_Ref_name


	bash_command = '/Users/abelvertesy/Mount_points/X_React/05.5.Reference_generation/03.gtf2fa.pl -in='+GTF_fnp+' -ref='+Genome_name+' > '+Out_Ref_name
	# bash_command = '/Users/abelvertesy/x_reactivation/analysis/DevCell_analysis/06.5.Reference_generation/03.gtf2fa.pl -in='+GTF_fnp+' -ref='+Genome_name+' > '+Out_Ref_name
	print bash_command
	call(bash_command, shell=True)

print "Sweet Ohio"
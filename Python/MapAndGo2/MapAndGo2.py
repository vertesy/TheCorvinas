#!/usr/bin/python
### This script automates the mapping procedure
### Put the raw sequence data in a folder, with this script
### The files should contain the minimum name of *_L00*_R*_001.fastq
### Run it with -> python MapAndGo.py
### For usage run -> python MapAndGo.py -help
# - Authors: Abel Vertesy and Thom de Hoog van Oudenaarden group, 02-03-2015
# - Maintained by Abel Vertesy, van Oudenaarden group, 25-04-2017

import glob
import sys
import os

# Default parameters
def_email = "<none>"
# def_email = "a.vertesy@hubrecht.eu"
def_CelSeqPrimerVersion = "2"
def_MaxHammingDist = "1"
def_ScriptsFolder = "/hpc/hub_oudenaarden/MapAndGo2/"
def_BWA_Folder = "/hpc/hub_oudenaarden/bin/software/bwa-0.7.10/"
def_bar = "/home/hub_oudenaarden/avertesy/var/"
def_bash_out = 'bash_files'
def_counts_out = 'count_files'
def_map_out = 'map_files'
def_logfiles_out = 'logs_and_errors'

def_refseq_human = "/hpc/hub_oudenaarden/gene_models/human_gene_models/hg19_mito/hg19_RefSeq_genes_clean_ERCC92_polyA_10_masked_Mito.fa"
def_refseq_mouse = "/hpc/hub_oudenaarden/gene_models/mouse_gene_models/mm10_eGFP_mito/mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked_eGFP_Mito.fa"
def_refseq_zebrafish = "/hpc/hub_oudenaarden/gene_models/zebrafish_gene_models/Danio_rerio_Zv9_ens74_extended3_genes_ERCC92.fa"
def_refseq_elegans = "/hpc/hub_oudenaarden/gene_models/cel_gene_models/Aggregate_1003_genes_sorted_oriented_ERCC92.fa"
def_refseq_briggsae = "/hpc/hub_oudenaarden/gene_models/cbr_gene_models/cb3_transcriptome_ERCC92.fa"

def_qsub = "no"
def_mem = "8"
def_time = "06:00:00"
def_core = "4"
def_unzip = "yes"
def_zip = "yes"

# Define necessary variables
argv_imput = {}
argv_valid = ["-help", "-email", "-qsub", "-bar", "-ref", "-mem", "-cores", "-time", "-CelSeqPrimerVersion", "-MaxHammingDist", "-ScriptsFolder", "-BWA_Folder", "-unzip", "-zip", "-bash_out", "-counts_out", "-map_out", "-logs_out"] # , "-cat_out"
bash_file_name = []
do_break = False
files = []
options = ["yes","no"]
params = {}

# Check if required default input is valid
if not def_qsub.lower() in options:
		sys.exit("\nThe default for the argument -qsub= is not supported!\n")
if not def_unzip.lower() in options:
		sys.exit("\nThe default for the argument -unzip= is not supported!\n")
if not def_zip.lower() in options:
		sys.exit("\nThe default for the argument -zip= is not supported!\n")

# Make help file
help_text = ( "\nThe script automates the mapping procedure by making bash files calling the cat, Concatenator.CELseq.py, bwa mem and the Tablator.CELseq.py scrips. \nIt works with .fastq files in the current working directory that have the minimal name of *_L00*_R*_*.fastq, \nand can submit them to the desired queue on the HPC. The following options are available for use:" +
"\n\n-help \n\n\tShows help!" +
"\n\n-qsub= \n\n\tThis determines if the generated bash file will be submitted on the queue. 	\n\tOptions: \"yes,\"no\". Default: \"no\"." +
"\n\n-email= \n\n\tThis specifies the email address to be used for the do_mapping.pl command. \n\tOption: [email address]. Default: no email." +
"\n\n-ref= \n\n\tThis specifies the refseq file to be used for the do_mappings_strand.pl command. \n\tOptions: [refseq path] or, \n\t\t \"human\" which uses: " + def_refseq_human + ", \n\t\t \"mouse\" which uses: " + def_refseq_mouse + ", \n\t\t \"zebrafish\" which uses: " + def_refseq_zebrafish + " \n\tAlso C elegans and briggsae default gene models are present under these names. Default: \"NONE\"." +
# "\n\n-bar= \n\n\tThis specifies the barcode file to be used for the extract_counts_rb.pl command. \n\tOptions: [barcode path]. Default: " + def_bar + "." +
"\n\n-CelSeqPrimerVersion= \n\n\t Using CELseq1 means read 1 starts with 8n cell barcode (cbc) + 4n UMI. CELseq2 means 6n UMI + 8 cbc \n\tOptions: [1 or 2]. Default: 1" +
"\n\n-MaxHammingDist= \n\n\t How many mismatches are allowed in a cell barcode. If there is a single barcode maximum MaxHammingDist away, it counts it. \n\tOptions: [integer<4]. Default: 1" +
"\n\n-ScriptsFolder= \n\n\t  Where are the MapAndGo2 sripts located.\n\tOptions: [folder name]. Default: /hpc/hub_oudenaarden/MapAndGo2/" +
"\n\n-BWA_Folder= \n\n\t Where is your version of BWA located\n\tOptions: [folder name]. Default: /hpc/hub_oudenaarden/bin/software/bwa-0.7.10/" +
"\n\n-mem= \n\n\t Memory to reserve \n\tOptions: [GB, integer]. Default: 5" +
"\n\n-time= \n\n\t Max estimated running time. If it runs longer, HPC kills the job \n\tOptions: [hours, HH:MM:SS]. Default: 12:00:00" +
"\n\n-cores= \n\n\t Number of cores reserved. The more cores you reserve the longer your job will wait to start, but then the faster it runs. \n\tOptions: [integer]. Default: 4" +
# "\n\n-cat_out= \n\n\tThis specifies the output folder for the cat command and makes the directory if it did not exist yet.\n\tOptions: [folder name]. Default: current working directory." +
"\n\n-bash_out= \n\n\tThis specifies the output folder for the bash files and makes the directory if it did not exist yet.\n\tOptions: [folder name]. Default: current working directory." +
"\n\n-map_out= \n\n\tThis specifies the output folder for the do_mapping.pl command and makes the directory if it did not exist yet.\n\tOptions: [folder name]. Default: current working directory." +
"\n\n-counts_out= \n\n\tThis specifies the output folder for the extract_counts.pl command and makes the directory if it did not exist yet.\n\tOptions: [folder name]. Default: current working directory." +
"\n\n-def_logfiles_out= \n\n\tThis specifies the output folder for Logfiles.\n\tOptions: [folder name]. Default: current working directory." +
"\n\n-unzip= \n\n\tThis determines if .fastq files to use are zipped and need to be unzipped. Options: \"yes\",\"no\". Default: \"no\"." +
"\n\n-zip= \n\n\tThis determines if the used .fastq files should be zipped after mapping. Options: \"yes\",\"no\". Default: \"no\"." +
"\n\nThe Defaults: \n\n\t" +
"\n\nExample: MapAndGo2.py -ref=mouse -qsub=no -unzip=yes -zip=no -mem=5 -time=12:00:00 -email=a.vertesy@hubrech.eu -CelSeqPrimerVersion=2 -cores=4 -MaxHammingDist=1 -bash_out='bash_files' -counts_out='count_files' -map_out='map_files'" +
"\n\nImportant: Be careful to check the bash files!\n" )

# Readout the command line arguments ----
for i in range(1, len(sys.argv)):
	if sys.argv[i].lower() == '-help':
		sys.exit(help_text)

	split = sys.argv[i].split("=")
	if (len(split) == 2 and split[1] != '' and split[0] in argv_valid):
			argv_imput[split[0]] = split[1]
	else:
		if do_break:
			break_message = break_message + "\nThe argument \"" + sys.argv[i] + "\" is not supported!"
		else:
			break_message = "\nThe argument \""+ sys.argv[i] + "\" is not supported!"
		do_break = True

if do_break:
	sys.exit(break_message + "\nFor help type: MapAndGo.py -help\n")

# Currently fixed (ignores what you put in)
params["-ScriptsFolder"] = def_ScriptsFolder
params["-BWA_Folder"] = def_BWA_Folder

# Process parameter for CelSeqPrimerVersion
if "-CelSeqPrimerVersion" in argv_imput:
	if argv_imput["-CelSeqPrimerVersion"] in ["1","2"]:
		params["-CelSeqPrimerVersion"] = argv_imput["-CelSeqPrimerVersion"]
	else:
		print "-CelSeqPrimerVersion should be either 1 or 2."
else:
	params["-CelSeqPrimerVersion"] =  def_CelSeqPrimerVersion

# Process parameter for MaxHammingDist
if "-MaxHammingDist" in argv_imput:
	if (int(argv_imput["-MaxHammingDist"]) <= 4):
		params["-MaxHammingDist"] = argv_imput["-MaxHammingDist"]
	else:
		print "-MaxHammingDist should be either a a number 1-4."
else:
	params["-MaxHammingDist"] =  def_MaxHammingDist

# Process parameter for barcode file
# if "-bar" in argv_imput:
# 	if  os.path.isfile(argv_imput["-bar"]):
# 		params["-bar"] = argv_imput["-bar"]
# 	else:
# 		sys.exit("\nThe barcode file does not exist!\n")
# else:
# 	if not os.path.isdir(def_bar):
# 		sys.exit("\nThe default barcodes file does not exist!\n")
# 	else:
# 		params["-bar"] =  def_bar

#Process parameter for bash_out
if "-bash_out" in argv_imput:
	bash_out = argv_imput["-bash_out"]
	if bash_out[(len(bash_out)-1)] == "/":
		bash_out = bash_out[0:(len(bash_out)-1)]
	params["-bash_out"] = os.getcwd() + "/" + bash_out
else:
	if def_bash_out[(len(def_bash_out)-1)] == "/":
		def_bash_out = def_bash_out[0:(len(def_bash_out)-1)]
	params["-bash_out"] = os.getcwd() + "/" + def_bash_out

#Process parameter for counts_out
if "-counts_out" in argv_imput:
	counts_out = argv_imput["-counts_out"]
	if counts_out[(len(counts_out)-1)] == "/":
		counts_out = counts_out[0:(len(counts_out)-1)]
	params["-counts_out"] = os.getcwd() + "/" + counts_out
else:
	if def_counts_out[(len(def_counts_out)-1)] == "/":
		def_counts_out = def_counts_out[0:(len(def_counts_out)-1)]
	params["-counts_out"] = os.getcwd() + "/" + def_counts_out

#Process parameter for map_out
if "-map_out" in argv_imput:
	map_out = argv_imput["-map_out"]
	if map_out[(len(map_out)-1)] == "/":
		map_out = map_out[0:(len(map_out)-1)]
	params["-map_out"] = os.getcwd() + "/" + map_out
else:
	if def_map_out[(len(def_map_out)-1)] == "/":
		def_map_out = def_map_out[0:(len(def_map_out)-1)]
	params["-map_out"] = os.getcwd() + "/" + def_map_out

#Process parameter for logs_out
if "-logs_out" in argv_imput:
	logs_out = argv_imput["-logs_out"]
	if logs_out[(len(logs_out)-1)] == "/":
		logs_out = logs_out[0:(len(logs_out)-1)]
	params["-logs_out"] = os.getcwd() + "/" + logs_out
else:
	if def_logfiles_out[(len(def_logfiles_out)-1)] == "/":
		def_logfiles_out = def_logfiles_out[0:(len(def_logfiles_out)-1)]
	params["-logs_out"] = os.getcwd() + "/" + def_logfiles_out

#Check if defaults folders exists
if not os.path.isdir(params["-bash_out"]):
		os.mkdir(params["-bash_out"])
if not os.path.isdir(params["-counts_out"]):
		os.mkdir(params["-counts_out"])
if not os.path.isdir(params["-map_out"]):
		os.mkdir(params["-map_out"])
# if not os.path.isdir(params["-logs_out"]):
# 		os.mkdir(params["-logs_out"])

#Set parameter for refseq file ------------------------------------------------------------------------------------------------------------
if "-ref" in argv_imput:
	if argv_imput["-ref"] == "human":
		if os.path.isfile(def_refseq_human):
			params["-ref"] = def_refseq_human
		else:
			sys.exit("\nThe default human refseq file is not found! Check if the file exist and if the right path is given!\n")

	elif argv_imput["-ref"] == "mouse":
		if os.path.isfile(def_refseq_mouse):
			params["-ref"] = def_refseq_mouse
		else:
			sys.exit("\nThe default mouse refseq file is not found! Check if the file exist and if the right path is given!\n")

	elif argv_imput["-ref"] == "zebrafish":
		if os.path.isfile(def_refseq_zebrafish):
			params["-ref"] = def_refseq_zebrafish
		else:
			sys.exit("\nThe default zebrafish refseq file is not found! Check if the file exist and if the right path is given!\n")

	elif argv_imput["-ref"] == "elegans":
		if os.path.isfile(def_refseq_elegans):
			params["-ref"] = def_refseq_elegans
		else:
			sys.exit("\nThe default C. elegans refseq file is not found! Check if the file exist and if the right path is given!\n")

	elif argv_imput["-ref"] == "briggsae":
		if os.path.isfile(def_refseq_briggsae):
			params["-ref"] = def_refseq_briggsae
		else:
			sys.exit("\nThe default C. briggsae refseq file is not found! Check if the file exist and if the right path is given!\n")

	else:
		if os.path.isfile(argv_imput["-ref"]):
					params["-ref"] = argv_imput["-ref"]
		else:
			sys.exit("\nThe refseq file is not found! Check if the file exist and if the right path is given! \n")
else:
	print "Specify mapping reference in '-ref' ! Either write human, mouse, zebrafish, elegans, briggsae, or provide the path"

#Set parameter for mem
if "-mem" in argv_imput:
	RequestedRAM = argv_imput["-mem"].lower()
	# sys.exit("\nMem has to be a whole number denoting the gigabytes requested: 10 \n")
else:
	RequestedRAM = def_mem
params["-mem"] = "#$ -l h_vmem=" + RequestedRAM + "G"

#Set parameter for time
if "-time" in argv_imput:
	runtime = argv_imput["-time"].lower()
	# sys.exit("\nTime has to be in this format: 12:00:00\n")
else:
	runtime = def_time
params["-time"] = "#$ -l h_rt=" + runtime

# Set parameter for threads
if "-cores" in argv_imput:
	params["-cores"] = argv_imput["-cores"]
else:
	params["-cores"] = def_core

#Set parameter for unzip
if "-unzip" in argv_imput:
		if argv_imput["-unzip"].lower() in options:
				params["-unzip"] = argv_imput["-unzip"].lower()
		else:
				sys.exit("\nOnly \"yes\" and \"no\" are supported for -unzip\n")
else:
		params["-unzip"] = def_unzip

#Set parameter for zip
if "-zip" in argv_imput:
	if argv_imput["-zip"].lower() in options:
		params["-zip"] = argv_imput["-zip"].lower()
	else:
		sys.exit("\nOnly \"yes\" and \"no\" are supported for -zip\n")
else:
	params["-zip"] = def_zip

#Process parameter for email
if "-email" in argv_imput:
	emailAddress = argv_imput["-email"]
else:
	emailAddress = def_email
params["-email"] = "\n#$ -M " + emailAddress

#Set parameter for qsub
if "-qsub" in argv_imput:
	if argv_imput["-qsub"].lower() in options:
		params["-qsub"] = argv_imput["-qsub"].lower()
	else:
		sys.exit("\nOnly \"yes\" and \"no\" are supported for -qsub\n")
else:
	params["-qsub"] = def_qsub

# Head of the bash files ------------------------------------------------------------------------------------------------------------------------------------------------
head_bash = "#! /bin/bash" + "\n#$ -cwd \n#$ -V \n" + params["-time"] +"\n" + params["-mem"] + params["-email"] +" \n#$ -m beas \n#$ -pe threaded "  + params["-cores"] + "\n\n"

# Get unique files in dir
if params["-unzip"] == "yes":
	gz=".gz"
	Files = glob.glob("*_L00*_R*_*.fastq.gz")
	if len(Files) == 0:
		print "\nThere are no *_L00*_R*_*.fastq.gz files in the current working directory with the minimal name of *_L00*_R*_*.fastq\nTrying to find unzipped files.\nparams['-unzip'] is set to NO."
		params["-unzip"] = "no"
		Files = glob.glob("*_R*_*.fastq")
		if len(Files) == 0:
			sys.exit("\nThere are no *_R*_*.fastq files in the current working directory with the minimal name of *_L00*_R*_*.fastq\n")
# Try unzipped files
else:
	gz=""
	Files = glob.glob("*_R*_*.fastq")
	if len(Files) == 0:
		sys.exit("\nThere are no *_R*_*.fastq files in the current working directory with the minimal name of *_L00*_R*_*.fastq\n")

for i in range(0, len(Files)):
	name = Files[i].split("_L00")

	if(len(name) == 2):
		files.append(name[0])

files = list(set(files))

#Check if all files exist
rs = ["1","2"] 				# We do paired end sequencing
ls = ["1","2","3","4"] 		# There are four lanes in a Illumina NextSeq run
missing_files = ""

for i in range(0, len(files)):
	for r in rs:
		for l in ls:
			path = files[i] + "_L00" + l  +"_R" + r + "_001.fastq"+gz
			if not os.path.exists(os.getcwd()+"/"+path):
				if missing_files == "":
					missing_files = path + " does not exist!"
				else:
					missing_files = missing_files + "\n" + path + " does not exist!"

if not missing_files == "":
		sys.exit("\n"+missing_files+"\n")

used_files = ""
for i in range(0, len(files)):
		if used_files == "":
				used_files = files[i]
		else:
				used_files = used_files + "\n\t\t\t\t\t" + files[i]

#Output variables
print "\n.fastq files recognised: 			" + str(len(Files))
print ".fastq files with unique index: 		" + str(len(files))
print "Unzip and use *.fastq.gz files:			" + params['-unzip']
print "Names of files of used:				" + used_files
print "Used RefSeq file: 				" + params['-ref']
# print "Used barcode file: 				" + params['-bar']
# print "Output dir for the concatenated .fastq files:	" + params['-cat_out'] + "/"
print "Output dir for bash files to be submitted: 	" + params['-bash_out'] + "/"
print "Output dir for sam files: 			" + params['-map_out'] + "/"
print "Output dir for count tables: 			" + params['-counts_out'] + "/"
print "Output dir for log & error files: 		" + params['-logs_out'] + "/"
print "The number of used cores: 			" + params['-cores']
print "Submission to queue: 				" + params['-qsub']
print "Memory requested: 				" + RequestedRAM
print "Time requested: 				" + runtime
print "Email: 						" + emailAddress
print "zip .fastq files afterwards:			" + params['-zip']

#Make unique filenames
print "\n\n\t\tIMPORTANT: Make sure to check the bash files:\n"

for i in range(0,len(files)):
	nr = 1
	file_name = params['-bash_out'] + "/map_"+files[i]+"_"+str(nr)+".sh"

	while os.path.exists(file_name):
		split = file_name.split("_")[len(file_name.split("_"))-1]
		nr = int(split.split(".sh")[0]) + 1
		file_name = params['-bash_out'] + "/map_"+files[i]+"_"+str(nr)+".sh"
	print "\t\t\t\t"+file_name

	bash_file_name.append(file_name)

print "\n\n"

# Make commands ------------------------------------------------------------------------------------------------------------------------
for i in range(0,len(files)):
	NameRoot = os.getcwd() +  "/" + files[i]
	ConcatRoot = params["-map_out"] + "/" + files[i]
	FastQ1_cat = ConcatRoot + "_R1_cat.fastq"
	FastQ2_cat = ConcatRoot + "_R2_cat.fastq"
	FastQ1_cbc = ConcatRoot + "_R2_cbc.fastq"
	FastQ1_sam = ConcatRoot + "_R2.sam"

	bash_file = open(bash_file_name[i], "w")
	if params["-unzip"] == "yes":
		unzip = "gunzip " + NameRoot + "_L00*_R*_*.fastq.gz" + "\n\n"
	else:
		print "... not unzipping, see first line!\n"
		unzip = ""

	cat_r1 = "cat " + NameRoot + "_L00*_R1* > " + FastQ1_cat + "\n\n"
	cat_r2 = "cat " + NameRoot + "_L00*_R2* > " + FastQ2_cat + "\n\n"

	MakeSingleEnd = params["-ScriptsFolder"] + "Concatenator.CELseq.py " + FastQ1_cat + " " + params["-CelSeqPrimerVersion"] + " " +  params["-MaxHammingDist"] + "\n\n"
	mapping = params["-BWA_Folder"] + "bwa mem -t " + params['-cores'] + " " + params['-ref'] + " " + FastQ1_cbc + " > " + FastQ1_sam + "\n\n"
	extract =  params["-ScriptsFolder"] + "Tablator.CELseq.py " + FastQ1_sam + " " + params["-CelSeqPrimerVersion"] + "\n\n"
	Copy = "mv " + params["-map_out"] + "/*.tsv " + params["-counts_out"] + "/\n\n\n"
	Mappability = "samtools flagstat " + FastQ1_sam +" >" +  params["-map_out"] + " flagstat.mappability.txt" + "\n\n"
	DownloadResults = "# rsync -avzP UMC: " + params["-counts_out"] + " " + "~/Downloads\n"

	if params['-zip'] == "yes":
		ZipUp = "\n \ngzip " + os.getcwd() + "/" + files[i] + "*.fastq"
		if params['-map_out'] != os.getcwd():
			ZipUp = ZipUp + "\n \ngzip " + params['-map_out'] + "/" + files[i] + "*.fastq"
	else:
		ZipUp = ""

	bash_text = head_bash + unzip + cat_r1 + cat_r2 + MakeSingleEnd + mapping +extract + ZipUp + Copy + Mappability + DownloadResults
	# bash_text = head_bash + MakeSingleEnd + mapping +extract + ZipUp + DownloadResults
	bash_file.write(bash_text)
	bash_file.close()

# Submit to queue ------------------------------------------------------------------------------------------------------------------------
if params['-qsub'] == 'yes':
	for i in range(0, len(files)):
		command = "qsub " + bash_file_name[i]
		os.system(command)
		os.system("qstat")

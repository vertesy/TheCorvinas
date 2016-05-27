#!/usr/bin/python
# ## Usage
# - Put the raw sequence data in a folder, put this script to folder where you store your scripts (e.g.: user/bin), and make it executable `chmod +X path/to/MapAndGo.py`
# - The files should contain the minimum name of *_L00*_R*_001.fastq (that is Illumina NextSeq's default output: R for read 1 or 2, L for one of the 4 lanes)
# - Test it by running in bash: `python path/to/MapAndGo.py -help`
# - Run it with your data: python `path/to/MapAndGo.py MapAndGo.py -ref=human -bar=cel-seq_barcodes.csv -bash_out=bash_files -cat_out=cat_files  -map_out=map_files -counts_out=count_files -email=x.y@hubrecht.eu`
# - Original author: Thom de Hoog, van Oudenaarden group, 02-03-2015
# - Modified and maintained by Abel Vertesy, van Oudenaarden group, 28-04-2016

import glob
import sys
import os

#Default parameters
default_bar = "/hpc/hub_oudenaarden/data/cel-seq_barcodes.csv"
default_bash_out = os.getcwd()
default_cat_out = os.getcwd()
default_counts_out = os.getcwd()
default_email = "<none>"
default_map_out = os.getcwd()
default_refseq_human = "/hpc/hub_oudenaarden/gene_models/human_gene_models/hg19_RefSeq_genes_clean_ERCC92_polyA_10_masked.fa"
default_refseq_mouse = "/hpc/hub_oudenaarden/gene_models/mouse_gene_models/mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked.fa"
default_refseq_zebrafish = "/hpc/hub_oudenaarden/gene_models/zebrafish_gene_models/Danio_rerio_Zv9_ens74_extended3_genes_ERCC92.fa"
default_qsub = "no"
default_mem = "5"
default_time = "24:00:00"
default_core = "8"
# default_qtype = "long"
default_unzip = "no"
default_zip = "no"


#Define necessary variables
argv_imput = {}
# argv_valid = ["-bar","-bash_out","-counts_out","-cat_out","-email","-help","-map_out","-ref","-qsub","-qtype","-unzip","-zip"]
argv_valid = ["-bar","-bash_out","-counts_out","-cat_out","-email","-help","-map_out","-ref","-qsub","-unzip","-zip", "-mem", "-time"]
bash_file_name = []
# cores = {"veryshort":"12","short":"9","medium":"7","long":"4","verylong":"1"}
do_break = False
files = []
options = ["yes","no"]
params = {}


#Check if defaults folders exists
if not os.path.isdir(default_bash_out):
		os.mkdir(default_bash_out)
if not os.path.isdir(default_cat_out):
		os.mkdir(default_cat_out)
if not os.path.isdir(default_counts_out):
		os.mkdir(default_counts_out)
if not os.path.isdir(default_map_out):
		os.mkdir(default_map_out)

#Check if required default input is valid
if not default_qsub.lower() in options:
		sys.exit("\nThe default for the argument -qsub= is not supported!\n")
# if not default_qtype.lower() in cores:
	# sys.exit("\nThe default for the argument -qtype= is not supported!\n")
if not default_unzip.lower() in options:
		sys.exit("\nThe default for the argument -unzip= is not supported!\n")
if not default_zip.lower() in options:
		sys.exit("\nThe default for the argument -zip= is not supported!\n")

#Make help file
line1 ="\nThe script automates the mapping procedure by making bash files containing the cat, do_mapping_strand.pl and extract_counts_rb.pl commands. \nIt works with .fastq files in the current working directory that have the minimal name of *_L00*_R*_*.fastq, \nand can submit them to the desired queue on the HPC. The following options are available for use:"
line2 ="\n\n-qsub= \n\n\tThis determines if the generated bash file will be submitted on the queue. 	\n\tOptions: \"yes,\"no\". Default: \"no\"."
# line3 ="\n\n-qtype= \n\n\tThis determines the type of queue type of queue to be used. \n\tOptions: \"veryshort\" (12 cores), \"short\" (9 cores), \"medium\" (7 cores), \n\t\t \"long\" (4 cores) and \"verylong\" (1 core). \n\tDefault: \"long\" (4 cores). \n"
line4 ="\n\n-ref= \n\n\tThis specifies the refseq file to be used for the do_mappings_strand.pl command. \n\tOptions: [refseq path] or, \n\t\t \"human\" which uses: " + default_refseq_human + ", \n\t\t \"mouse\" which uses: " + default_refseq_mouse + ", \n\t\t \"zebrafish\" which uses: " + default_refseq_zebrafish + " \n\tDefault: \"mouse\"."
line5 ="\n\n-bar= \n\n\tThis specifies the barcode file to be used for the extract_counts_rb.pl command. \n\tOptions: [barcode path]. Default: " + default_bar + "."
line6 ="\n\n-cat_out= \n\n\tThis specifies the output folder for the cat command and makes the directory if it did not exist yet.\n\tOptions: [folder name]. Default: current working directory."
line7 ="\n\n-bash_out= \n\n\tThis specifies the output folder for the bash files and makes the directory if it did not exist yet.\n\tOptions: [folder name]. Default: current working directory."
line8 ="\n\n-map_out= \n\n\tThis specifies the output folder for the do_mapping.pl command and makes the directory if it did not exist yet.\n\tOptions: [folder name]. Default: current working directory."
line9 ="\n\n-counts_out= \n\n\tThis specifies the output folder for the extract_counts.pl command and makes the directory if it did not exist yet.\n\tOptions: [folder name]. Default: current working directory."
line10 ="\n\n-email= \n\n\tThis specifies the email address to be used for the do_mapping.pl command. \n\tOption: [email address]. Default: no email."
line11 =" \n\n-unzip= \n\n\tThis determines if .fastq files to use are zipped and need to be unzipped. Options: \"yes\",\"no\". Default: \"no\"."
line12 ="\n\n-zip= \n\n\tThis determines if the used .fastq files should be zipped after mapping. Options: \"yes\",\"no\". Default: \"no\"."
line13 ="\n\n-help \n\n\tShows help!"
line14 ="\n\nSummary of the defaults: \n\n\t"
# line15 ="\n\nExample: make_bash.py -qsub=no -qtype=medium -ref=human -bar=cel-seq_barcodes.csv -bash_out=bash_files -cat_out=cat_files  -map_out=map_files -counts_out=count_files -email=t.hoog@hubrecht.eu"
line15 ="\n\nExample: MapAndGo.py -qsub=no -ref=human -bar=cel-seq_barcodes.csv -bash_out=bash_files -cat_out=cat_files  -map_out=map_files -counts_out=count_files -email=x.y@hubrecht.eu"
line16 ="\n\nImportant: Be careful to check the bash files!\n"
# help_text = line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8 + line9 + line10 + line11 + line12 + line13 + line14 + line15 + line16
help_text = line1 + line2 + line4 + line5 + line6 + line7 + line8 + line9 + line10 + line11 + line12 + line13 + line14 + line15 + line16


#Readout the command line arguments
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


#Process parameter for barcode file
if "-bar" in argv_imput:
	if  os.path.isfile(argv_imput["-bar"]):
		params["-bar"] = argv_imput["-bar"]
	else:
		sys.exit("\nThe barcode file does not exist!\n")
else:
	if not os.path.isfile(default_bar):
		sys.exit("\nThe default barcodes file does not exist!\n")
	else:
		params["-bar"] =  default_bar


#Process parameter for bash_out
if "-bash_out" in argv_imput:
	bash_out = argv_imput["-bash_out"]
	if bash_out[(len(bash_out)-1)] == "/":
		argv_imput["-bash_out"] = bash_out[0:(len(bash_out)-1)]

	params["-bash_out"] = os.getcwd() + "/" + argv_imput["-bash_out"]
	if not os.path.isdir(argv_imput["-counts_out"]):
		os.mkdir(params["-bash_out"])
else:
	if default_bash_out[(len(default_bash_out)-1)] == "/":
		default_bash_out = default_bash_out[0:(len(default_bash_out)-1)]
	params["-bash_out"] = default_bash_out

#Process parameter for cat_out
if "-cat_out" in argv_imput:
		cat_out = argv_imput["-cat_out"]
		if cat_out[(len(cat_out)-1)] == "/":
				argv_imput["-cat_out"] = cat_out[0:(len(cat_out)-1)]

		params["-cat_out"] = os.getcwd() + "/" + argv_imput["-cat_out"]
		if not os.path.isdir(argv_imput["-cat_out"]):
				os.mkdir(params["-cat_out"])
else:
		if default_cat_out[(len(default_cat_out)-1)] == "/":
				default_cat_out = default_cat_out[0:(len(default_cat_out)-1)]
		params["-cat_out"] = default_cat_out

#Process parameter for counts_out
if "-counts_out" in argv_imput:
	counts_out = argv_imput["-counts_out"]
	if counts_out[(len(counts_out)-1)] == "/":
		argv_imput["-counts_out"] = counts_out[0:(len(counts_out)-1)]

	params["-counts_out"] = os.getcwd() + "/" + argv_imput["-counts_out"]
	if not os.path.isdir(argv_imput["-counts_out"]):
		os.mkdir(params["-counts_out"])
else:
	if default_counts_out[(len(default_counts_out)-1)] == "/":
		default_counts_out = default_counts_out[0:(len(default_counts_out)-1)]
	params["-counts_out"] = default_counts_out


#Process parameter for email
if "-email" in argv_imput:
	params["-email"] = "\n#$ -M " + argv_imput["-email"]
	params["-email_fb"] = argv_imput["-email"]
else:
	params["-email"] = "\n#$ -M " + default_email
	params["-email_fb"] = default_email


#Process parameter for map_out
if "-map_out" in argv_imput:

	map_out = argv_imput["-map_out"]

	if map_out[(len(map_out)-1)] == "/":
		argv_imput["-map_out"] = map_out[0:(len(map_out)-1)]

	params["-map_out"] = os.getcwd() + "/" +  argv_imput["-map_out"]
	if not os.path.isdir(argv_imput["-map_out"]):
		os.mkdir(params["-map_out"])
else:
	if default_map_out[(len(default_map_out)-1)] == "/":
		default_map_out = default_map_out[0:(len(default_map_out)-1)]
	params["-map_out"] = default_map_out


#Set parameter for refseq file
if "-ref" in argv_imput:
	if argv_imput["-ref"] == "human":
		if os.path.isfile(default_refseq_human):
			params["-ref"] = default_refseq_human

		else:
			sys.exit("\nThe default human refseq file is not found! Check if the file exist and if the right path is given!\n")

	elif argv_imput["-ref"] == "mouse":
		if os.path.isfile(default_refseq_mouse):
			params["-ref"] = default_refseq_mouse

		else:
			sys.exit("\nThe default mouse refseq file is not found! Check if the file exist and if the right path is given!\n")

	elif argv_imput["-ref"] == "zebrafish":
		if os.path.isfile(default_refseq_zebrafish):
			params["-ref"] = default_refseq_zebrafish

		else:
			sys.exit("\nThe default zebrafish refseq file is not found! Check if the file exist and if the right path is given!\n")

	else:
		if os.path.isfile(argv_imput["-ref"]):
					params["-ref"] = argv_imput["-ref"]
		else:
			sys.exit("\nThe refseq file is not found! Check if the file exist and if the right path is given! \n")
else:
	if not os.path.isfile(default_refseq_mouse):
		sys.exit("\nThe default refseq file is not found! Check if the file exist and if the right path is given!\n")
	else:
		params["-ref"] =  default_refseq_mouse


#Set parameter for qsub
if "-qsub" in argv_imput:
	if argv_imput["-qsub"].lower() in options:
		params["-qsub"] = argv_imput["-qsub"].lower()
	else:
		sys.exit("\nOnly \"yes\" and \"no\" are supported for -qsub\n")
else:
	params["-qsub"] = default_qsub


# #Set parameter for qtype
# if "-qtype" in argv_imput:
# 	if argv_imput["-qtype"].lower() in cores:
# 		params["-qtype"] = argv_imput["-qtype"]
# 		params["-cores"] = cores[argv_imput["-qtype"]]
# 	else:
# 		sys.exit("\nOnly \"veryshort\", \"short\", \"medium\", \"long\" and \"verylong\" are supported for -qtype\n")
# else:
# 	params["-qtype"] = default_qtype
# 	params["-cores"] = cores[default_qtype]


#Set parameter for mem
if "-mem" in argv_imput:
	params["-mem"] = argv_imput["-mem"].lower()
	# sys.exit("\nMem has to be a whole number denoting the gigabytes requested: 10 \n")
else:
	params["-mem"] = "#$ -l h_vmem=" + default_mem + "G"

#Set parameter for time
if "-time" in argv_imput:
	params["-time"] = argv_imput["-time"].lower()
	# sys.exit("\nTime has to be in this format: 12:00:00\n")
else:
	params["-time"] = "#$ -l h_rt=" + default_time

# Set parameter for threads
if "-cores" in argv_imput:
	params["-cores"] = argv_imput["-cores"]
else:
	params["-cores"] = default_core

#Set parameter for unzip
if "-unzip" in argv_imput:
		if argv_imput["-unzip"].lower() in options:
				params["-unzip"] = argv_imput["-unzip"].lower()
		else:
				sys.exit("\nOnly \"yes\" and \"no\" are supported for -unzip\n")
else:
		params["-unzip"] = default_unzip

#Set parameter for zip
if "-zip" in argv_imput:
	if argv_imput["-zip"].lower() in options:
		params["-zip"] = argv_imput["-zip"].lower()
	else:
		sys.exit("\nOnly \"yes\" and \"no\" are supported for -zip\n")
else:
	params["-zip"] = default_zip


#Head of the bash files
# head_bash = "#! /bin/bash \n#$ -q " + params["-qtype"] + "\n#$ -cwd \n#$ -V " + params["-email"] +" \n#$ -m beas \n#$ -pe threaded " + params["-cores"]
head_bash = "#! /bin/bash" + "\n#$ -cwd \n#$ -V \n" + params["-time"] +"\n" + params["-mem"] + params["-email"] +" \n#$ -m beas \n#$ -pe threaded "  + params["-cores"]

#Get unique files in dir
if params["-unzip"] == "yes":
		gz = ".gz"
else:
		gz = ""

#Get unique files in dir
dir = glob.glob("*_L00*_R*_*.fastq"+gz)

if len(dir) == 0:
	sys.exit("\nThere are no files in the current working directory with the minimal name of *_L00*_R*_*.fastq\n")

for i in range(0, len(dir)):
	name = dir[i].split("_L00")

	if(len(name) == 2):
		files.append(name[0])

files = list(set(files))

#Check if all files exist
rs = ["1","2"]
ls = ["1","2","3","4"]
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
				used_files = used_files +  "\n\t\t\t\t\t" + files[i]

#Output variables
print "\n.fastq files recognised: 		" + str(len(dir))
print ".fastq files with unique index: 	" + str(len(files))
print "Unzip and use *.fastq.gz files:		" + params['-unzip']
print "Names of files of used:				 " + used_files
print "Used RefSeq file: 			" + params['-ref']
print "Used barcode file: 			" + params['-bar']
print "Output directory for the cat files:	" + params['-cat_out'] + "/"
print "Output directory for the bash files: 	" + params['-bash_out'] + "/"
print "Output directory for the mapping files: " + params['-map_out'] + "/"
print "Output directory for the counts files: 	" + params['-counts_out'] + "/"
# print "The used queue is: 			" + params['-qtype']
print "The number of used cores: 		" + params['-cores']
print "Submission to queue: 			" + params['-qsub']
print "Memory requested: 			" + params['-mem']
print "Time requested: 			" + params['-time']
print "Email: 					" + params['-email_fb']
print "zip .fastq files afterwards:		" + params['-zip']
print "\nImportant: Make sure to check the bash files!\n"

#Make unique filenames
for i in range(0,len(files)):
	nr = 1
	file_name = params['-bash_out'] + "/map_"+files[i]+"_"+str(nr)+".sh"

	while os.path.exists(file_name):
		split = file_name.split("_")[len(file_name.split("_"))-1]
		nr = int(split.split(".sh")[0]) + 1
		file_name = params['-bash_out'] + "/map_"+files[i]+"_"+str(nr)+".sh"

	bash_file_name.append(file_name)

#Make commands
for i in range(0,len(files)):
	bash_file = open(bash_file_name[i], "w")

	if params["-unzip"] == "yes":
		unzip = "gunzip " + os.getcwd() + "/" + files[i] + "_L00*_R*_*.fastq.gz" + "\n \n"
	else:
		unzip = ""

	cat_r1 = "cat " + os.getcwd() + "/" + files[i] + "_L00*_R1* > " + params["-cat_out"] + "/" + files[i] + "_R1_cat.fastq"
	cat_r2 = "cat " + os.getcwd() + "/" + files[i] + "_L00*_R2* > " + params["-cat_out"] + "/" + files[i] + "_R2_cat.fastq"
	# map = "do_mappings_strand.pl -r=" + params['-ref'] + " -f1=" + params["-cat_out"] + "/" + files[i] + "_R1_cat.fastq -f2=" + params["-cat_out"] + "/" + files[i] + "_R2_cat.fastq -out="+ files[i] + " -outdir=" + params['-map_out'] + " -t=" + params['-cores'] + " -uniq=1 -i=0 -cel=1 -fstr=1 -bar=" + params['-bar'] + " -rb > " + files[i] + ".log1 2> " + files[i] + ".log2"
	map = "do_mappings_strand.pl -r=" + params['-ref'] + " -f1=" + params["-cat_out"] + "/" + files[i] + "_R1_cat.fastq -f2=" + params["-cat_out"] + "/" + files[i] + "_R2_cat.fastq -out="+ files[i] + " -outdir=" + params['-map_out'] + " -t=" + params['-cores'] + " -uniq=1 -i=0 -cel=1 -fstr=1 -bar=" + params['-bar'] + " -rb > " + files[i] + ".log1 2> " + files[i] + ".log2"
	extract = "extract_counts_rb.pl -in=" + params['-map_out'] + "/" + files[i] + ".cout.csv -outc=" + params['-counts_out'] + "/" + files[i] + ".coutc.csv -outb=" + params['-counts_out'] + "/" + files[i] + ".coutb.csv -outt=" + params['-counts_out'] + "/" + files[i] + ".coutt.csv"

	if params['-zip'] == "yes":
		zip = "\n \ngzip " + os.getcwd() + "/" + files[i] + "*.fastq"
		if params['-cat_out'] != os.getcwd():
			zip = zip + "\n \ngzip " + params['-cat_out'] + "/" + files[i] + "*.fastq"
	else:
		zip = ""

	bash_text = head_bash + "\n \n" + unzip + cat_r1 + "\n \n" + cat_r2 + "\n \n" + map + "\n \n" + extract + zip
	bash_file.write(bash_text)
	bash_file.close()

#Submit to queue
if params['-qsub'] == 'yes':
	for i in range(0, len(files)):
		command = "qsub " + bash_file_name[i]
		os.system(command)

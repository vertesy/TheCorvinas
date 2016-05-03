### MapAndGo.py

This script is a wrapper over Dominic's wrapper around `BWA, samtools` and a transcript-count table generator.

MapAndGo is essentially generates an executable file that is list of bash commands, with specification ready so that it can be submitted to the queue on UMC's HPC.

## What it does

It creates a bash file (.sh) with 4 commands that are executed after each other.

1. `cat` Merge lanes 1-4 for read-pair 1
-  `cat` Merge lanes 1-4 for read-pair 2
-  `do_mappings_strand.pl` 
	- Align to reference using `BWA` → generating an `.sai` file.
	- Generate `.sam` file with `samtools`
-  ` extract_counts_rb.pl`

## Example output of MapAndGo.py

     #! /bin/bash
     #$ -cwd
     #$ -V
     #$ -l h_rt=24:00:00
     #$ -l h_vmem=10G
     #$ -M a.vertesy@hubrecht.eu
     #$ -m beas
     
     cat /hpc/hub_oudenaarden/Abel/sp1/play/lib1_L00*_R1* > /hpc/hub_oudenaarden/Abel/sp1/play/cat_files/lib1_R1_cat.fastq
     
     cat /hpc/hub_oudenaarden/Abel/sp1/play/lib1_L00*_R2* > /hpc/hub_oudenaarden/Abel/sp1/play/cat_files/lib1_R2_cat.fastq
     
     do_mappings_strand.pl -r=/hpc/hub_oudenaarden/gene_models/mouse_gene_models/mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked.fa -f1=/hpc/hub_oudenaarden/Abel/sp1/play/cat_files/lib1_R1_cat.fastq -f2=/hpc/hub_oudenaarden/Abel/sp1/play/cat_files/lib1_R2_cat.fastq -out=lib1 -outdir=/hpc/hub_oudenaarden/Abel/sp1/play/map_files -t=8 -uniq=1 -i=0 -cel=1 -fstr=1 -bar=/hpc/hub_oudenaarden/data/cel-seq_barcodes.csv -rb > lib1.log1 2> lib1.log2
     
     extract_counts_rb.pl -in=/hpc/hub_oudenaarden/Abel/sp1/play/map_files/lib1.cout.csv -outc=/hpc/hub_oudenaarden/Abel/sp1/play/count_files/lib1.coutc.csv -outb=/hpc/hub_oudenaarden/Abel/sp1/play/count_files/lib1.coutb.csv -outt=/hpc/hub_oudenaarden/Abel/sp1/play/count_files/lib1.coutt.csv
     


### Overview of the parameters for Dominic's script.


     do_mappings_strand.pl
      -r=/hpc/hub_oudenaarden/gene_models/mouse_gene_models/mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked.fa
      -f1=/hpc/hub_oudenaarden/Abel/sp1/play/cat_files/lib1_R1_cat.fastq
      -f2=/hpc/hub_oudenaarden/Abel/sp1/play/cat_files/lib1_R2_cat.fastq
      -out=lib1
      -outdir=/hpc/hub_oudenaarden/Abel/sp1/play/map_files
      -t=8
      -uniq=1
      -i=0
      -cel=1
      -fstr=1
      -bar=/hpc/hub_oudenaarden/data/cel-seq_barcodes.csv
      -rb > lib1.log1 2> lib1.log2
     
     extract_counts_rb.pl
      -in=/hpc/hub_oudenaarden/Abel/sp1/play/map_files/lib1.cout.csv
      -outc=/hpc/hub_oudenaarden/Abel/sp1/play/count_files/lib1.coutc.csv
      -outb=/hpc/hub_oudenaarden/Abel/sp1/play/count_files/lib1.coutb.csv
      -outt=/hpc/hub_oudenaarden/Abel/sp1/play/count_files/lib1.coutt.csv
      
## Usage
1. Put the raw sequence data in a folder, put this script to folder where you store your scripts (e.g.: user/bin), and make it executable `chmod +X path/to/MapAndGo.py`
- The files should contain the minimum name of *_L00*_R*_001.fastq (that is Illumina NextSeq's default output: R for read 1 or 2, L for one of the 4 lanes)
- Test it by running in bash: `python path/to/MapAndGo.py -help`
- Run it with your data: python `path/to/MapAndGo.py MapAndGo.py -ref=human -bar=cel-seq_barcodes.csv -bash_out=bash_files -cat_out=cat_files  -map_out=map_files -counts_out=count_files -email=x.y@hubrecht.eu`
- Original author: Thom de Hoog, van Oudenaarden group, 02-03-2015
- Modified and maintained by Abel Vertesy, van Oudenaarden group, 28-04-2016



## Final outputs of the mapping

`coutc.csv`	Raw read counts  
`coutb.csv` UMI-corrected read counts ()  `coutt.csv` Estimated original transcript counts, based on how many .

Use `coutt` for downstream analysis.
`coutc / coutb` = oversequenceing

look at the `.sout` for mapping statistics.
### MapAndGo.py


MapAndGo is a makefile that generates an executable file that is list of bash commands, with specification ready so that it can be submitted to the queue on UMC's HPC.

## What it does

It creates a bash file (.sh) with 4 commands that are executed after each other.

1. Merge lanes 1-4 for read 1
-  Merge lanes 1-4 for read 2
-  

## Example output

```

     #! /bin/bash
     #$ -cwd
     #$ -V
     #$ -l h_rt=24:00:00
     #$ -l h_vmem=10G
     #$ -M a.vertesy@hubrecht.eu
     #$ -m beas
     
     # cat /hpc/hub_oudenaarden/Abel/sperm/play/SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12_L00*_R1* > /hpc/hub_oudenaarden/Abel/sperm/play/cat_files/SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12_R1_cat.fastq
     
     # cat /hpc/hub_oudenaarden/Abel/sperm/play/SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12_L00*_R2* > /hpc/hub_oudenaarden/Abel/sperm/play/cat_files/SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12_R2_cat.fastq
     
     do_mappings_strand.pl -r=/hpc/hub_oudenaarden/gene_models/mouse_gene_models/mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked.fa -f1=/hpc/hub_oudenaarden/Abel/sperm/play/cat_files/SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12_R1_cat.fastq -f2=/hpc/hub_oudenaarden/Abel/sperm/play/cat_files/SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12_R2_cat.fastq -out=SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12 -outdir=/hpc/hub_oudenaarden/Abel/sperm/play/map_files -t=8 -uniq=1 -i=0 -cel=1 -fstr=1 -bar=/hpc/hub_oudenaarden/data/cel-seq_barcodes.csv -rb > SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12.log1 2> SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12.log2
     
     extract_counts_rb.pl -in=/hpc/hub_oudenaarden/Abel/sperm/play/map_files/SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12.cout.csv -outc=/hpc/hub_oudenaarden/Abel/sperm/play/count_files/SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12.coutc.csv -outb=/hpc/hub_oudenaarden/Abel/sperm/play/count_files/SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12.coutb.csv -outt=/hpc/hub_oudenaarden/Abel/sperm/play/count_files/SvdB-6e1h-pl1-lib1-Diss1h_AHVCTTBGXX_S12.coutt.csv```


## Usage
- Put the raw sequence data in a folder, put this script to folder where you store your scripts (e.g.: user/bin), and make it executable `chmod +X path/to/MapAndGo.py`
- The files should contain the minimum name of *_L00*_R*_001.fastq (that is Illumina NextSeq's default output: R for read 1 or 2, L for one of the 4 lanes)
- Test it by running in bash: `python path/to/MapAndGo.py -help`
- Run it with your data: python `path/to/MapAndGo.py MapAndGo.py -ref=human -bar=cel-seq_barcodes.csv -bash_out=bash_files -cat_out=cat_files  -map_out=map_files -counts_out=count_files -email=x.y@hubrecht.eu`
- Original author: Thom de Hoog, van Oudenaarden group, 02-03-2015
- Modified and maintained by Abel Vertesy, van Oudenaarden group, 28-04-2016

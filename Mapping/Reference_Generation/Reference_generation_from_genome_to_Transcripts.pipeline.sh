#!/usr/bin/env bash
'''
### Reference_generation_from_genome_to_Transcripts.pipeline.sh
- Clean up and reformat the genome by *00.genome_fa_reformatter_1chr_per_line.py*
- Convert UCSC to GTF by *01.ucsc2gtf.pl*
- Merge of all exons in all isofrom by *02.merge_isoforms_gtf.pl*
- Convert the genome from normal [multifasta] genome.fa > genome.refomat.fa with: 1 chr / 1 line
- Extract sequences from genome.fa by *03.gtf2fa.pl*
- Add ERCC SpikeIn sequences
- Optional: Replace poly-A stretches by poly N-s by *04.mask_polyA.pl*
'''

genome_path='/Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/genomes/'

scripts_path='/Users/abelvertesy/Github_repos/TheCorvinas/bash/Reference_Generation/'
cd $scripts_path

# 2.-3. Reformat and clean the genome
$scripts_path'00.genome_fa_reformatter_1chr_per_line.py'

# Transform ucsc it into gtf format
$scripts_path'01.ucsc2gtf.pl' -in='hg38.RefSeq.ucsc' -out='02.hg38.Refseq.gtf' -m='02.hg38.Refseq.tsv'
'maybe new line after last'
# gtf 1 line per exon, multiple isoforms for the same gene follow each other: each marked by 'transcript' in col3

#  Merge of all exons in all isofrom
$scripts_path'02.merge_isoforms_gtf.pl' -in='02.hg38.Refseq.gtf' -cl='02.hg38.Refseq.tsv' -out='03.hg38.Refseq_genes.gtf'
# -cl= key value: gene name >> all isoforms
# -out= no more different isoforms


# extract sequences from genome.fa
### genome has to be in 1 line / chr format
$scripts_path'03.gtf2fa.pl' -in=$scripts_path'03.hg38.Refseq_genes.gtf' -ref=$genome_path'/hg38.reformat.clean.fa' > 'hg38.REF.Refseq_genes.fa'


# Check if the NR of genes and the names are correct
grep '^>' 'hg38.REF.Refseq_genes.fa' | wc -l
grep '^>' 'hg38.REF.Refseq_genes.fa'
head 'hg38.REF.Refseq_genes.fa'


# Add ERCC SpikeIn sequences
cd $genome_path
ERCC_fnp='/Users/abelvertesy/Github_repos/TheCorvinas/bash/Reference_Generation/05.5.Reference_generation/ERCC/ERCC92_reformat.fa'
cat 'hg38.REF.Refseq_genes.fa' $ERCC_fnp > 'hg38.REF.Refseq_genes_ERCC.fa'


# Replace poly-A stretches by poly N-s [BWA counts as mismatch]
$scripts_path'04.mask_polyA.pl' -in='hg38.REF.Refseq_genes_ERCC.fa' -n=15 -out='hg38.REF.Refseq_genes_ERCC_polyN.fa'

# 'Make an index with bwa in the next script'



# 'Reference GENCODE -----------'
# $scripts_path'03.gtf2fa.pl' -in='/Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/Data_analysis/UCSC_GTF_tables/hg38/Gencode/03.hg38_GencodeV19_genes.gtf' -ref='/Users/abelvertesy/Downloads/00_zacc/hg38_reformat.fa' > '/Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/Data_analysis/UCSC_GTF_tables/hg38/Gencode/hg38.REF.GencodeV19.fa'
# cd /Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/Data_analysis/UCSC_GTF_tables/hg38/Gencode/
# cat 'hg38.REF.GencodeV19.fa' $ERCC_fnp > 'hg38.REF.GencodeV19_genes_ERCC.fa'
# $scripts_path'04.mask_polyA.pl' -in='hg38.REF.GencodeV19_genes_ERCC.fa' -n=15 -out='hg38.REF.GencodeV19_genes_ERCC_polyN.fa'
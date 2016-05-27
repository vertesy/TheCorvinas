#!/usr/bin/env bash
'''
### Reference_generation_from_genome_to_Transcripts.pipeline.sh
The numbers correspond to the steps in the wiki.
- 2.-3. Clean up and reformat the genome by *00.genome_fa_reformatter_1chr_per_line.py*
- 5. Convert UCSC to GTF by *01.ucsc2gtf.pl*
- 6. Merge of all exons in all isofrom by *02.merge_isoforms_gtf.pl*
- 7. Extract sequences from genome.fa by *03.gtf2fa.pl*
	- Check if the NR of genes and the names are correct
- 8. Add ERCC SpikeIn sequences
- 9. Optional: Replace poly-A stretches by poly N-s by *04.mask_polyA.pl*
'''

genome_path='/Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/genomes/'

scripts_path='/Users/abelvertesy/Github_repos/TheCorvinas/bash/Reference_Generation/'
cd $scripts_path

# 2.-3. Reformat and clean the genome
$scripts_path'00.genome_fa_reformatter_1chr_per_line.py'

# 5. Transform ucsc it into gtf format
$scripts_path'01.ucsc2gtf.pl' -in='hg38.RefSeq.ucsc' -out='02.hg38.Refseq.gtf' -m='02.hg38.Refseq.tsv'
'maybe new line after last'
# gtf 1 line per exon, multiple isoforms for the same gene follow each other: each marked by 'transcript' in col3

# 6. Merge of all exons in all isofrom
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


# 8. Add ERCC SpikeIn sequences
cd $genome_path
ERCC_fnp='/Users/abelvertesy/Github_repos/TheCorvinas/bash/Reference_Generation/05.5.Reference_generation/ERCC/ERCC92_reformat.fa'
cat 'hg38.REF.Refseq_genes.fa' $ERCC_fnp > 'hg38.REF.Refseq_genes_ERCC.fa'


# 9. Replace poly-A stretches by poly N-s [BWA counts as mismatch]
$scripts_path'04.mask_polyA.pl' -in='hg38.REF.Refseq_genes_ERCC.fa' -n=15 -out='hg38.REF.Refseq_genes_ERCC_polyN.fa'

# 'Make an index with bwa in the next script'


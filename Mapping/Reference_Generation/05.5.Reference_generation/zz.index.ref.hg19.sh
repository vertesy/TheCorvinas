echo "bwa index -a bwtsw /hpc/hub_oudenaarden/Abel/genomes/hg_19/hg19_reformat.fa " | qsub -q medium -M a.vertesy@hubrecht.eu -m beas -pe threaded 6
echo "bwa index -a bwtsw /hpc/hub_oudenaarden/Abel/genomes/hg_19/GencodeV19/hg19.REF.GencodeV19.fa " | qsub -q medium -M a.vertesy@hubrecht.eu -m beas -pe threaded 6
echo "bwa index -a bwtsw /hpc/hub_oudenaarden/Abel/genomes/hg_19/GencodeV19/hg19.REF.GencodeV19_genes_ERCC.fa " | qsub -q medium -M a.vertesy@hubrecht.eu -m beas -pe threaded 6
echo "bwa index -a bwtsw /hpc/hub_oudenaarden/Abel/genomes/hg_19/GencodeV19/hg19.REF.GencodeV19_genes_ERCC_polyN.fa " | qsub -q medium -M a.vertesy@hubrecht.eu -m beas -pe threaded 6
echo "bwa index -a bwtsw /hpc/hub_oudenaarden/Abel/genomes/hg_19/RefSeq/hg19.REF.Refseq_genes.fa " | qsub -q medium -M a.vertesy@hubrecht.eu -m beas -pe threaded 6
echo "bwa index -a bwtsw /hpc/hub_oudenaarden/Abel/genomes/hg_19/RefSeq/hg19.REF.Refseq_genes_ERCC.fa " | qsub -q medium -M a.vertesy@hubrecht.eu -m beas -pe threaded 6
echo "bwa index -a bwtsw /hpc/hub_oudenaarden/Abel/genomes/hg_19/RefSeq/hg19.REF.Refseq_genes_ERCC_polyN.fa " | qsub -q medium -M a.vertesy@hubrecht.eu -m beas -pe threaded 6
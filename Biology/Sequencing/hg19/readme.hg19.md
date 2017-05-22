## Reference_Stats 23868 transcripts in hg19

`/hpc/hub_oudenaarden/gene_models/human_gene_models/hg19_mito/hg19_RefSeq_genes_clean_ERCC92_polyA_10_masked_Mito.fa`

### Count bases per transcript

`cd /hpc/hub_oudenaarden/gene_models/human_gene_models/hg19_mito/`
`awk '{ print length($0); }' hg19_RefSeq_genes_clean_ERCC92_polyA_10_masked_Mito.fa > Length.txt`

`grep ">" hg19_RefSeq_genes_clean_ERCC92_polyA_10_masked_Mito.fa > ListOfTranscripts.mm10.txt`

### Process

```

TranscriptLength.hg19 = read.simple.tsv.named.vector("/Users/abelvertesy/Github_repos/TheCorvinas/Mapping/Reference_Stats/hg19/TranscriptLength.hg19.tsv")
OutDir = "/Users/abelvertesy/Github_repos/TheCorvinas/Mapping/Reference_Stats/hg19/"
whist(TranscriptLength.hg19, breaks = 50)

log10_TranscriptLength.hg19 = log10(TranscriptLength.hg19)
whist(log10_TranscriptLength.hg19, breaks = 50, vline = log10(c(1977,2234)))
round(10^(summary(log10_TranscriptLength.hg19)))

```
#### Summary

    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     21    1334    2449    1986    4064  117220 
     
![](log10_TranscriptLength.hg19.hist.pdf)
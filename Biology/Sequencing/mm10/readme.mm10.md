## Reference_Stats


### Count bases per transcript

`awk '{ print length($0); }' mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked_eGFP_Mito.fa > tmp.txt`

`grep ">" mm10_RefSeq_genes_clean_ERCC92_polyA_10_masked_eGFP_Mito.fa > ListOfTranscripts.mm10.txt`


### Process

```
TranscriptLength.mm10 = read.simple.tsv.named.vector("/Users/abelvertesy/Github_repos/TheCorvinas/Mapping/Reference_Stats/mm10/TranscriptLength.mm10.tsv")
OutDir = "/Users/abelvertesy/Github_repos/TheCorvinas/Mapping/Reference_Stats/mm10/"
whist(TranscriptLength.mm10, breaks = 50)

log10_TranscriptLength.mm10 = log10(TranscriptLength.mm10)
whist(log10_TranscriptLength.mm10, breaks = 50)
round(10^(summary(log10_TranscriptLength.mm10)))

```
#### Summary

    Min. 1st Qu.  Median    Mean    3rd Qu.    Max. 
     21    1178    2234       1977    3673      104713 
     
![](log10_TranscriptLength.mm10.hist.pdf)


### Freq

#### Quantiles  

     5%  95%   
     0.19 0.32 

#### Code

```
BaseFrequencies.mm10 = read.simple.tsv("/Users/abelvertesy/Github_repos/TheCorvinas/Mapping/Reference_Stats/mm10/BaseFrequencies.mm10.tsv")
OutDir = "/Users/abelvertesy/Github_repos/TheCorvinas/Mapping/Reference_Stats/mm10/"

pname=  "Basedistr"
pdfA4plot_on(pname = pname, rows = 3, cols = 2)

Bases = c("A","C","G","T" )
  for (L in 1:l(Bases)) {
    base=Bases[L]; print(base)

    freqz = BaseFrequencies.mm10[,base]
    linez = iround(c(mean(freqz), median(freqz)))
    whist(freqz, breaks = 50, vline = linez, col=L+1, savefile = F, main = base)
    print (quantile(freqz, probs = c(.05,.95)))
  }

Distribution = splitByCol(BaseFrequencies.mm10, f = Bases)
wvioplot_list(Distribution, savefile = F)
pdfA4plot_off()
```

![](Basedistr.pdf)
####  Stats

| | A 	| C 	| G 	| T  |
| ---| ---| ---| ---| --- |
| quantile_5% 	| 19 % 	| 19 % 	| 19 % 	| 19 %  |
| mean 	| 25.4 % 	| 24.6 % 	| 24.9 % 	| 25 %  |
| quantile_95 	| 32 % 	| 31 % 	| 31 % 	| 32 %  |

## Filtering done by RaceID

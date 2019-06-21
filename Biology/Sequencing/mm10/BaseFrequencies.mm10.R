
BaseFrequencies.mm10 = read.simple.tsv("~/GitHub/TheCorvinas/Mapping/Reference_Stats/mm10/BaseFrequencies.mm10.tsv")
OutDir = "~/GitHub/TheCorvinas/Mapping/Reference_Stats/mm10/"

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


# -----

path_of_report = "~/GitHub/TheCorvinas/Mapping/Reference_Stats/mm10/readme.mm10.md"

Stats = rbind(
  "quantile_5%" = percentage_formatter(apply(BaseFrequencies.mm10[,1:4], 2,quantile, probs = .05)),
  "mean" = percentage_formatter(apply(BaseFrequencies.mm10[,1:4], 2,mean)),
  "quantile_95" = percentage_formatter(apply(BaseFrequencies.mm10[,1:4], 2,quantile, probs = .95))
)
colnames(Stats)=Bases[1:4]
MarkDown_Table_writer_DF_RowColNames(Stats)

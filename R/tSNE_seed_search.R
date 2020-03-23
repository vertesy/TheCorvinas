
# source("~/GitHub/TheCorvinas/R/tSNE_seed_search.R")

# try (source ('~/GitHub/CodeAndRoll/CodeAndRoll.R'),silent= F)

OutDir = kollapse("~/Google_Drive/Spermatogenesis_Data/sp1_analysis_new/tSNE_seed_Search_",signif(runif(1)))
setup_MarkdownReports(OutDir = OutDir)

seedz =  c(211, -211, iround (runif(50, min = -10000, max = 10000) ));seedz
write.simple.vec(seedz)

s=1
for (s in 1:seedz) {
  sx = seedz[s]
  sc <- comptsne2(sc,rseed=sx)
  plottsne2(sc, final=F)
  wplot_save_this(plotname = kollapse("tSNE_seed_",sx))
}


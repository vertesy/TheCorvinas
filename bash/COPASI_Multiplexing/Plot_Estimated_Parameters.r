######################################################################
# Plot_Estimated_Parameters.r
######################################################################
# source ("/Users/abelvertesy/Dokumentumok/Tanulas/HD/HD_04/Thesis/Presentations_Thesis/After_Thesis_submission/2016_Modeling/Q6/Plot_Estimated_Parameters.r")

# This scripts read in a Fits*Parameters dataframe and plots histograms and pairwise correlations among them. The file should have row and column names.

# Functions ------------------------
require(MarkdownReports)
source ('/Users/abelvertesy/Github_repos/TheCorvinas/R/CodeAndRoll.R')

cormethod = "spearman"
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, method = cormethod) {
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y, method = method)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}


# Setup ------------------------
OutDir = "/Users/abelvertesy/Dokumentumok/Tanulas/HD/HD_04/Thesis/Presentations_Thesis/After_Thesis_submission/2016_Modeling/q10/"
setup_MarkdownReports(OutDir = OutDir, scriptname = "Plot_Estimated_Parameters.r", append = F)

# Read in the files ------------------------
infile = "/Users/abelvertesy/Dokumentumok/Tanulas/HD/HD_04/Thesis/Presentations_Thesis/After_Thesis_submission/2016_Modeling/q10/FitResults_q10_BestFits.tsv"
# infile = "/Users/abelvertesy/Dokumentumok/Tanulas/HD/HD_04/Thesis/Presentations_Thesis/After_Thesis_submission/2016_Modeling/Q6/FitResults_best.tsv"

Estimates = read.simple.tsv(infile)
NrParams = ncol(Estimates)
NrFits = nrow(Estimates)
Params = colnames(Estimates)

pdfA4plot_on(pname = "ParameterHistograms")
for (p in Params) {
  print(p)
  hist(Estimates[,p], main = p, breaks =NrFits/2, col="seagreen")
  # whist_dfCol(Estimates, colName = p, breaks =50)
}
pdfA4plot_off()

pdfA4plot_on(pname = "ParameterHistograms_log10")
for (p in Params) {
  print(p)
  hist(log10(Estimates[,p]), main = p, breaks =NrFits/2, col="seagreen")
  # whist_dfCol(Estimates, colName = p, breaks =50)
}
pdfA4plot_off()

# Correlation plots ------------------------------------------------------------------------------------------------------------

par("pch" =18)
par("col" =rgb(0,0,0,.3))
ccc=(val2col(-Estimates$ObjValue))
  
cormethod = "spearman"
pname ="Pairwise_correlation_of_parameters"
pairs(Estimates, lower.panel=panel.smooth, upper.panel=panel.cor, main = pname) # , pch =".", cex.labels = .5
wplot_save_this(plotname = pname, w=15, h=15, mdlink = T)
help(pairs)

cormethod = "pearson"
pname ="Pairwise_correlation_of_parameters_log10_pearson"
pairs(x = log10(Estimates), lower.panel=panel.smooth, upper.panel=panel.cor, main = pname) # , cex.labels = .5
wplot_save_this(plotname = pname, w=15, h=15, mdlink = T)

cormethod = "spearman"
pname ="Pairwise_correlation_of_parameters_log10_spearman"
pairs(x = log10(Estimates), lower.panel=panel.smooth, upper.panel=panel.cor, main = pname) # , cex.labels = .5
wplot_save_this(plotname = pname, w=15, h=15, mdlink = T)

par("col" ="black")

# Color by goodness of fit (smallest objective vales are red) ------------------------------------------------------

cormethod = "spearman"
pname ="Pairwise_correlation_of_parameters_RedHot"
pairs(Estimates, col=ccc, lower.panel=panel.smooth, main = pname) # , pch =".", cex.labels = .5
wplot_save_this(plotname = pname, w=15, h=15, mdlink = T)
help(pairs)

cormethod = "pearson"
pname ="Pairwise_correlation_of_parameters_log10_pearson_RedHot"
pairs(x = log10(Estimates), col=ccc, lower.panel=panel.smooth, main = pname) # , cex.labels = .5
wplot_save_this(plotname = pname, w=15, h=15, mdlink = T)

cormethod = "spearman"
pname ="Pairwise_correlation_of_parameters_log10_spearman_RedHot"
pairs(x = log10(Estimates), col=ccc, lower.panel=panel.smooth, main = pname) # , cex.labels = .5
wplot_save_this(plotname = pname, w=15, h=15, mdlink = T)




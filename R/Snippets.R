######################################################################
# Description
######################################################################
# source ('')
# rm(list=ls(all.names = TRUE)); try(dev.off(), silent = T)

i=1
# Combinations.all.vs.all.R  ------------------------

CN =1:3; i=1


# With itself, and both A, B and B,A
PariwiseCombinations = data.table::CJ(CN,CN); idim(PariwiseCombinations)
for (i in 1:NROW(PariwiseCombinations) ) {  CPL = as.character(PariwiseCombinations[i,]) }

# WithOUT itself, and only A, B NO B,A pairs
PariwiseCombinations = t(combinat::combn(CN, m = 2)); idim(PariwiseCombinations)
for (i in 1:NROW(PariwiseCombinations) ) {  CPL = PariwiseCombinations[i,] }




# Hist w normal distribution ------------------------
set.seed(21768134)
ccc = sample(richColors(ncol(Protz.NullNorm)+n)[-(1:n)] )

pdfA4plot_on(pname = "Histograms.Normdata", rows = 4, cols=3)
i=2
for (i in 1:ncol(Protz.NullNorm) ) {
  P =colnames(Protz.NullNorm)[i]
  g = na.omit.strip(Protz.NullNorm[,P])
  hist(g, xlab = P, breaks = 20, main = P, col=ccc[i], sub="Matched (mean, var) normal distribution overlaid")
  curve(dnorm(x, mean=mean(g), sd=sd(g)), lwd=2, lty=2,add=TRUE, yaxt="n")
} #for
pdfA4plot_off()


# Metadata ------------------------

# Parameters ------------------------

# Read In ------------------------

# QC ------------------------

# Plot ------------------------



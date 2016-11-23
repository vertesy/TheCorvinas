######################################################################
# SurfacePatterns
######################################################################
# source ('/Users/abelvertesy/Github_repos/TheCorvinas/R/SurfacePatterns.R')
# rm(list=ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
try (source ('/Users/abelvertesy/Github_repos/TheCorvinas/R/CodeAndRoll.R'),silent= F)
try (source("/Users/abelvertesy/Github_repos/Spermatogenesis/Mapping/functions_sp.r") , silent= F)
source("/Users/abelvertesy/Github_repos/hmC/RaceID/RaceID2_class_ext.2016.05.10.R") ## load class definition and functions
source ('/Users/abelvertesy/Github_repos/TheCorvinas/R/DatabaseLinke.r')

# Setup ------------------------

setup_MarkdownReports(OutDir = "/Users/abelvertesy/Google_Drive/Zacc/SurfacePatterns", scriptname = "SurfacePatterns")
# Metadata ------------------------

CexSizes = c(14, 15)
Xpoz = seq(from=15, to = 90, by = 15)
Ypoz = seq(from=.15, to = .9, by = .15)

# Parameters ------------------------

x=cbind(Xpoz, Ypoz)
plot(x, cex=CexSizes)
# Read In ------------------------


# QC ------------------------

# Plot ------------------------


######################################################################
# Install R packages
######################################################################

install.packages("tidyverse")
install.packages("Seurat")

install.packages("clipr")
install.packages("tictoc")
install.packages("doMC")
install.packages("tictoc")
install.packages("pheatmap")
install.packages("biomaRt")
# If you don't have it
install.packages("openxlsx")
install.packages("plotrix")
install.packages("princurve")
install.packages("HGNChelper")
install.packages("R.utils")
# install.packages("limma")
install.packages("SoupX")

install.packages("colorRamps")
install.packages("data.table")
install.packages("readr")
install.packages("sm")
install.packages("gplots")
install.packages("gtools")
install.packages("RColorBrewer")
install.packages("vioplot")
install.packages("VennDiagram")


# BioConductor---------------------------------------------------------------------------------
install.packages("BiocManager")
BiocManager::install("schex")
BiocManager::install("biomaRt")
BiocManager::install("STRINGdb")
BiocManager::install("sparseMatrixStats")
BiocManager::install("limma")

# BiocManager::install("GO.db")
# BiocManager::install("WGCNA")

update.packages(ask = F)

# Github---------------------------------------------------------------------------------
install.packages("devtools")
require("devtools")

devtools::install_github(repo = "jalvesaq/colorout", upgrade = F)
"maybe requires X11"




remotes::install_github(repo = "vertesy/Stringendo", upgrade = F)
remotes::install_github(repo = "vertesy/CodeAndRoll2", upgrade = F)
remotes::install_github(repo = "vertesy/ReadWriter", upgrade = F)
remotes::install_github(repo = "vertesy/MarkdownHelpers", upgrade = F)
remotes::install_github(repo = "vertesy/Markdownreports", upgrade = F)
remotes::install_github(repo = "vertesy/ggExpress", upgrade = F)
remotes::install_github(repo = "vertesy/Seurat.utils", upgrade = F)

remotes::install_github(repo = "vertesy/UVI.tools", upgrade = F)
remotes::install_github(repo = "vertesy/Connectome.tools", upgrade = F)

# Less important ones
remotes::install_github(repo = "vertesy/DataInCode", upgrade = F)
remotes::install_github(repo = "vertesy/DatabaseLinke.R", upgrade = F)




## Else   -------------------------------------------------------------------------------------------------

#
# # MULTI-seq
# devtools::install_github('chris-mcginnis-ucsf/MULTI-seq', force = TRUE)
# require("deMULTIplex") # This is MULTI-seq
#
# BiocManager::install("ShortRead")
#
# BiocManager::install("DropletUtils")
#
#
#
#
# library(devtools)
# install_github("immunogenomics/harmony")
# devtools::install_github("immunogenomics/presto")
#
# devtools::install_local("/Users/abel.vertesy/Downloads/harmony")

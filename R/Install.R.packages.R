######################################################################
# A collection of custom R functions
######################################################################



install.packages("tidyverse")
install.packages("Seurat")
install.packages("clipr")
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

# BioConductor---------------------------------------------------------------------------------
install.packages("BiocManager")
BiocManager::install("schex")
BiocManager::install("biomaRt")
BiocManager::install("STRINGdb")
BiocManager::install("sparseMatrixStats")
BiocManager::install("limma")

update.packages(ask = F)

# Github---------------------------------------------------------------------------------
install.packages("devtools")
require("devtools")

devtools::install_github(repo = "jalvesaq/colorout", upgrade = F)
"maybe requires X11"


devtools::install_github(repo = "vertesy/Stringendo", upgrade = F)
devtools::install_github(repo = "vertesy/CodeAndRoll2", upgrade = F)
devtools::install_github(repo = "vertesy/ReadWriter", upgrade = F)
devtools::install_github(repo = "vertesy/MarkdownHelpers", upgrade = F)
devtools::install_github(repo = "vertesy/ggExpress", upgrade = F)
devtools::install_github(repo = "vertesy/Markdownreports", upgrade = F)
devtools::install_github(repo = "vertesy/Seurat.utils", upgrade = F)

# Less important ones
devtools::install_github(repo = "vertesy/DataInCode", upgrade = F)
devtools::install_github(repo = "vertesy/DatabaseLinke.R", upgrade = F)




## Setup   -------------------------------------------------------------------------------------------------
install.packages("BiocManager")
install.packages("tidyverse")
install.packages("devtools")
# install.packages("bookdown")
install.packages("sm")
install.packages("gplots")
install.packages("gtools")

install.packages("RColorBrewer")
install.packages("colorRamps")
install.packages("vioplot")
install.packages("VennDiagram")
# install.packages("rmarkdown")
# install.packages("commonmark")
install.packages("clipr")
install.packages("data.table")
install.packages("pheatmap")
# install.packages("revealjs")
install.packages("readr")
install.packages("Seurat")
# install.packages("robustbase")
# devtools::install_github(repo = "vertesy/MarkdownReportsDev")
# devtools::install_github("csgillespie/roxygen2Comment")



require("devtools")
devtools::install_github(repo = "vertesy/Stringendo", upgrade = F)
devtools::install_github(repo = "vertesy/CodeAndRoll2", upgrade = F)
devtools::install_github(repo = "vertesy/ReadWriter", upgrade = F)
devtools::install_github(repo = "vertesy/MarkdownHelpers", upgrade = F)
devtools::install_github(repo = "vertesy/ggExpress", upgrade = F)
devtools::install_github(repo = "vertesy/Markdownreports", upgrade = F)
devtools::install_github(repo = "vertesy/Seurat.utils", upgrade = F)
devtools::install_github(repo = "vertesy/DataInCode", upgrade = F)
devtools::install_github(repo = "vertesy/DatabaseLinke.R", upgrade = F)


######################################################################
# Install R packages
######################################################################

if(!require("tidyverse")) install.packages("tidyverse")
if(!require("Seurat")) install.packages("Seurat")

if(!require("clipr")) install.packages("clipr")
if(!require("tictoc")) install.packages("tictoc")
if(!require("doMC")) install.packages("doMC")
if(!require("tictoc")) install.packages("tictoc")
if(!require("pheatmap")) install.packages("pheatmap")
if(!require("biomaRt")) install.packages("biomaRt")
# If you don't have it
if(!require("openxlsx")) install.packages("openxlsx")
if(!require("plotrix")) install.packages("plotrix")
if(!require("princurve")) install.packages("princurve")
if(!require("HGNChelper")) install.packages("HGNChelper")
if(!require("R.utils")) install.packages("R.utils")
if(!require("SoupX")) install.packages("SoupX")

if(!require("colorRamps")) install.packages("colorRamps")
if(!require("data.table")) install.packages("data.table")
if(!require("readr")) install.packages("readr")
if(!require("sm")) install.packages("sm")
if(!require("gplots")) install.packages("gplots")
if(!require("gtools")) install.packages("gtools")
if(!require("RColorBrewer")) install.packages("RColorBrewer")
if(!require("vioplot")) install.packages("vioplot")
if(!require("VennDiagram")) install.packages("VennDiagram")
if(!require("devtools")) install.packages("devtools")



# BioConductor---------------------------------------------------------------------------------

'Error: Bioconductor version '3.16' requires R version '4.2'; use'
'`BiocManager::install(version = '3.17')` with R version 4.3; see'

# install.packages("BiocManager")
if( !require("biomaRt") ) BiocManager::install("biomaRt")
if( !require("sparseMatrixStats") ) BiocManager::install("sparseMatrixStats")
if( !require("limma") ) BiocManager::install("limma")
if( !require("scales") ) BiocManager::install("scales")
if( !require("vroom") ) BiocManager::install("vroom")

# if( !require("ggcorrplot") ) BiocManager::install("ggcorrplot")
if( !require("schex") ) BiocManager::install("schex")
if( !require("STRINGdb") ) BiocManager::install("STRINGdb")
if( !require("EnhancedVolcano") ) BiocManager::install("EnhancedVolcano")
if( !require("ggcorrplot") ) BiocManager::install("ggcorrplot")


# BiocManager::install("GO.db")
# BiocManager::install("WGCNA")


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
require(ggExpress)
remotes::install_github(repo = "vertesy/Seurat.utils", upgrade = F)
require(Seurat.utils)

remotes::install_github(repo = "vertesy/UVI.tools", upgrade = F, auth_token = "")
remotes::install_github(repo = "vertesy/Connectome.tools", upgrade = F, auth_token = "")

# Less important ones
remotes::install_github(repo = "vertesy/DataInCode", upgrade = F)
remotes::install_github(repo = "vertesy/DatabaseLinke.R", upgrade = F)



## Else   -------------------------------------------------------------------------------------------------
update.packages(ask = F)




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


require (biomaRt)
"http://bioconductor.org/packages/release/bioc/html/biomaRt.html"

# #### Set up biomart ####
ensembl <- useMart("ensembl")

"Mart object ocnnected to ENS database"
listMarts (); "to find out waht marts you can use"

"select a dataset form the mart object"
listDatasets () "to find out waht ds you can use"

ensembl_mouse <- useDataset(ensembl, dataset = "mmusculus_gene_ensembl")
ensembl_human <- useDataset(ensembl, dataset = "hsapiens_gene_ensembl")

#### Select info from dataset homology (only genes included that have a mgi-symbol and a hgnc-symbol) ####
# get homology list
hom <- getBM(attributes=c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"),mart=ensembl_mouse)
getGene(); "is another func"


hom <- hom[hom[,2] != "",]
names(hom) <- c("mouse_ensembl_id","human_ensembl_id")

"Which version??? figure out"

"attributes and filters are the 2 main things"
listAttributes(ensembl_human) # 1300 attribs!
listFilters(ensembl_human)

# get mouse gene symbols
mouse_symbol <- getBM(attributes=c("ensembl_gene_id","mgi_symbol","external_gene_name"),filters = c("ensembl_gene_id"), values = hom[,1], mart=ensembl_mouse)
names(mouse_symbol) <- c("mouse_ensembl_id","mgi_symbol","other_symbol")

"Maybe better  to get ensembl ID first"
"you cannot flter on more that 3 external references, like HGNC symbols"
# get human gene symbols
human_symbol <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol","external_gene_name"),filters = c("ensembl_gene_id"), values = hom[,2], mart=ensembl_human)
names(human_symbol) <- c("human_ensembl_id","hngc_symbol","other_symbol")


getBM( mart = ensembl_human, attributes = "ensembl_gene_id", filters = "hgnc_symbol", values = c("INS"))
"@ mart you provide the dataset, not the Mart"

# remove duplicated gene IDs (sometime two names for one gene id = take only name in external resource)
mouse_symbol <- mouse_symbol[mouse_symbol[,2] == mouse_symbol[,3],]
human_symbol <- human_symbol[human_symbol[,2] == human_symbol[,3],]

# remove homologous genes that dont have a gene symbol
hom <- hom[hom[,1] %in% mouse_symbol[,1],]
hom <- hom[hom[,2] %in% human_symbol[,1],]

# make lists of homologous genes with IDs and Symbols
homology <- merge(hom[,1:2],mouse_symbol[,1:2])
homology <- merge(homology, human_symbol[,c(1,2)])


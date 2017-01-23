######################################################################
# Parse links to databases from your list of gene symbols
######################################################################
# source ('/Users/abelvertesy/Github_repos/TheCorvinas/R/DatabaseLinke.R')

######################################################################
## If you encounter a bug, or something doesn't work, please let me know by raising an issue on Github/vertesy/TheCorvinas
## https://github.com/vertesy/TheCorvinas/issues/new?milestone=DatabaseLinke.R
######################################################################

# vector_of_gene_symbols = c("Oct4", "Dazl")

# User Setup ----------------------------------------------------------------------
BashScriptLocation = "/Users/abelvertesy/bin/run.sh"
# ALT USAGE link_fromToClilpboard = toClipboard(link_uniprot_mice(fromClipboard.as_vec(), writeOut = F))

# Static part of Query links ------------------------
HGNC_symbol_search = "http://www.genenames.org/cgi-bin/gene_search?search="
wikipedia = "http://en.wikipedia.org/w/index.php?search="

ensembl_multispecies = c("http://www.ensembl.org/Multi/Search/Results?q=",";site=ensembl")
grc37 = c("http://grch37.ensembl.org/Human/Search/Results?q=", ";site=ensembl;facet_feature_type=Gene;facet_species=Human")
grc38 = c("http://www.ensembl.org/Human/Search/Results?q=", ";site=ensembl;facet_feature_type=Gene;facet_species=Human")

grc_mm38 = c("http://www.ensembl.org/Mouse/Search/Results?q=", ";site=ensembl;facet_feature_type=;facet_species=Mouse")
grc_Zebra = c("http://www.ensembl.org/Search/Results?q=", ";site=ensembl;facet_feature_type=;facet_species=Zebrafish")


uniprot_mouse = c('http://www.uniprot.org/uniprot/?query=organism%3A"Mus+musculus+[10090]"+',"&sort=score")
# I left out (Mouse) part because it messes up with markdown, it still works.
uniprot_human = c('http://www.uniprot.org/uniprot/?query=organism%3A"Homo+sapiens+%28Human%29+[9606]"+',"&sort=score")
uniprot_zebra = c('http://www.uniprot.org/uniprot/?query=organism%3A"Danio+rerio+(Zebrafish)+(Brachydanio+rerio)+[7955]"+',"&sort=score")


# HELP: http://string-db.org/help/faq/#how-do-i-link-to-string. Find Species ID by: http://www.ncbi.nlm.nih.gov/taxonomy
STRING = "http://string-db.org/newstring_cgi/show_network_section.pl?identifier="
STRING_mouse_suffix = "&species=10090"
STRING_human_suffix = "&species=9606"
STRING_elegans_suffix = "&species=6239"

wormbase_search_prefix = "https://www.wormbase.org/search/gene/"

PUBMED_search_prefix = "https://www.ncbi.nlm.nih.gov/pubmed/?term="

# HGNC links ------------------------------------------------------------------------------------------------
link_HGNC <- function (vector_of_gene_symbols, writeOut = T, Open=!writeOut) { # Parse HGNC links to your list of gene symbols
	links = paste0(HGNC_symbol_search, vector_of_gene_symbols)
	if (writeOut) {
		bash_commands = paste0("open ", links)
		write.simple.append("", ManualName = BashScriptLocation)
		write.simple.append(bash_commands, ManualName = BashScriptLocation)
	} else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}

# ENSEMBL Links --------------------------------------------------------------------------------------------------------------------------------

link_ensembl_zebra <- function (vector_of_gene_symbols, writeOut = T, Open=!writeOut) { # Parse the latest ensembl (GRC38) links to your list of gene symbols
  links = paste0(grc_Zebra[1], vector_of_gene_symbols, grc_Zebra[2])
  if (writeOut) {
    bash_commands = paste0("open ", links)
    write.simple.append("", ManualName = BashScriptLocation)
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}

link_ensembl_mice <- function (vector_of_gene_symbols, writeOut = T, Open=!writeOut) { # Parse the latest ensembl (GRC38) links to your list of gene symbols
  links = paste0(grc_mm38[1], vector_of_gene_symbols, grc_mm38[2])
  if (writeOut) {
    bash_commands = paste0("open ", links)
    write.simple.append("", ManualName = BashScriptLocation)
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}

link_ensembl_mice <- function (vector_of_gene_symbols, writeOut = T, Open=!writeOut) { # Parse the latest ensembl (GRC38) links to your list of gene symbols
  links = paste0(grc_mm38[1], vector_of_gene_symbols, grc_mm38[2])
  if (writeOut) {
    bash_commands = paste0("open ", links)
    write.simple.append("", ManualName = BashScriptLocation)
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}

link_ensembl <- function (vector_of_gene_symbols, writeOut = T, Open=!writeOut) { # Parse the latest ensembl (GRC38) links to your list of gene symbols
  links = paste0(ensembl_multispecies[1], vector_of_gene_symbols, ensembl_multispecies[2])
  if (writeOut) {
    bash_commands = paste0("open ", links)
    write.simple.append("", ManualName =BashScriptLocation )
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}

link_ensembl.grc37 <- function (vector_of_gene_symbols, writeOut = T, Open=!writeOut) { # Parse ensembl GRC37 links to your list of gene symbols
	links = paste0(grc37[1], vector_of_gene_symbols, grc37[2])
	if (writeOut) {
		bash_commands = paste0("open ", links)
		write.simple.append("", ManualName = BashScriptLocation)
		write.simple.append(bash_commands, ManualName = BashScriptLocation)
	} else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}

## UNIPROT Links --------------------------------------------------------------------------------------------------------------------------------

link_uniprot_mice <- function (vector_of_gene_symbols, writeOut = F, Open=!writeOut) { # Parse the latest UNIPROT links to your list of gene symbols
  links = paste0(uniprot_mouse[1], vector_of_gene_symbols, uniprot_mouse[2] )
  if (writeOut) {
    bash_commands = paste0("open ", links)
    write.simple.append("", ManualName = BashScriptLocation)
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}

link_uniprot_human <- function (vector_of_gene_symbols, writeOut = F, Open=!writeOut) { # Parse the latest UNIPROT links to your list of gene symbols
  links = paste0( uniprot_human[1], vector_of_gene_symbols, uniprot_human[2] )
  if (writeOut) {
    bash_commands = paste0("open ", links)
    write.simple.append("", ManualName = BashScriptLocation)
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}

link_uniprot_zebrafish <- function (vector_of_gene_symbols, writeOut = F, Open=!writeOut) { # Parse the latest UNIPROT links to your list of gene symbols
  links = paste0( uniprot_zebra[1], vector_of_gene_symbols, uniprot_zebra[2] )
  if (writeOut) {
    bash_commands = paste0("open ", links)
    write.simple.append("", ManualName = BashScriptLocation)
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}


# SRING links ------------------------------------------------------------------------
link_String <- function (vector_of_gene_symbols, organism="mouse", writeOut = T, Open=!writeOut) { # Parse STRING protein interaction database links to your list of gene symbols. "organism" can be mouse, human or NA
  suffix = if (is.na(organism)) { "" } 
  else if (organism== "elegans") { STRING_elegans_suffix } 
  else if (organism== "mouse") { STRING_mouse_suffix } 
  else if (organism== "human") { STRING_human_suffix }
    links = paste0( STRING, vector_of_gene_symbols, suffix )
  if (writeOut) {
    bash_commands = paste0("open '", links, "'")
    write.simple.append("", ManualName = BashScriptLocation)
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}


# PUBMED links ------------------------------------------------------------------------

link_pubmed <- function (vector_of_gene_symbols, additional_terms = "", writeOut = T, Open=!writeOut) { # Parse PUBMED database links to your list of gene symbols. "additional_terms" can be any vector of strings that will be searched for together with each gene.
  links = paste0( PUBMED_search_prefix, vector_of_gene_symbols, additional_terms )
  if (writeOut) {
    bash_commands = paste0("open '", links, "'")
    write.simple.append("", ManualName = BashScriptLocation)
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}

# wormbase links ------------------------------------------------------------------------

link_wormbase <- function (vector_of_gene_symbols, writeOut = T, Open=!writeOut) { # Parse wormbase database links to your list of gene symbols. "additional_terms" can be any vector of strings that will be searched for together with each gene.
  links = paste0( wormbase_search_prefix, vector_of_gene_symbols)
  if (writeOut) {
    bash_commands = paste0("open '", links, "'")
    write.simple.append("", ManualName = BashScriptLocation)
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}

# Wikipedia links ------------------------------------------------------------------------

link_wikipedia <- function (vector_of_gene_symbols, writeOut = T, Open=!writeOut) { # Parse wormbase database links to your list of gene symbols. "additional_terms" can be any vector of strings that will be searched for together with each gene.
  links = paste0( wikipedia, vector_of_gene_symbols)
  if (writeOut) {
    bash_commands = paste0("open '", links, "'")
    write.simple.append("", ManualName = BashScriptLocation)
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}

# CGC links (worms mutant database) ------------------------------------------------------------------------

worm_CGC_prefix = c("http://www.cgc.cbs.umn.edu/search.php?st=","&field=all&exst=&exfield=all")

link_CGC <- function (vector_of_gene_symbols, writeOut = T, Open=!writeOut) { # Parse CGC links (worms mutant database).
  links = paste0( worm_CGC_prefix[1], vector_of_gene_symbols, worm_CGC_prefix[2] )
  if (writeOut) {
    bash_commands = paste0("open '", links, "'")
    write.simple.append("", ManualName = BashScriptLocation)
    write.simple.append(bash_commands, ManualName = BashScriptLocation)
  } else if (Open) { for (l in links) browseURL(l) }	else { return(links) }
}



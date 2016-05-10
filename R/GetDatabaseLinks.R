######################################################################
# Parse links to databases from your list of human gene symbols
######################################################################
# source ('/Users/abelvertesy/Github_repos/TheCorvinas/R/GetDatabaseLinks.R')



# Functions ------------------------
HGNC_symbol_search = "http://www.genenames.org/cgi-bin/gene_search?search="
grc37 = c("http://grch37.ensembl.org/Human/Search/Results?q=", ";site=ensembl;facet_feature_type=Gene;facet_species=Human")
grc38 = c("http://www.ensembl.org/Human/Search/Results?q=", ";site=ensembl;facet_feature_type=Gene;facet_species=Human")

grc_mm38 = c("http://www.ensembl.org/Mouse/Search/Results?q=", ";site=ensembl;facet_feature_type=;facet_species=Mouse")



# Functions ------------------------
link_HGNC <- function (vector_of_gene_symbols, writeOut = T) { # Parse HGNC links to your list of gene symbols
	links = paste0(HGNC_symbol_search, vector_of_gene_symbols)
	if (writeOut) {
		bash_commands = paste0("open ", links)
		write.simple.append("", ManualName = "/Users/abelvertesy/bin/run.sh")
		write.simple.append(bash_commands, ManualName = "/Users/abelvertesy/bin/run.sh")
	} else { return(links) }
}


link_ensembl_mice <- function (vector_of_gene_symbols, writeOut = T) { # Parse the latest ensembl (GRC38) links to your list of gene symbols
  links = paste0(grc_mm38[1], vector_of_gene_symbols, grc_mm38[2])
  if (writeOut) {
    bash_commands = paste0("open ", links)
    write.simple.append("", ManualName = "/Users/abelvertesy/bin/run.sh")
    write.simple.append(bash_commands, ManualName = "/Users/abelvertesy/bin/run.sh")
  } else { return(links) }
}

link_ensembl <- function (vector_of_gene_symbols, writeOut = T) { # Parse the latest ensembl (GRC38) links to your list of gene symbols
  links = paste0(grc38[1], vector_of_gene_symbols, grc38[2])
  if (writeOut) {
    bash_commands = paste0("open ", links)
    write.simple.append("", ManualName = "/Users/abelvertesy/bin/run.sh")
    write.simple.append(bash_commands, ManualName = "/Users/abelvertesy/bin/run.sh")
  } else { return(links) }
}


link_ensembl.grc37 <- function (vector_of_gene_symbols, writeOut = T) { # Parse ensembl GRC37 links to your list of gene symbols
	links = paste0(grc37[1], vector_of_gene_symbols, grc37[2])
	if (writeOut) {
		bash_commands = paste0("open ", links)
		write.simple.append("", ManualName = "/Users/abelvertesy/bin/run.sh")
		write.simple.append(bash_commands, ManualName = "/Users/abelvertesy/bin/run.sh")
	} else { return(links) }
}




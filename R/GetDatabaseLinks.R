######################################################################
# Parse links to databases from your list of human gene symbols
######################################################################
# source ('/Users/abelvertesy/Github_repos/TheCorvinas/R/GetDatabaseLinks.R')



# Static part of Query links ------------------------
HGNC_symbol_search = "http://www.genenames.org/cgi-bin/gene_search?search="
grc37 = c("http://grch37.ensembl.org/Human/Search/Results?q=", ";site=ensembl;facet_feature_type=Gene;facet_species=Human")
grc38 = c("http://www.ensembl.org/Human/Search/Results?q=", ";site=ensembl;facet_feature_type=Gene;facet_species=Human")

grc_mm38 = c("http://www.ensembl.org/Mouse/Search/Results?q=", ";site=ensembl;facet_feature_type=;facet_species=Mouse")


grc_Zebra = c("http://www.ensembl.org/Mouse/Search/Results?q=", ";site=ensembl;facet_feature_type=;facet_species=Zebrafish")


uniprot_mouse = c('http://www.uniprot.org/uniprot/?query=organism%3A"Mus+musculus+[10090]"+',"&sort=score")
# I left out (Mouse) part because it messes up with markdown, it still works.
uniprot_human = c('http://www.uniprot.org/uniprot/?query=organism%3A"Homo+sapiens+%28Human%29+[9606]"+',"&sort=score")


STRING = "http://string-db.org/newstring_cgi/show_network_section.pl?identifier="

# HELP: http://string-db.org/help/faq/#how-do-i-link-to-string
STRING_mouse_suffix = "&species=10090"
STRING_human_suffix = "&species=9606"

# Functions ------------------------
link_HGNC <- function (vector_of_gene_symbols, writeOut = T) { # Parse HGNC links to your list of gene symbols
	links = paste0(HGNC_symbol_search, vector_of_gene_symbols)
	if (writeOut) {
		bash_commands = paste0("open ", links)
		write.simple.append("", ManualName = "/Users/abelvertesy/bin/run.sh")
		write.simple.append(bash_commands, ManualName = "/Users/abelvertesy/bin/run.sh")
	} else { return(links) }
}

# ENSEMBL Links --------------------------------------------------------------------------------------------------------------------------------

link_ensembl_zebra <- function (vector_of_gene_symbols, writeOut = T) { # Parse the latest ensembl (GRC38) links to your list of gene symbols
  links = paste0(grc_mm38[1], vector_of_gene_symbols, grc_Zebra[2])
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

## UNIPROT Links --------------------------------------------------------------------------------------------------------------------------------

link_uniprot_mice <- function (vector_of_gene_symbols, writeOut = T) { # Parse the latest UNIPROT links to your list of gene symbols
  links = paste0(uniprot_mouse[1], vector_of_gene_symbols, uniprot_mouse[2] )
  if (writeOut) {
    bash_commands = paste0("open ", links)
    write.simple.append("", ManualName = "/Users/abelvertesy/bin/run.sh")
    write.simple.append(bash_commands, ManualName = "/Users/abelvertesy/bin/run.sh")
  } else { return(links) }
}

link_uniprot_human <- function (vector_of_gene_symbols, writeOut = T) { # Parse the latest UNIPROT links to your list of gene symbols
  links = paste0( uniprot_human[1], vector_of_gene_symbols, uniprot_human[2] )
  if (writeOut) {
    bash_commands = paste0("open ", links)
    write.simple.append("", ManualName = "/Users/abelvertesy/bin/run.sh")
    write.simple.append(bash_commands, ManualName = "/Users/abelvertesy/bin/run.sh")
  } else { return(links) }
}


link_String <- function (vector_of_gene_symbols, organism="mouse", writeOut = T) { # Parse STRING protein interaction database links to your list of gene symbols. "organism" can be mouse, human or NA
  suffix = if (is.na(organism)) { "" } else if (organism== "mouse") { STRING_mouse_suffix } else if (organism== "human") { STRING_human_suffix }
    links = paste0( STRING, vector_of_gene_symbols, suffix )
  if (writeOut) {
    bash_commands = paste0("open '", links, "'")
    write.simple.append("", ManualName = "/Users/abelvertesy/bin/run.sh")
    write.simple.append(bash_commands, ManualName = "/Users/abelvertesy/bin/run.sh")
  } else { return(links) }
}

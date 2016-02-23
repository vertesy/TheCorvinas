
# RNA-seq specific R functions -----------------------------------------------------------------------------------------------------
# source("/Users/abelvertesy/TheCorvinas/R/RNA_seq_specific_functions.r")


chop_chr_from_gene_name <- function  (name, splitcharacter = "__") { # Chop the chromosome ending!
	strsplit(name, splitcharacter)[[1]][1] }

write.simple.append.vcf  <- function(input_df, extension='vcf', ManualName ="", ... ){ # use: arg1 data, arg's... strings to be concatenated to yield the path and FILE NAME
	fname = kollapse (...) ; if (nchar (fname) < 2 ) { fname = substitute(input_df) }
	if (nchar(ManualName)) {FnP = kollapse(ManualName)} else  { FnP = FnP_parser (fname, extension) }
	write.table (input_df, file = FnP, sep = "\t", row.names = F, col.names = F, quote=FALSE, append = T )
} # fun

fix_missing_entries <- function (complete_vec, partial_vec) { # if there are some categories missing, by creating a table from a vector, you can add the missing categories
	nr_cat = length (complete_vec)
	fixed_vec=rep(NA, nr_cat); names (fixed_vec) = names (complete_vec)
	for (n in names(complete_vec)) {
			fixed_vec[n] = partial_vec[n]
			partial_vec[n]
			if ( is.na(partial_vec[n]) ) {fixed_vec[n] = 0 }
	} # for
	return (fixed_vec)
}

wbarplot_cellID <-  function(variable, col ="gold1", ...) { # in ... you can pass on ANY plotting parameter exc SUB, MAIN!!!!
	plotname = kollapse(substitute(variable),"-", trunk[i], print =F)
	FnP = kollapse (OutDir,"/", plotname, "-", trunk[i], ".barplot.pdf", print=F)
	cexNsize = 0.7/abs (log10 (length(variable)) ); cexNsize = min (cexNsize, 1)
	barplot (variable, ..., main= plotname, col=col, las=2, cex.names = cexNsize,
					 sub = paste ("mean:", iround(mean(variable, na.rm=T)),  "CV:", percentage_formatter(cv(variable)) ) )
	dev.copy2pdf (file=FnP, width=w, height=h )
}
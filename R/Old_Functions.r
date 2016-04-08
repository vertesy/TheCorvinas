
## _Old_Functions.r
# source("/Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/Data_analysis/X_inact/Scripts_Xreact/zz_Old_versions/_Old_Functions.r")


# old <- `[`
'If the rownames are retained, the subselected column and row remains the original type of the object, eg df. hist () does not take df-s for instance'
# `[` <- function(...) { old(..., drop=FALSE) }
'Convenience in the original [ function, objects are simplified df[,col] -> vector.'
# `[` <- old



normalize_per_column <- function (m, factor) { # Normalize each column by its by a corresponding value. It is solely t(t(matrix)/factor) with a couple of error checks.
	any_print("Range of normalization factors:" , range(factor))
	stopifnot(dim(m)[2] == l(factor))
	stopifnot(any(is.finite(factor)))
	if( any(factor == 0) ) { print("Your normalization factor contains 0-s: NA-s will be produced.")}
	t(t(m)/factor)
}


secondDerivative_discrete <- function (x) { # Approximate the "discrete 2nd derivative" or "Finite difference" by central differences. Source: https://stackoverflow.com/questions/4471993/compute-the-elbow-for-a-curve-automatically-and-mathematically?lq=1
	secondDerivative = rep(NA,l(x))
	for(i in 2:l(x-1)){
		secondDerivative[i] = x[i+1] + x[i-1] - 2 * x[i]
	}
	return(secondDerivative)
}




read.simple2 <-  function(header=FALSE, ...) {
	pfn = kollapse (...) # merge path and filename
	read_in = unlist(read.table( pfn ,stringsAsFactors=FALSE, header=header ) )
	return(read_in)
}

help.cast <-  function() any_print( 'acast (long_data, rows ~ columns, colum-of-observations, filter_column/subset)')


write.simple.wRowNames = write.simple.tsv

# kollapse <- function(...,print =T) {
# 	if (print==T) {print (paste(c(...), sep="",collapse="")) }
# 	paste(c(...), sep="",collapse="")
# }

# setup_logging file, path and modification date
setup_logging <- function  (append=TRUE, ...) {
	Log_PnF <- kollapse (...)
	write (kollapse("                   Created: ",date() ), Log_PnF , append=append)
	assign ("Log_PnF",Log_PnF, envir = .GlobalEnv)
}

# # setup_logging file, path and modification date
# setup_logging2 <- function  (fname, append=T) {
# 	if ( exists('OutDir') ) { path = OutDir } else { path = getwd() ; any_print ("OutDir not defined !!!") }
# 	Log_PnF <- kollapse (path,'/',fname,'.log')
# 	write (kollapse("                   Modified: ",date() ), Log_PnF , append=append)
# 	assign ("Log_PnF",Log_PnF, envir = .GlobalEnv)
# }




# plot to screen and save as .pdf
wplot.nonCol <-  function(variable, plotname = substitute(variable), ..., w=7, h=7 ) {
	FnP = FnP_parser (plotname, 'plot.pdf')
	plot (variable, ..., main=plotname)
	dev.copy2pdf (file=FnP, width=w, height=h )
}


# plot with fitted line
wplot2 <-  function(variable, extension =".pdf", plotname = substitute(variable), col=rgb(0,0,0,75,maxColorValue=100), w=7, h=7, pch='.') {
	FnP = FnP_parser (plotname, 'plot.pdf')
	regr=lm(variable[,2]~variable[,1])
	subt=as.character(round(unlist(regr[1]), digits=2) )
	plot (variable, col=col, main=plotname, sub=paste (subt, collapse="   "), pch = pch, cex=1)
	abline(1,1, col=3, lty=1)
	abline(regr, col =2)
	suppressWarnings(rug(jitter(variable[,1]), side=1, col=rgb(100,100,100,50,maxColorValue=255)))
	suppressWarnings(rug(jitter(variable[,2]), side=2, col=rgb(100,100,100,50,maxColorValue=255)))
	dev.copy2pdf (file=FnP, width=w, height=h )
}

# plot with fitted line AND fixed axis
wplot3_fixAx <-  function(variable, extension =".pdf", plotname = substitute(variable), col=rgb(0,0,0,75,maxColorValue=100), w=7, h=7, fixed_axes=numeric(0), pch='.') {
	FnP = FnP_parser (plotname, 'plot.pdf')
	regr=lm(variable[,2]~variable[,1])
	subt=as.character(round(unlist(regr[1]), digits=2) )
	if (length (fixed_axes))  {
		plot (c(0,fixed_axes[1]), c(0,fixed_axes[2]), type="n", main=plotname,
			  sub=paste (subt, collapse="   "), xlab="REF depth", ylab="ALT depth")
		points (variable, col=col, pch = pch, cex=1)
	} else {
		plot (variable, col=col, main=plotname, sub=paste (subt, collapse="   "), pch = pch, cex=1)
	}
	abline(1,1, col=3, lty=1)
	abline(regr, col =2)
	suppressWarnings(rug(jitter(variable[,1]), side=1, col=rgb(100,100,100,50,maxColorValue=255)))
	suppressWarnings(rug(jitter(variable[,2]), side=2, col=rgb(100,100,100,50,maxColorValue=255)))
	dev.copy2pdf (file=FnP, width=w, height=h )
}


# plot with fitted line AND fixed axis
# Mat_DP[Filters_Nouvelle == 'ChildNotHet']
wplot4_fixAx <-  function(variable, extension =".pdf", plotname = substitute(variable), col=rgb(0,0,0,75,maxColorValue=100), w=7, h=7, fixed_axes=numeric(0), pch='.') {
	FnP = FnP_parser (plotname, 'plot.pdf')
	regr=lm(variable[,2]~variable[,1])
	subt=as.character(round(unlist(regr[1]), digits=2) )
	if (length (fixed_axes))  {
		plot (c(0,fixed_axes[1]), c(0,fixed_axes[2]), type="n", main=plotname,
			  sub=paste (subt, collapse="   "), xlab="PAT depth", ylab="MAT depth")
		points (variable, col=col, pch = pch, cex=1)
	} else {
		plot (variable, col=col, main=plotname, sub=paste (subt, collapse="   "), pch = pch, cex=1)
	}
	abline(1,1, col=3, lty=1)
	abline(regr, col =2)
	suppressWarnings(rug(jitter(variable[,1]), side=1, col=rgb(100,100,100,50,maxColorValue=255)))
	suppressWarnings(rug(jitter(variable[,2]), side=2, col=rgb(100,100,100,50,maxColorValue=255)))
	dev.copy2pdf (file=FnP, width=w, height=h )
}


wbarplot.nonCol <-  function(variable, plotname = substitute(variable), ..., w=7, h=7) {
	# in ... you can pass on ANY plotting parameter!!!!
	FnP = FnP_parser (plotname, 'barplot.pdf')
	barplot (variable, ..., main=plotname, las=2)
	dev.copy2pdf (file=FnP, width=w, height=h )
}

wimage <-  function(variable,  xlab ="", ylab ="", plotname = substitute(variable), ..., w=7, h=7, mdlink =F) {
	fname = kollapse (plotname, '.heatmap')
	dimMat= dim(variable)
	image (x =1:dimMat[1], y =1:dimMat[2], z=variable,
		   # col = HeatMapCol_BWR(10),
		   col = HeatMapCol_RedBlackBlue(10),
		   main=plotname, sub = kollapse ("Range of values: ", paste(signif(range(variable, na.rm=T),3), collapse = ", ") )
		   , xlab=xlab, ylab=ylab , las=2
		   , xaxp = c( range(1,dimMat[1] ), dimMat[1]-1)
		   , yaxp = c( range(1,dimMat[2] ), dimMat[2]-1)
	)
	dev.copy2pdf (file=FnP_parser (fname, 'pdf'), width=w, height=h )
	if (mdlink) { 	MarkDown_Img_Logger_PDF_and_PNG (fname_wo_ext = fname) }# put a markdown image link if the log file exists
}

simplify_categories <-   function(category_vec, replaceit , to ) { # Replace every instance of a value with another value provided
	# could be extended to set of needles & set of replacements (with lenght check, 2>1 is also OK)
	matches  = which(category_vec %in% replaceit); any_print(l(matches), "instances of", replaceit, "are replaced by", to)
	category_vec[matches] =  to
	return(category_vec)
}


as.numeric.convertString <- function (vec) { # Converts any vector into a numeric vector, and puts the original character values into the names of the new vector, unless it already has names. Useful for coloring a plot by categories, name-tags, etc.
	numerified_vec = as.numeric(as.factor(vec))
	if (!is.null(names(vec))) {names (numerified_vec) = names (vec)}
	return(numerified_vec)
}

# WO as.factor!
as.numeric.wNames <- function (vec) { # Converts any vector into a numeric vector, and puts the original character values into the names of the new vector, unless it already has names. Useful for coloring a plot by categories, name-tags, etc.
	numerified_vec = as.numeric((vec))
	if (!is.null(names(vec))) {names (numerified_vec) = names (vec)}
	return(numerified_vec)
}



# df_col_to_vector_w_names <- function (OneColDF, WhichDimension ="col") {
# 	if (WhichDimension == "col") {	n = rownames(OneColDF)} # get the rownames for that column of the DF
# 	if (WhichDimension == "row") {	n = colnames(OneColDF)} # get colnames...
# 	OneColDF = as.vector(unlist(OneColDF));	names(OneColDF) = n 				# add names
# 	return (OneColDF)
# 	assign (substitute(OneColDF), OneColDF, envir = .GlobalEnv)
# }
df_col_to_vector_w_names = as.named.vector

stdError = sem


shannon.entropy <- function(p) {
	if (min(p) < 0 || sum(p) <= 0)
		return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
}


whist.nonCol <-  function(variable, plotname = substitute(variable), ..., w=7, h=7) {
	FnP = FnP_parser (plotname, 'hist.pdf')
	hist (variable, ..., main=plotname)
	dev.copy2pdf (file=FnP, width=w, height=h )
}

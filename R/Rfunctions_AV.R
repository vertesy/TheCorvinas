######################################################################
# Custom made functions for R
# - in categories
# - source it when processing the data

######################################################################
# source ('/Users/abelvertesy/TheCorvinas/R/Rfunctions_AV.R')
# cp .Rfunctions_Base.09.R /Users/abelvertesy/Dropbox_at_open/Dropbox/X_reactivation/scripts/Abel/Rfunctions_Base.10.R

# CHAPTERS:
### quick help / interpretatio
### File handling [read & write]
### Math $ stats
### Printing and Strings
### Generic
### Plotting and Graphics
### RNA-seq specific

# quick help / interpretatio  -------------------------------------------------------------------------------------------------

help.cast <-  function() any_print( 'acast (long_data, rows ~ columns, colum-of-observations, filter_column/subset)')
l=length

sortbyitsnames  <-  function(vec) {vec[order(names(vec) )]}

stopif  <-  function(condition, message ="") { if(condition) {any_print (message); stop()} }

### MarkDownLogg.R Library-------------------------------------------------------------------------------------------------
source("/Users/abelvertesy/TheCorvinas/R/MarkDownLogger_functions.R")

### File handling [read & write] -------------------------------------------------------------------------------------------------

# old <- `[`
'If the rownames are retained, the subselected column and row remains the original type of the object, eg df. hist () does not take df-s for instance'
# `[` <- function(...) { old(..., drop=FALSE) }
'Convenience in the original [ function, objects are simplified df[,col] -> vector.'
# `[` <- old


attach_w_rownames <- function (df_w_dimnames) {
	if(!is.null(rownames(df_w_dimnames)) & !is.null(colnames(df_w_dimnames))) {
		namez= rownames(df_w_dimnames)
		any_print("Now directly available in the workspace:      ", colnames(df_w_dimnames))
		attach (df_w_dimnames)
		for (n in colnames(df_w_dimnames)) {
			x=get(n); names(x) = namez
			assign (n,x,envir =.GlobalEnv) } # for
	} else { print ("ERROR: the DF does not have some of the dimnames!")}
}

Color_Check <- function (...) {
	Numbers  = c(...)
	barplot (rep(10,length(Numbers)), col =Numbers, xlab = paste (Numbers, collapse="") )
}


Clipboard_Copy <- function (x, sep="\t", header=FALSE, row.names=FALSE, col.names =F) {
	write.table(x, pipe("pbcopy"), sep=sep, row.names=row.names, col.names =col.names)
}

Clipboard_Paste <- function ( sep="\t", header=F) {
	return (read.table(pipe("pbpaste"), sep=sep, header=header, stringsAsFactors =F))
}

Clipboard_Paste_vec <- function ( sep="\t", header=F) {
	return (as.vector(unlist(read.table(pipe("pbpaste"), sep=sep, header=header, stringsAsFactors =F))))
}

Clipboard_Paste_num_vec <- function ( sep="\t", header=F) {
	return (as.numeric(unlist(read.table(pipe("pbpaste"), sep=sep, header=header, stringsAsFactors =F))))
}


FnP_parser <-  function(fname, ext_wo_dot) {
	# parses Filename & Path for output file
	if ( exists('OutDir') ) { path = OutDir } else { path = getwd() ; any_print ("OutDir not defined !!!") }
	if (hasArg(ext_wo_dot) ) { FnP = kollapse (path,"/", fname, ".", ext_wo_dot)
	} else { 					FnP = kollapse (path,"/", fname) }
}

read.simple.vec <-  function(...) {
	pfn = kollapse (...) # merge path and filename
	read_in = as.vector(unlist(read.table( pfn ,stringsAsFactors=F )) )
	any_print(length (read_in), "elements")
	return(read_in);
}

read.simple <-  function(...) {
	pfn = kollapse (...) # merge path and filename
	read_in = unlist(read.table( pfn ,stringsAsFactors=F ) )
	return(read_in)
}

read.simple2 <-  function(header=FALSE, ...) {
	pfn = kollapse (...) # merge path and filename
	read_in = unlist(read.table( pfn ,stringsAsFactors=FALSE, header=header ) )
	return(read_in)
}

read.simple_char_list <-  function(...) {
	pfn = kollapse (...) # merge path and filename
	read_in = unlist(read.table( pfn ,stringsAsFactors=F ) )
	any_print ("New variable head: ",what(read_in))
	return(read_in)
}

read.simple.table <-  function(...,rownames=NULL, colnames=T) {
	# default: header defines colnames, no rownames. For rownames give the col nr. with rownames, eg. 1
	pfn = kollapse (...) # merge path and filename
	read_in = read.table( pfn ,stringsAsFactors=FALSE, sep="\t", row.names= rownames, header=colnames )
	any_print ("New variable dim: ",dim(read_in))
	return(read_in)
}

read.simple.tsv <-  function(...) {
# for excel style data: rownames in col1, headers SHIFTED
	pfn = kollapse (...)
	read_in = read.delim( pfn ,stringsAsFactors=FALSE, sep="\t", row.names=1, header=T )
	# HEADER SHOULD start like \t colname1 \t ...
	any_print ("New variable dim: ",dim(read_in))
	return(read_in)
}

read.simple.tsv.named.vector <-  function(...) {
	# for excel style named vectors, names in col1, headers SHIFTED
	pfn = kollapse (...)
	read_in = read.delim( pfn ,stringsAsFactors=FALSE, sep="\t", row.names=1, header=T )
	# HEADER SHOULD start like \t colname1 \t ...
	rn = row.names(read_in)
	read_in =  as.vector(unlist(read_in));
	names(read_in) = rn
	any_print ("New vectors length is: ",length(read_in))
	return(read_in)
}

write.simple  <- function(input_df, extension='tsv', ManualName ="", ...  ) {
	# use: arg1 data, arg's... strings to be concatenated to yield the path and FILE NAME
	fname = kollapse (...) ; if (nchar (fname) < 2 ) { fname = substitute(input_df) }
	if (nchar(ManualName)) {FnP = kollapse(ManualName)} else  { FnP = FnP_parser (fname, extension) }
	write.table (input_df, file = FnP, sep = "\t", row.names = F, col.names = T, quote=FALSE  )
	any_print ("Length: ", length(input_df))
} # fun

write.simple.vec  <- function(input_vec, extension='vec', ManualName ="", ... ) {
	# use: arg1 data, arg's... strings to be concatenated to yield the path and FILE NAME
	fname = kollapse (...) ; if (nchar (fname) < 2 ) { fname = substitute(input_vec) }
	if (nchar(ManualName)) {FnP = kollapse(ManualName)} else  { FnP = FnP_parser (fname, extension) }
	write.table (input_vec, file = FnP, sep = "\t", row.names = F, col.names = F, quote=FALSE  )
	any_print ("Length: ", length(input_vec))
} # fun

write.simple.tsv  <- function(input_df, extension='tsv', ManualName ="", ... ) {
	# ROW & COL names, + blank 1st colname, fname is optional [var names is used instead]
	fname = kollapse (..., print = F); if (nchar (fname) < 2 ) { fname = substitute(input_df) }
	if (nchar(ManualName)) {FnP = kollapse(ManualName)} else  { FnP = FnP_parser (fname, extension) }
	write.table (input_df, file = FnP, sep = "\t", row.names = T, col.names = NA, quote=FALSE  )
	any_print ("Dim: ", dim(input_df))
} # fun

# If col.names = NA and row.names = TRUE a blank column name is added, which is the convention used for CSV files to be read by spreadsheets.

# write.simple.wRowNames = write.simple.tsv

write.simple.append  <- function(input_df, extension='tsv', ManualName ="", ... ) {
	# NO ROW names or else,  fname is optional [var names is used instead]
	fname = kollapse (...) ; if (nchar (fname) < 2 ) { fname = substitute(input_df) }
	if (nchar(ManualName)) { FnP = kollapse(ManualName)} else  { FnP = FnP_parser (fname, extension) }
	write.table (input_df, file = FnP, sep = "\t", row.names = F,col.names = F, quote=FALSE, append=T  )
} # fun

# InCodeDataFormat_num InCodeDataFormat_txt
inline_vec_char <- function(char_vector) {	Clipboard_Copy(print(paste("c( '", paste (char_vector, collapse =  "', '"),  "')", collapse = "", sep=""), quote = F)); print(" Copied to Clipboard") }
inline_vec_num <- function(num_vector)  { Clipboard_Copy(print(paste("c( ", paste (num_vector, collapse =  ", "),  " )", collapse = "", sep=""), quote = F)); print(" Copied to Clipboard") }
inline_list_char <- function(char_list) {
	print ("list(", quote = F)
	for (l in 1: length(list)) {
		print(paste("c( '", paste (char_list[[l]], collapse =  "', '"),  "')", collapse = "", sep=""), quote = F)
	};	print (")", quote = F)
}

# Math $ stats -------------------------------------------------------------------------------------------------

which_names <- function (named_Vec) {return(names(which(named_Vec)))}

as.named.vector <- function (df_col, WhichDimNames = 1) {
	# use RowNames: WhichDimNames = 1 , 2: use ColNames
	# !!! might require drop=F in subsetting!!! eg: df_col[,3, drop=F]
	namez = dimnames(df_col)[[WhichDimNames]]
	if (is.list(df_col) & !is.data.frame(df_col)) {namez = names(df_col)}
	vecc = as.vector(unlist (df_col))
	names (vecc)= namez
	return (vecc)
}

as.numeric.convertString <- function (vec) {
	numerified_vec = as.numeric(as.factor(vec))
	if (!is.null(names(cell_IDs))) {names (numerified_vec) = names (vec)}
	return(numerified_vec)
}

as.numeric.wNames <- function (vec) {
	numerified_vec = as.numeric(vec)
	if (!is.null(names(vec))) {names (numerified_vec) = names (vec)}
	return(numerified_vec)
}

as.logical.wNames <- function (vec) {
	numerified_vec = as.logical(vec)
	if (!is.null(names(vec))) {names (numerified_vec) = names (vec)}
	return(numerified_vec)
}

as.character.wNames <- function (vec) {
	char_vec = as.character(vec)
	if (!is.null(names(vec))) {names (char_vec) = names (vec)}
	return(char_vec)
}

# df_col_to_vector_w_names <- function (OneColDF, WhichDimension ="col") {
# 	if (WhichDimension == "col") {	n = rownames(OneColDF)} # get the rownames for that column of the DF
# 	if (WhichDimension == "row") {	n = colnames(OneColDF)} # get colnames...
# 	OneColDF = as.vector(unlist(OneColDF));	names(OneColDF) = n 				# add names
# 	return (OneColDF)
# 	assign (substitute(OneColDF), OneColDF, envir = .GlobalEnv)
# }
df_col_to_vector_w_names = as.named.vector


idim  <- function (any_object) {
	if (is.null(dim(any_object))) {print(l(any_object))}
	else {	print(dim(any_object))	}
}


idimnames  <- function (any_object) {
	if (!is.null(dimnames(any_object))) 	{ print(dimnames(any_object)) }
	else if (!is.null(colnames(any_object))) { any_print("colnames:", colnames(any_object))	}
	else if (!is.null(rownames(any_object))) { any_print("rownames:", rownames(any_object))	}
	else if (!is.null(names(any_object))) { any_print("names:", names(any_object))	}
}

top_indices <- function(x, n = 3, top = T){
	# for equal values, it maintains the original order
	head( order(x, decreasing =top), n )
}

stdError <- function(x) sd(x)/sqrt(length(x))
sem = stdError

cv <- function(x) ( sd(x, na.rm=T)/mean(x, na.rm=T) ) # change this!

fano <- function(x) ( var(x, na.rm=T)/mean(x, na.rm=T) )

modus <- function(x) {
	x= unlist(na.exclude(x))
	ux <- unique(x)
	tab <- tabulate(match(x, ux));
	ux[tab == max(tab)]
}

gm_mean <- function(x, na.rm=TRUE){ 	exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) }

# movingAve <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
movingAve <- function(x, oneSide) { y = NULL
	for (i in oneSide:l(x)) {
		y[i] = mean( x[ (i-oneSide):(i+oneSide) ] )
	}; 	return (y)
}

movingSEM <- function(x, oneSide) { y = NULL
	for (i in oneSide:l(x)) {
		y[i] = stdError( x[ (i-oneSide):(i+oneSide) ] )
	}; 	return (y)
}

imovingSEM <- function(x, oneSide = 5) { y = NULL
for (i in 1:l(x)) {
	oneSideDynamic = min(i-1,oneSide, l(x)-i); oneSideDynamic
	# any_print(i, " : ", oneSideDynamic)
	indexx = (i-oneSideDynamic):(i+oneSideDynamic);indexx
	y[i] = stdError( x[ indexx ] )
}; 	return (y)
}

remove_outliers <- function(x, na.rm = TRUE, ..., probs = c(.05, .95)) {
  qnt <- quantile(x, probs=probs, na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

colSort <- function(data, ...) sapply(data, sort, ...)

colMedians <- function(mat,na.rm=TRUE) return(apply(mat,2,median, na.rm=na.rm))
rowMedians <- function(mat,na.rm=TRUE) return(apply(mat,1,median, na.rm=na.rm))

colCV <- function(mat) return(apply(mat,2,cv ) )

rowMin <- function(x) {
	code = paste("x[,",1:(NCOL(x)),"]",sep="",collapse=",")	# Construct a call pmin(x[,1],x[,2],...x[,NCOL(x)])
	code = paste("pmin(",code,")"); 	return(eval(parse(text=code)))
}

rowMax <- function(x) {
	code = paste("x[,",1:(NCOL(x)),"]",sep="",collapse=",")
	code = paste("pmax(",code,")");	return(eval(parse(text=code)))
}

colMin <- function(x) {
	code = paste("x[",1:(NCOL(x)),",]",sep="",collapse=",")
	code = paste("pmin(",code,")");	return(eval(parse(text=code)))
}

colMax <- function(x) {
	code = paste("x[",1:(NCOL(x)),",]",sep="",collapse=",")
	code = paste("pmax(",code,")");	return(eval(parse(text=code)))
}

pc_TRUE <- function (logical_vector, percentify =T) {
	out = sum(logical_vector, na.rm=T) / length(logical_vector)
	if (percentify) {out = percentage_formatter (out) }
	return(out)
	}

rescale <- function (vec, from=0, upto=100) {
	# linear transformation to a given range of values
	vec = vec-min(vec, na.rm = T)
	vec = vec*((upto-from)/max(vec, na.rm = T))
	vec = vec+ from
	return (vec)
} # fun

filter_survival_length <- function (length_new, length_old, prepend ="") {
	pc = percentage_formatter(length_new/length_old)
	llprint (prepend, pc, " of ",length_old," entries make through the filter")
}

filter_HP <- function (vector, threshold, prepend ="", survival=F) {
	survivors = vector>threshold
	pc = percentage_formatter(sum(survivors)/length(survivors))
	if (file.exists(Log_PnF) ) {	llprint (prepend, pc, " or ", sum(survivors), " of ",length(vector)," entries in ", substitute (vector)," fall above a threshold value of: ", threshold)
		} else { any_print  (pc, " of ",length(vector)," entries in ", substitute (vector)," fall above a threshold value of: ", threshold, "NOT LOGGED") }
	if (survival) {return (sum(survivors)/length(survivors))} else if (!survival) { return (survivors) }
}

filter_LP <- function (vector, threshold, prepend ="", survival=F) {
	survivors = vector<threshold
	pc = percentage_formatter(sum(survivors)/length(survivors))
	if (file.exists(Log_PnF) ) {	llprint (prepend, pc, " or ", sum(survivors), " of ",length(vector)," entries in ", substitute (vector)," fall below a threshold value of: ", threshold )
		} else { any_print  (pc, " of ",length(vector)," entries in ", substitute (vector)," fall below a threshold value of: ", threshold, "NOT LOGGED") }
	if (survival) {return (sum(survivors)/length(survivors))} else if (!survival) { return (survivors) }
}

filter_MidPass <- function (vector, HP_threshold, LP_threshold, prepend ="", survival=FALSE, EdgePass = F) {
	# LP <= x < HP
	survivors = ( vector >= HP_threshold & vector < LP_threshold); keyword = "between"; relation = " <= x < "
	if (EdgePass) {survivors = ( vector <= HP_threshold | vector > LP_threshold); keyword = "outside"; relation = " >= x OR x > " }
	pc = percentage_formatter(sum(survivors)/length(survivors))
	Texxt = kollapse(pc, " or ", sum(survivors), " of ",length(vector)," entries in ", substitute (vector)," fall ", keyword, " the thresholds: ", HP_threshold, relation, LP_threshold, print = F)
	if (file.exists(Log_PnF) ) {	llprint (prepend, Texxt)
	} else { any_print  (Texxt, "NOT LOGGED") }
	if (survival) {return (sum(survivors)/length(survivors))} else if (!survival) { return (survivors) }
}


shannon.entropy <- function(p) {
	if (min(p) < 0 || sum(p) <= 0)
		return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
}


# Printing and Strings  -------------------------------------------------------------------------------------------------
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

eval_parse_kollapse <- dyn_var_caller <- function( ... ){
	substitute(eval(parse(text=kollapse( ... , print=F))))
}

substrRight <- function(x, n){
	substr(x, nchar(x)-n+1, nchar(x))
}

lm_equation_formatter <- function(lm) {
	eq = (lm$coefficients);
	kollapse ("Intercept:", eq[1], " Slope:", eq[2]);
}

na.omit.strip <- function (vec) {
	if (is.data.frame(vec)) {
		if ( min(dim(vec)) > 1 ) { any_print(dim(vec), "dimensional array is converted to a vector.") }
		vec = unlist(vec) }
	clean = na.omit(vec)
	attributes(clean)$na.action <- NULL
	return(clean)
}

inf.omit <- function (vec) {
	if (is.data.frame(vec)) {
		if ( min(dim(vec)) > 1 ) { any_print(dim(vec), "dimensional array is converted to a vector.") }
		vec = unlist(vec) }
	clean = vec[is.finite(vec)]
	# attributes(clean)$na.action <- NULL
	return(clean)
}



zero.omit <- function (vec) {
	v2= vec[vec!=0]
	any_print("range: ", range(v2))
	if ( !is.null(names(vec)) ) {names(v2) = names(vec)[vec!=0]}
	return(v2)
}

# Generic -------------------------------------------------------------------------------------------------

most_frequent_elements <- function(thingy, topN=10) {
	tail(sort(table(thingy, useNA = "ifany")), topN)
}

what <- function(x, printme=0) {
	# it can print the first "printme" elements
	any_print (is (x),"; nr. of elements:", length (x))
	if (is.numeric (x) ) 		{ any_print ("min&max:", range(x) ) } else {print ("Not numeric")}
	if ( length(dim(x) ) > 0 ) 	{ any_print ("Dim:", dim (x) )	}
	if ( printme>0) 			{ any_print ("Elements:", x[0:printme] )	}
	head (x)
}

lookup <- function(needle, haystack, exact =TRUE, report = FALSE) {
	ls_out = as.list( c(ln_needle = length(needle), ln_haystack = length(haystack), ln_hits = "",  hit_poz = "", hits = "") )
	Findings = numeric(0)
	ln_needle = length(needle)
	if (exact) {
		for (i in 1:ln_needle) {			Findings= c(Findings, which(haystack == needle[i]) )		} # for
	} else {
		for (i in 1:ln_needle) {			Findings = c(Findings,grep(needle[i], haystack,  ignore.case = T, perl = FALSE))		} # for
	} # exact or partial match
	ls_out$'hit_poz' = Findings
	ls_out$'ln_hits' = length(Findings)
	ls_out$'hits' = haystack[Findings]
	if (l(Findings)) {	ls_out$'nonhits' = haystack[-Findings]
	} else { 			ls_out$'nonhits' = haystack	}
	if (report){
		llprint (length(Findings), "/", ln_needle, '(', percentage_formatter(length(Findings)/ln_needle)
					 , ") of", substitute(needle), "were found among", length(haystack), substitute(haystack), "." )
		if (length(Findings)) { llprint( substitute(needle),"findings: ",paste ( haystack[Findings], sep= " " ) ) }
	} else { any_print(length(Findings), "Hits:", haystack[Findings]) } # if (report)
	return (ls_out)
} # lookup fun # lookup (needle, haystack); needle = c('a', 'd', 'c', 'v') ; haystack = c('as', 'ff', 'c', 'v', 'x')

# Plotting and Graphics -----------------------------------------------------------------------------------------------------

# HeatMapCol_BGR <- colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
# HeatMapCol_BWR <- colorRampPalette(c("blue", "white", "red"), bias=1)
# HeatMapCol_RedBlackGreen <- colorRampPalette(c("red", "black", "green"), bias=1)

val2col<-function(z, zlim, col = rev(heat.colors(12)), breaks){
	if(!missing(breaks)){
		if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
	}
	if(missing(breaks) & !missing(zlim)){
		breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
	}
	if(missing(breaks) & missing(zlim)){
		zlim <- range(z, na.rm=TRUE)
		zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
		zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
		breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
	}
	colorlevels <- col[((as.vector(z)-breaks[1])/(range(breaks)[2]-range(breaks)[1]))*(length(breaks)-1)+1] # assign colors to heights for each point
	colorlevels
}
# https://stackoverflow.com/questions/8717669/heat-map-colors-corresponding-to-data-in-r


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	stopifnot (length(x) == length(y) & length(y) ==length(lower) & length(lower) == length(upper))
	if (l(dim(y)) > 1 ) { # if a matrix
				arrows(as.vector(x),as.vector(y+upper), as.vector(x), as.vector(y-lower), angle=90, code=3, length=length, ...)
	} else if (l(dim(y)) == 1) { 	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...) }
}

barplot.label <- function(x, y, labels, bottom = F, relpos_top =.9, relpos_bottom =.1 ,...){
	stopifnot (length(x) == length(y))
	if (bottom) { y = rep (relpos_bottom * max(y, na.rm=T), length(x))} # if put labels at the foot
	if (l(dim(x)) > 1 ) { # if a matrix
		text(as.vector(x),as.vector(y * relpos_top), labels = as.vector(labels), ...)
	} else if (l(dim(x)) == 1) { 	text( (x), (y), labels = (labels), ...) }
}


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


whist.nonCol <-  function(variable, plotname = substitute(variable), ..., w=7, h=7) {
	FnP = FnP_parser (plotname, 'hist.pdf')
	hist (variable, ..., main=plotname)
	dev.copy2pdf (file=FnP, width=w, height=h )
}


whist_dfCol <-  function(df, colName, col ="gold", ..., w=7, h=7) {
	# name the file  by naming the variable! Cannot be used with dynamically called variables [e.g. call vectors within a loop]
	stopifnot(colName %in% colnames(df))
	variable = unlist(df[,colName])
	stopifnot(length (variable) >1 )
	plotname =  paste(substitute(df),'__', colName, sep="")
	FnP = FnP_parser (plotname, 'hist.pdf')
	if ( !is.numeric(variable)) { variable = table (variable) ;
							cexNsize = 0.7/abs (log10 (length(variable)) ); cexNsize = min (cexNsize, 1)
							barplot (variable, ..., main=plotname, col=col, las=2, cex.names = cexNsize,
							 sub = paste ("mean:", iround(mean(variable, na.rm=T)),  "CV:", percentage_formatter(cv(variable)) ) )
	} else {
		zz=hist (variable, ..., plot=F)
		hist (variable, ..., main=plotname, col=col, las=2,
					sub = paste ("mean:", iround(mean(zz$counts)),  "median:", iround(median(zz$counts)) ) )
	} # if is.numeric
	dev.copy2pdf (file=FnP, width=w, height=h )
}


wbarplot.nonCol <-  function(variable, plotname = substitute(variable), ..., w=7, h=7) {
	# in ... you can pass on ANY plotting parameter!!!!
	FnP = FnP_parser (plotname, 'barplot.pdf')
	barplot (variable, ..., main=plotname, las=2)
	dev.copy2pdf (file=FnP, width=w, height=h )
}


wbarplot_dfCol <-  function(df,colName, col ="gold1", w=7, h=7, ...) {
	stopifnot(colName %in% colnames(df))
	variable = unlist(df[,colName])
	stopifnot(length (variable) >1 )
	plotname =  paste(substitute(df),'__', colName, sep="")
	FnP = FnP_parser (plotname, 'barplot.pdf')
	cexNsize = 0.7/abs (log10 (length(variable)) ); cexNsize = min (cexNsize, 1)
	barplot (variable, ..., main=plotname, col=col, las=2, cex.names = cexNsize,
					 sub = paste ("mean:", iround(mean(variable, na.rm=T)),  "CV:", percentage_formatter(cv(variable)) ) )
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
HeatMapCol_RedBlackBlue <- colorRampPalette(c("red", "black", "green"), bias=1)


# Read and write plotting functions READ -------------------------------------
# rw-funcitons cannot call the w-functions, there would be too much things lost: rw writes where the file comes from, w writes in Outdir & current dir

rwplot <- function(FnP, ..., w=7, h=7) {
	print ('file without header, and file path between ""')
	variable=read.simple  (FnP);
	fname=gsub (".*/", "",FnP);  trunk = gsub ("\\..*", "",fname)
	plot (variable, ..., main=trunk )
	dev.copy2pdf (file=kollapse(FnP,".plot.pdf"), width=w, height=h )
}

rwscatterplot <- function(FnP, ..., w=7, h=7) {
	print ('file without header, and file path between ""')
	variable=read.simple.table  (FnP);
	fname=gsub (".*/", "",FnP);  trunk = gsub ("\\..*", "",fname)
	line = lm(variable[,2]~variable[,1])
	plot (variable, ..., main=trunk, sub = lm_equation_formatter (line) )
	abline(line)
	dev.copy2pdf (file=kollapse(FnP,".scatter.pdf"), width=w, height=h )
}

rwboxplot <- function(FnP, col ="gold1", ..., w=7, h=7) {
	print ('inputfile is a dataframe, with header')
	variable=read.simple.table (FnP);
	fname=gsub (".*/", "",FnP);  trunk = gsub ("\\..*", "",fname)
	boxplot (variable, ..., main=trunk, col =col, las=2)
	dev.copy2pdf (file=kollapse(FnP,".boxplot.pdf"), width=w, height=h )
}

rwhist <- function(FnP, col ="gold1", ..., w=7, h=7) {
	# print ('file without header, and file path between ""')
	variable=read.simple (FnP);
	fname=gsub (".*/", "",FnP);  plotname = gsub ("\\..*", "",fname)
	if ( length (variable) > 0 ) {
		if ( !is.numeric(variable)) { variable = table (variable) ; wbarplot (variable); print ("FILENAME IS: VARIABLE.BARPLOT.PDF") }
		else { 	hist (variable, ..., plotname=trunk, col =col, las=2)
				dev.copy2pdf (file=kollapse(FnP,".hist.pdf"), width=w, height=h )
		} # if is.numeric
	} # if non empty
}

rwbarplot <- function(FnP, col ="gold1", ..., w=7, h=7) {
	print ('file without header, and file path between ""')
	variable=read.simple (FnP);
	fname=gsub (".*/", "",FnP);  trunk = gsub ("\\..*", "",fname)
	barplot (variable, ..., main=trunk, col =col, las=2)
	dev.copy2pdf (file=kollapse(FnP,".barplot.pdf"), width=w, height=h )
}

# RNA-seq specific -----------------------------------------------------------------------------------------------------

# Chop the chromosome ending!
chop_chr_from_gene_name <- function  (name, splitcharacter = "__") { 	strsplit(name, splitcharacter)[[1]][1] }

write.simple.append.vcf  <- function(input_df, extension='vcf', ManualName ="", ... ){
	# use: arg1 data, arg's... strings to be concatenated to yield the path and FILE NAME
	fname = kollapse (...) ; if (nchar (fname) < 2 ) { fname = substitute(input_df) }
	if (nchar(ManualName)) {FnP = kollapse(ManualName)} else  { FnP = FnP_parser (fname, extension) }
	write.table (input_df, file = FnP, sep = "\t", row.names = F, col.names = F, quote=FALSE, append = T )
} # fun

fix_missing_entries <- function (complete_vec, partial_vec) {
	# if there are some categories missing, by creating a table from a vector, you can add the missing categories
	nr_cat = length (complete_vec)
	fixed_vec=rep(NA, nr_cat); names (fixed_vec) = names (complete_vec)
	for (n in names(complete_vec)) {
			fixed_vec[n] = partial_vec[n]
			partial_vec[n]
			if ( is.na(partial_vec[n]) ) {fixed_vec[n] = 0 }
	} # for
	return (fixed_vec)
}

wbarplot_cellID <-  function(variable, col ="gold1", ...) {
	# in ... you can pass on ANY plotting parameter exc SUB, MAIN!!!!
	plotname = kollapse(substitute(variable),"-", trunk[i], print =F)
	FnP = kollapse (OutDir,"/", plotname, "-", trunk[i], ".barplot.pdf", print=F)
	cexNsize = 0.7/abs (log10 (length(variable)) ); cexNsize = min (cexNsize, 1)
	barplot (variable, ..., main= plotname, col=col, las=2, cex.names = cexNsize,
					 sub = paste ("mean:", iround(mean(variable, na.rm=T)),  "CV:", percentage_formatter(cv(variable)) ) )
	dev.copy2pdf (file=FnP, width=w, height=h )
}

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
sort.mat <- function (df, colname_in_df = 1, decrease = F, na_last = T) {
	# ALTERNATIVE: dd[with(dd, order(-z, b)), ]
	# https://stackoverflow.com/questions/1296646/how-to-sort-a-dataframe-by-columns-in-r
	if (length(colname_in_df)>1) { print ("cannot handle multi column sort") }
	else {df[ order(df[,colname_in_df], decreasing = decrease, na.last = na_last), ]}
}

pc_in_total_of_match <- function (vec_or_table, category, NA_omit=T) {
	## percentage of a certain value within a vec_or_table
	if (is.table(vec_or_table)) { vec_or_table[category]/sum(vec_or_table, na.rm=NA_omit) }
	else { # if (is.vector(vec_or_table))
		if (NA_omit){
			if (sum(is.na(vec_or_table))) { vec_or_table = na.omit(vec_or_table); any_print (sum(is.na(vec_or_table)), 'NA are omitted from the vec_or_table of:',length(vec_or_table))}
			"Not wokring complelety : if NaN is stored as string, it does not detect it"
			}
		sum (vec_or_table==category) /  length (vec_or_table)
	} # else: is vector
} # fun

table_fixed_categories  <- function (vector, categories_vec) {
	# this function fills up the table with categories that might not occur in your vector, but are relevant.
	if ( !is.vector(vector)) {print (is(vector[]))}
	table (factor(unlist(vector), levels = categories_vec))
}

# X-react project -----------------------------------------------------------------------------------------------------

# For PNR plots
calculate_MRR <- function (vector) {
	# assume 0 = REF
	return(percentage_formatter(sum(vector == 0)/sum(vector == 1)))
}

calculate_totBR <- function (vector) {
	# assume 0 = REF
	return(percentage_formatter(sum(vector > 0.02 & vector <0.98)/length(vector)))
}

ecdf_Abel <- function (distribution, test_values=F) {
	DisLen = length(distribution)
	PCtiles = numeric(DisLen)
	for (i in 1:DisLen) {
		PCtiles[i] = sum(distribution < distribution[i]) / DisLen
	}
# 	hist(PCtiles, breaks=20)
# 	hist(distribution, breaks=20)
return(PCtiles)
}


vector_to_matrix_fillup  <- function (vector, HowManyTimes=3, IsItARow = T) {
	matt = matrix(ccc,nrow = l(ccc),ncol = HowManyTimes)
	if ( !IsItARow ) {matt = t(matt)}
	return(matt)
}
matrix_from_vector = vector_to_matrix_fillup


rowNameMatrix <- function (mat_w_dimnames) {
	matrix(rep(rownames(mat_w_dimnames), ncol(mat_w_dimnames) ),nrow = nrow(mat_w_dimnames),ncol = ncol(mat_w_dimnames))
}

colNameMatrix <- function (mat_w_dimnames) {
	x = rep(colnames(mat_w_dimnames), nrow(mat_w_dimnames) )
	t(matrix(x, nrow = ncol(mat_w_dimnames), ncol = nrow(mat_w_dimnames)))
}

simplify_categories <-   function(category_vec, replaceit , to ) {
	matches  = which(category_vec %in% replaceit); any_print(l(matches), "instances of", replaceit, "are replaced by", to)
	category_vec[matches] =  to
	return(category_vec)
}

percentile2value  <-  function(vector, percentile = 0.95, FirstValOverPercentile =T) {
	index = percentile * l(vector)
	if (FirstValOverPercentile){ index = ceiling(index)
	} else {index = floor(index) }
	value = sort(vector)[index]
	return (value)
}



# LIST FUNCTIONS ------------------------------------------------------------------------------------------------------------------------------------
as.listalike <-   function(vec, list_wannabe) {
	# convert a vector to a list with certain dimensions, taken from the list it wanna resemble
	stopifnot(length(vec) == length(unlist(list_wannabe)))
	list_return = list_wannabe
	past =0
	for (v in 1:length (list_wannabe)) {
		lv = length (list_wannabe[[v]])
		list_return[[v]] = vec[(past+1):(past+lv)]
		past = past+lv
	} # for
	list_return

}

reorder.list <- function (L, namesOrdered) {
	Lout = list(NA)
	for (x in 1:length(namesOrdered)) { Lout[[x]] = L[[namesOrdered[x] ]]  }
	if(length(names(L))) { names(Lout) = namesOrdered }
	return (Lout)
}

range.list <- function (L, namesOrdered) {
	return(range(unlist(L), na.rm=T))
}

intermingle2lists  <- function (L1, L2) {
	stopifnot(length(L1) == length(L2) )
	Lout = list(NA)
	for (x in 1:(2*length(L1)) ) {
		print (x)
		if (x  %% 2) {	Lout[[x]] = L1[[((x+1)/2)]]; names(Lout)[x] = names(L1)[((x+1)/2)]
		} else { 		Lout[[x]] = L2[[(x)/2]]; names(Lout)[x] = names(L2)[(x)/2]			}
	} # for
	return(Lout)
}

as.list.df.by.row <- function (dtf, na.omit =T, zero.omit =F, omit.empty = F) {
	# omit.empty for the listelments; na.omit and zero.omit are applied on elements inside each list elements
	outList = as.list(as.data.frame(t( dtf ) ) )
	if (na.omit){		outList =  lapply(outList, na.omit.strip)	}
	if (zero.omit){ 	outList =  lapply(outList, zero.omit) }
	if (omit.empty) { 	outList = outList[(lapply(outList, length))>0] }
	print(str(outList,vec.len = 2))
	return(outList)
}

as.list.df.by.col <- function (dtf, na.omit =T, zero.omit =F, omit.empty = F) {
	# omit.empty for the listelments; na.omit and zero.omit are applied on elements inside each list elements
	outList = as.list(dtf )
	if (na.omit){		outList =  lapply(outList, na.omit.strip)	}
	if (zero.omit){		outList =  lapply(outList, zero.omit)	}
	if (omit.empty) { 	outList = outList[(lapply(outList, length))>0] }
	print(str(outList,vec.len = 2))
	return(outList)
}

# -----------------------------------------------------------------------------------------------------

pdf.options(title= paste0('Copyright Abel Vertesy ',Sys.Date()))
# old <- `[`
# `[` <- function(...) { old(..., drop=FALSE) }



symdiff <- function(x, y, ...) { # Quasy symmetric difference of any number of vectors
	big.vec <- c(x, y, ...)
	ls = list(x, y, ...); if ( l(ls) >2) {print("# Not Mathematically correct, but logical for n>2 vectors: https://en.wikipedia.org/wiki/Symmetric_difference#Properties")}
	names(ls) = paste ("Only in", as.character(match.call()[-1]))
	duplicates <- big.vec[duplicated(big.vec)]
	lapply(ls, function (x) setdiff (x, duplicates))
}

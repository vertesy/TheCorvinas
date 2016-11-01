######################################################################
# A collection of custom R functions
######################################################################
# source ('/Users/abelvertesy/Github_repos/TheCorvinas/R/CodeAndRoll.R')
## If something is not found:
# source("/Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/Data_analysis/X_inact/Scripts_Xreact/zz_Old_versions/_Old_Functions.r")
# source("/Users/abelvertesy/MarkdownReports/MarkDownLogg.R")
## For RNA-seq specific functions, call:
# source("/Users/abelvertesy/TheCorvinas/R/RNA_seq_specific_functions.r")

### CHAPTERS:
# -  File handling, export, import [read & write]
# 	- Clipboard interaction (OS X)
# 	- Reading files in
# 	- Writing files out
# -  Vector operations
# 	- Vector filtering
# -  Matrix operations
# -  List operations
# -  Set operations
# -  Math $ stats
# -  String operations
# -  Plotting and Graphics
# -  Read and write plotting functions READ
# -  Generic
# -  Plots
# -  New additions


## Setup   -------------------------------------------------------------------------------------------------
# pdf.options(title= paste0('Copyright Abel Vertesy ',Sys.Date())) # Setup to your own name
debuggingState(on=FALSE)
l=length
oo <- function () {toClipboard(OutDir); print("OutDir is copied to the Clipbiard")}


### Load the MarkdownReports Library -------------------------------------------------------------------------------------------------
# source("/Users/abelvertesy/Github_repos/MarkdownReports/MarkdownReports/R/MarkdownReports.R")
require("MarkdownReports")
"Depends: gtools"

## File handling, export, import [read & write] -------------------------------------------------------------------------------------------------

### Clipboard interaction -------------------------------------------------------------------------------------------------
toClipboard <- function(x, sep="\t", header=FALSE, row.names=FALSE, col.names =F) { # Copy an R-object to your clipboard on OS X.
	write.table(x, pipe("pbcopy"), sep=sep, row.names=row.names, col.names =col.names, quote = F)
}

fromClipboard <- function( sep="\t", header=F) { # Paste data from your clipboard (e.g. a table from Excel) into R, parse it to a code-snippet defining an R data frame on OS X.
	return (read.table(pipe("pbpaste"), sep=sep, header=header, stringsAsFactors =F))
}

fromClipboard.as_vec <- function( sep="\t", header=F) { # Paste a list of numbers from your clipboard (e.g. from Excel) into R, parse it to a code-snippet defining an R vector on OS X.
	return (as.vector(unlist(read.table(pipe("pbpaste"), sep=sep, header=header, stringsAsFactors =F))))
}

fromClipboard.as_num_vec <- function( sep="\t", header=F) { # Paste a list of strings from your clipboard (e.g. from Excel) into R, parse it to a numeric R vector on OS X.
  return (as.numeric(unlist(read.table(pipe("pbpaste"), sep=sep, header=header, stringsAsFactors =F))))
}

fromClipboard.as_named_vec <- function( sep="\t", header=F) { # Paste a list of strings from your clipboard (e.g. from Excel) into R, parse it to a numeric R vector on OS X.
  tbl = read.table(pipe("pbpaste"), sep=sep, header=header, stringsAsFactors =F)
  vecc = tbl[ ,2]
  names(vecc) = tbl[ ,1]
  print("Names should eb in column 1, data in column 2, no header row.")
  return (vecc)
}

inline_vec.char <- function(char_vector) {	# Paste data into your code easily. Take a character vector, parse it to a code-snippet defining an R character vector, and copy back to the Clipboard.
	toClipboard(print(paste("c( '", paste (char_vector, collapse =  "', '"),  "')", collapse = "", sep=""), quote = F)); print(" Copied to Clipboard")
}

inline_vec.num <- function(num_vector) {	# Paste data into your code easily. Take a numeric vector, parse it to a code-snippet defining an R character vector, and copy back to the Clipboard.
	toClipboard(print(paste("c( ", paste (num_vector, collapse =  ", "),  " )", collapse = "", sep=""), quote = F)); print(" Copied to Clipboard")
}

inline_named_vec <- function(num_vector) {	# Paste data into your code easily. Take a numeric vector, parse it to a code-snippet defining an R character vector, and copy back to the Clipboard.
  toClipboard(    print(paste("c( ", paste (paste0('"', names(num_vector),'"'),"=", num_vector, collapse =  ", "),  " )", collapse = "", sep=""), quote = F)    )
  print(" Copied to Clipboard")
}

inline_list_char <- function(char_list) {	# Paste data into your code easily. Take a list of character vectors, parse it to a code-snippet defining an R list, and copy back to the Clipboard.
	print ("list(", quote = F)
	for (l in 1: length(list)) {
		print(paste("c( '", paste (char_list[[l]], collapse =  "', '"),  "')", collapse = "", sep=""), quote = F)
	};	print (")", quote = F)
}

inline_vec.char.from_Clipboard <- function() {	# Paste data into your code easily. Take a list of strings from your clipboard, parse it to a code-snippet defining an R character vector, and copy back to the Clipboard.
	toClipboard(print(paste("c( '", paste (fromClipboard.as_vec(), collapse =  "', '"),  "')", collapse = "", sep=""), quote = F)); print(" Copied from & to Clipboard")
}

inline_vec.num.from_Clipboard <- function() {	# Paste data into your code easily. Take a list of numbers from your clipboard, parse it to a code-snippet defining an R numeric vector, and copy back to the Clipboard.
	toClipboard(print(paste("c( ", paste (fromClipboard.as_num_vec(), collapse =  ", "),  " )", collapse = "", sep=""), quote = F)); print(" Copied from Clipboard")
}

### Reading files in -------------------------------------------------------------------------------------------------
FnP_parser <- function(fname, ext_wo_dot) { # Parses the full path from the filename & location of the file.
	if ( exists('OutDir') ) { path = OutDir } else { path = getwd() ; any_print ("OutDir not defined !!!") }
	if (hasArg(ext_wo_dot) ) { FnP = kollapse (path,"/", fname, ".", ext_wo_dot)
	} else { 					FnP = kollapse (path,"/", fname) }
}

read.simple.vec <- function(...) {  # Read each line of a file to an element of a vector (read in new-line separated values, no header!).
	pfn = kollapse (...) # merge path and filename
	read_in = as.vector(unlist(read.table( pfn , stringsAsFactors=F, sep = "\n" )) )
	any_print(length (read_in), "elements")
	return(read_in);
}

read.simple <- function(...) { # It is essentially read.table() with file/path parsing.
	pfn = kollapse (...) # merge path and filename
	read_in = read.table( pfn , stringsAsFactors=F)
	return(read_in)
}

read.simple_char_list <- function(...) { # Read in a file.
	pfn = kollapse (...) # merge path and filename
	read_in = unlist(read.table( pfn , stringsAsFactors=F ) )
	any_print ("New variable head: ",what(read_in))
	return(read_in)
}

read.simple.table <- function(..., colnames=T ) { # Read in a file. default: header defines colnames, no rownames. For rownames give the col nr. with rownames, eg. 1 The header should start with a TAB / First column name should be empty.
	pfn = kollapse (...) # merge path and filename
	read_in = read.table( pfn , stringsAsFactors=FALSE, sep="\t", header=colnames )
	any_print ("New variable dim: ", dim(read_in))
	return(read_in)
}

read.simple.tsv <- function(...) { # Read in a file with excel style data: rownames in col1, headers SHIFTED. The header should start with a TAB / First column name should be empty.
	pfn = kollapse (...) # merge path and filename
	read_in = read.delim( pfn , stringsAsFactors=FALSE, sep="\t", row.names=1, header=T )
	any_print ("New variable dim: ",dim(read_in))
	return(read_in)
}

read.simple.tsv.named.vector <- function(...) { # Read in a file with excel style named vectors, names in col1, headers SHIFTED. The header should start with a TAB / First column name should be empty.
	pfn = kollapse (...) # merge path and filename
	read_in = read.delim( pfn , stringsAsFactors=FALSE, sep="\t", row.names=1, header=T )
	rn = row.names(read_in)
	read_in =  as.vector(unlist(read_in));
	names(read_in) = rn
	any_print ("New vectors length is: ",length(read_in))
	return(read_in)
}

### Writing files out -------------------------------------------------------------------------------------------------

write.simple <- function(input_df, extension='tsv', ManualName ="", o = F,...  ) { # Write out a matrix-like R-object to a file with as tab separated values (.tsv). Your output filename will be either the variable's name. The output file will be located in "OutDir" specified by you at the beginning of the script, or under your current working directory. You can pass the PATH and VARIABLE separately (in order), they will be concatenated to the filename.
	fname = kollapse (...) ; if (nchar (fname) < 2 ) { fname = substitute(input_vec) }
	if (nchar(ManualName)) {FnP = kollapse(ManualName)} else  { FnP = FnP_parser (fname, extension) }
	write.table (input_df, file = FnP, sep = "\t", row.names = F, col.names = T, quote=FALSE  )
	if (o) { system(paste0("open ", FnP), wait = F) }
	any_print ("Length: ", length(input_df))
} # fun

write.simple.vec <- function(input_vec, extension='vec', ManualName ="", o = F, ... ) { # Write out a vector-like R-object to a file with as newline separated values (.vec). Your output filename will be either the variable's name. The output file will be located in "OutDir" specified by you at the beginning of the script, or under your current working directory. You can pass the PATH and VARIABLE separately (in order), they will be concatenated to the filename.
	fname = kollapse (...) ; if (nchar (fname) < 2 ) { fname = substitute(input_vec) }
	if (nchar(ManualName)) {FnP = kollapse(ManualName)} else  { FnP = FnP_parser (fname, extension) }
	write.table (input_vec, file = FnP, sep = "\t", row.names = F, col.names = F, quote=FALSE  )
	any_print ("Length: ", length(input_vec))
	if (o) { system(paste0("open ", FnP), wait = F) }
} # fun

write.simple.tsv <- function(input_df, extension='tsv', ManualName ="", o = F, ... ) { # Write out a matrix-like R-object WITH ROW- AND COLUMN- NAMES to a file with as tab separated values (.tsv). Your output filename will be either the variable's name. The output file will be located in "OutDir" specified by you at the beginning of the script, or under your current working directory. You can pass the PATH and VARIABLE separately (in order), they will be concatenated to the filename.
	fname = kollapse (..., print = F); if (nchar (fname) < 2 ) { fname = substitute(input_df) }
	if (nchar(ManualName)) {FnP = kollapse(ManualName)} else  { FnP = FnP_parser (fname, extension) }
	write.table (input_df, file = FnP, sep = "\t", row.names = T, col.names = NA, quote=FALSE  )
	printme = if(l(dim(input_df))) paste0("Dim: ", dim(input_df) ) else paste0("Length (of your vector): ", l(input_df) )
	any_print (printme)
	if (o) { system(paste0("open ", FnP), wait = F) }
} # fun
# If col.names = NA and row.names = TRUE a blank column name is added, which is the convention used for CSV files to be read by spreadsheets.

write.simple.append <- function(input_df, extension='tsv', ManualName ="", o = F, ... ) { # Append an R-object WITHOUT ROWNAMES, to an existing .tsv file of the same number of columns. Your output filename will be either the variable's name. The output file will be located in "OutDir" specified by you at the beginning of the script, or under your current working directory. You can pass the PATH and VARIABLE separately (in order), they will be concatenated to the filename.
	fname = kollapse (...) ; if (nchar (fname) < 2 ) { fname = substitute(input_df) }
	if (nchar(ManualName)) { FnP = kollapse(ManualName)} else  { FnP = FnP_parser (fname, extension) }
	write.table (input_df, file = FnP, sep = "\t", row.names = F,col.names = F, quote=FALSE, append=T  )
	if (o) { system(paste0("open ", FnP), wait = F) }
} # fun


## Vector operations -------------------------------------------------------------------------------------------------

as.named.vector <- function(df_col, WhichDimNames = 1) { # Convert a dataframe column or row into a vector, keeping the corresponding dimension name.
	# use RowNames: WhichDimNames = 1 , 2: use ColNames
	# !!! might require drop=F in subsetting!!! eg: df_col[,3, drop=F]
	# df_col[which(unlist(lapply(df_col, is.null)))] = "NULL" # replace NULLs - they would fall out of vectors - DOES not work yet
	namez = dimnames(df_col)[[WhichDimNames]]
	if (is.list(df_col) & !is.data.frame(df_col)) {namez = names(df_col)}
	vecc = as.vector(unlist (df_col))
	names (vecc)= namez
	return (vecc)
}

as.numeric.wNames <- function(vec) { # Converts any vector into a numeric vector, and puts the original character values into the names of the new vector, unless it already has names. Useful for coloring a plot by categories, name-tags, etc.
	numerified_vec = as.numeric(as.factor(vec))
	if (!is.null(names(vec))) {names (numerified_vec) = names (vec)}
	return(numerified_vec)
}

as.logical.wNames <- function(vec) { # Converts your input vector into a logical vector, and puts the original character values into the names of the new vector, unless it already has names.
	numerified_vec = as.logical(vec)
	if (!is.null(names(vec))) {names (numerified_vec) = names (vec)}
	return(numerified_vec)
}

as.character.wNames <- function(vec) { # Converts your input vector into a character vector, and puts the original character values into the names of the new vector, unless it already has names.
	char_vec = as.character(vec)
	if (!is.null(names(vec))) {names (char_vec) = names (vec)}
	return(char_vec)
}

rescale <- function(vec, from=0, upto=100) { # linear transformation to a given range of values
	vec = vec-min(vec, na.rm = T)
	vec = vec*((upto-from)/max(vec, na.rm = T))
	vec = vec+ from
	return (vec)
} # fun

flip_value2name <- function(named_vector, NumericNames =F) { # Flip the values and the names of a vector with names
  if (! is.null(names(named_vector))) {
    newvec = names(named_vector)
    if (NumericNames) {  newvec = as.numeric(names(named_vector))     }
    names(newvec) = named_vector
  } else {llprint("Vector without names!", head(named_vector))}
  if (any(duplicated(named_vector))) {llprint("New names contain duplicated elements",head(named_vector[which(duplicated(named_vector))]))  }
  if (any(duplicated(newvec))) {llprint("Old names contained duplicated elements",head(newvec[which(duplicated(newvec))]))  }
  return(newvec)
}

value2name_flip = flip_value2name
# sortbyitsnames <- function(vec_or_list) { # Sort a vector by the alphanumeric order of its names (instead of its values).
# 	print("THIS FUCNTION MAKES MISTAKES WITH DUPLICATE NAMES")
# 	if (is.vector(vec_or_list) & !is.list(vec_or_list)) {  vec[gtools::mixedsort(names(vec_or_list) )]
# 	} else if (is.list(vec_or_list)) {	reorder.list(L = (vec_or_list), namesOrdered = gtools::mixedsort(names(vec_or_list))) }
# 	}


sortbyitsnames <- function(vec_or_list) { # Sort a vector by the alphanumeric order of its names (instead of its values).
	xx = names(vec_or_list)
	names(xx) = 1:l(vec_or_list)
	order = as.numeric(names(gtools::mixedsort(xx)))
	vec_or_list[order]
}

### Vector filtering  -------------------------------------------------------------------------------------------------

which_names <- function(named_Vec) { # Return the names where the input vector is TRUE. The input vector is converted to logical.
	return(names(which(as.logical.wNames(named_Vec)))) }


na.omit.strip <- function(vec) {  # Omit NA values from a vector and return a clean vector without any spam.
	if (is.data.frame(vec)) {
		if ( min(dim(vec)) > 1 ) { any_print(dim(vec), "dimensional array is converted to a vector.") }
		vec = unlist(vec) }
	clean = na.omit(vec)
	attributes(clean)$na.action <- NULL
	return(clean)
}

inf.omit <- function(vec) { # Omit infinite values from a vector.
	if (is.data.frame(vec)) {
		if ( min(dim(vec)) > 1 ) { any_print(dim(vec), "dimensional array is converted to a vector.") }
		vec = unlist(vec) }
	clean = vec[is.finite(vec)]
	# attributes(clean)$na.action <- NULL
	return(clean)
}

zero.omit <- function(vec) { # Omit zero values from a vector.
	v2= vec[vec!=0]
	any_print("range: ", range(v2))
	if ( !is.null(names(vec)) ) {names(v2) = names(vec)[vec!=0]}
	return(v2)
}

pc_TRUE <- function(logical_vector, percentify =T) { # Percentage of true values in a logical vector, parsed as text (useful for reports.)
	out = sum(logical_vector, na.rm=T) / length(logical_vector)
	if (percentify) {out = percentage_formatter (out) }
	return(out)
	}

pc_in_total_of_match <- function(vec_or_table, category, NA_omit=T) { # Percentage of a certain value within a vector or table.
	if (is.table(vec_or_table)) { vec_or_table[category]/sum(vec_or_table, na.rm=NA_omit) }
	else { # if (is.vector(vec_or_table))
		if (NA_omit){
			if (sum(is.na(vec_or_table))) { vec_or_table = na.omit(vec_or_table); any_print (sum(is.na(vec_or_table)), 'NA are omitted from the vec_or_table of:',length(vec_or_table))}
			"Not wokring complelety : if NaN is stored as string, it does not detect it"
			}
		sum (vec_or_table==category) /  length (vec_or_table)
	} # else: is vector
} # fun

filter_survival_length <- function(length_new, length_old, prepend ="") { # Parse a sentence reporting the % of filter survival.
	pc = percentage_formatter(length_new/length_old)
	llprint (prepend, pc, " of ",length_old," entries make through the filter")
}

remove_outliers <- function(x, na.rm = TRUE, ..., probs = c(.05, .95)) { # Remove values that fall outside the trailing N % of the distribution.
	qnt <- quantile(x, probs=probs, na.rm = na.rm, ...)
	H <- 1.5 * IQR(x, na.rm = na.rm)
	y <- x
	y[x < (qnt[1] - H)] <- NA
	y[x > (qnt[2] + H)] <- NA
	y
}

simplify_categories <-  function(category_vec, replaceit , to ) { # Replace every entry that is found in "replaceit", by a single value provided by "to"
	matches  = which(category_vec %in% replaceit); any_print(l(matches), "instances of", replaceit, "are replaced by", to)
	category_vec[matches] =  to
	return(category_vec)
}

## Matrix operations -------------------------------------------------------------------------------------------------
colSort <- function(data, ...) { # Sort each column of a numeric matrix / data frame.
	sapply(data, sort, ...) }

colMedians <- function(mat,na.rm=TRUE) { # Calculates the median of each column of a numeric matrix / data frame.
	return(apply(mat,2,median, na.rm=na.rm)) }

rowMedians <- function(mat,na.rm=TRUE) { # Calculates the median of each row of a numeric matrix / data frame.
	return(apply(mat,1,median, na.rm=na.rm)) }

colCV <- function(mat) return(apply(mat,2,cv ) ) # Calculates the CV of each column of a numeric matrix / data frame.

rowMin <- function(x) apply(x, 1, min) # Calculates the minimum of each row of a numeric matrix / data frame.

rowMax <- function(x) apply(x, 1, max) # Calculates the maximum of each row of a numeric matrix / data frame.
	
colMin <- function(x) apply(x, 2, min) # Calculates the minimum of each column of a numeric matrix / data frame.

colMax <- function(X) apply(x, 2, max) # Calculates the maximum of each column of a numeric matrix / data frame.

sort.mat <- function(df, colname_in_df = 1, decrease = F, na_last = T) { # Sort a matrix. ALTERNATIVE: dd[with(dd, order(-z, b)), ]. Source: https://stackoverflow.com/questions/1296646/how-to-sort-a-dataframe-by-columns-in-r
	if (length(colname_in_df)>1) { print ("cannot handle multi column sort") }
	else {df[ order(df[,colname_in_df], decreasing = decrease, na.last = na_last), ]}
}

rowNameMatrix <- function(mat_w_dimnames) {  # Create a copy of your matrix, where every entry is replaced by the corresponding row name. Useful if you want to color by row name in a plot (where you have different number of NA-values in each row).
	matrix(rep(rownames(mat_w_dimnames), ncol(mat_w_dimnames) ),nrow = nrow(mat_w_dimnames),ncol = ncol(mat_w_dimnames))
}

colNameMatrix <- function(mat_w_dimnames) { # Create a copy of your matrix, where every entry is replaced by the corresponding column name. Useful if you want to color by column name in a plot (where you have different number of NA-values in each column).
	x = rep(colnames(mat_w_dimnames), nrow(mat_w_dimnames) )
	t(matrix(x, nrow = ncol(mat_w_dimnames), ncol = nrow(mat_w_dimnames)))
}

matrix_from_vector <- function(vector, HowManyTimes=3, IsItARow = T) { # Create a matrix from values in a vector repeated for each column / each row. Similar to rowNameMatrix and colNameMatrix.
	matt = matrix(vector,nrow = l(vector),ncol = HowManyTimes)
	if ( !IsItARow ) {matt = t(matt)}
	return(matt)
}

matrix_from_dimnames  <- function (rownames, colnames, fill = NA) { # Set up an empty matrix from two vectors that will define row- and column-names.
  mm = matrix(data = fill, nrow = length(rownames), ncol = length(colnames), dimnames = list(rownames, colnames))
}

colsplit <- function(df, f) { # split a data frame by a factor corresponding to columns.
  ListOfDFs = NULL
  levelz = unique(f)
  for (i in 1:l(levelz)) {    ListOfDFs[[i]] = df[ , which(f== levelz[i]) ]  }
  return(ListOfDFs)
}
splitByCol = colsplit


median_normalize <- function(mat) { # normalize each column to the median of the columns
  cs = colSums(mat, na.rm = T)
  norm_mat = (t(t(mat) / cs)) * median(cs)
}

## List operations -------------------------------------------------------------------------------------------------
list.wNames <- function(...){ # create a list with names from ALL variables you pass on to the function
	l = list(...)
	names(l) = as.character(match.call()[-1])
	return(l)
}

as.list.df.by.row <- function(dtf, na.omit =T, zero.omit =F, omit.empty = F) { # Split a dataframe into a list by its columns. omit.empty for the listelments; na.omit and zero.omit are applied on entries inside each list element.
	outList = as.list(as.data.frame(t( dtf ) ) )
	if (na.omit){		outList =  lapply(outList, na.omit.strip)	}
	if (zero.omit){ 	outList =  lapply(outList, zero.omit) }
	if (omit.empty) { 	outList = outList[(lapply(outList, length))>0] }
	print(str(outList,vec.len = 2))
	return(outList)
}

as.list.df.by.col <- function(dtf, na.omit =T, zero.omit =F, omit.empty = F) { # oSplit a dataframe into a list by its rows. omit.empty for the listelments; na.omit and zero.omit are applied on entries inside each list element.
	outList = as.list(dtf )
	if (na.omit){		outList =  lapply(outList, na.omit.strip)	}
	if (zero.omit){		outList =  lapply(outList, zero.omit)	}
	if (omit.empty) { 	outList = outList[(lapply(outList, length))>0] }
	print(str(outList,vec.len = 2))
	return(outList)
}

reorder.list <- function(L, namesOrdered) { # reorder elements of lists in your custom order of names / indices.
	Lout = list(NA)
	for (x in 1:length(namesOrdered)) { Lout[[x]] = L[[namesOrdered[x] ]]  }
	if(length(names(L))) { names(Lout) = namesOrdered }
	return (Lout)
}

range.list <- function(L, namesOrdered) { # range of values in whole list
	return(range(unlist(L), na.rm=T))
}

intermingle2lists <- function(L1, L2) { # Combine 2 lists (of the same length) so that form every odd and every even element of a unified list. Useful for side-by-side comparisons, e.g. in wstripchart_list().
	stopifnot(length(L1) == length(L2) )
	Lout = list(NA)
	for (x in 1:(2*length(L1)) ) {
		if (x  %% 2) {	Lout[[x]] = L1[[((x+1)/2)]]; names(Lout)[x] = names(L1)[((x+1)/2)]
		} else { 		Lout[[x]] = L2[[(x)/2]]; names(Lout)[x] = names(L2)[(x)/2]			}
	} # for
	return(Lout)
}

intermingle2vec <- function(V1, V2) { # Combine 2 vectors (of the same length) so that form every odd and every even element of a unified vector.
  stopifnot(length(V1) == length(V2) )
  LEN = length(c(V1, V2))
  
  Vout = rep(NA, LEN)
  Vout[seq(1,LEN, by=2)] = V1
  Vout[seq(2, LEN, by=2)] = V2
  return(Vout)
}

as.listalike <-  function(vec, list_wannabe) { # convert a vector to a list with certain dimensions, taken from the list it wanna resemble
	stopifnot(length(vec) == length(unlist(list_wannabe)))
	list_return = list_wannabe
	past =0
	for (v in 1:length (list_wannabe)) {
		lv = length (list_wannabe[[v]])
		list_return[[v]] = vec[(past+1):(past+lv)]
		past = past+lv
	} # for
	return(list_return)
}

list2df_presence <- function(yalist, entries_list = F, matrixfill = "") { # Convert a list to a full dataframe, summarizing the presence or absence of elements
  if( is.null(names(yalist)) ) {names(yalist) = 1:length(yalist)}
  
  rown = unique(unlist(yalist))
  coln =  names(yalist)
  mm = matrix_from_dimnames(rown, coln, fill = matrixfill)
  entries_list = lapply(yalist, names)
  
  for (i in 1:length(yalist)) {
    print(i)
    le = unlist(yalist[i])
    names(le) = unlist(entries_list[i])
    
    list_index = which( le  %in% rown)
    m_index = which( rown %in% le)
    mm[ m_index,i] = names(le[list_index])
  }
  return(mm)  
}



list_to_fullDF <- function(ll){ # convert a list to a full numeric data matrix. Designed for occurence counting, think tof table()
  entrytypes = unique(names(unlist(ll)))
  mat =matrix(0, ncol = l(ll), nrow = l(entrytypes))
  colnames(mat) = 1:l(ll);   rownames(mat) = sort(entrytypes)
  for (i in 1:l(ll)) {
    entries = ll[[i]]
    mat[names(entries) ,i] = entries
  }
  mat
}


## Set operations -------------------------------------------------------------------------------------------------

symdiff <- function(x, y, ...) { # Quasy symmetric difference of any number of vectors
	big.vec <- c(x, y, ...)
	ls = list(x, y, ...); if ( l(ls) >2) {print("# Not Mathematically correct, but logical for n>2 vectors: https://en.wikipedia.org/wiki/Symmetric_difference#Properties")}
	names(ls) = paste ("Only in", as.character(match.call()[-1]))
	duplicates <- big.vec[duplicated(big.vec)]
	lapply(ls, function(x) setdiff (x, duplicates))
}

## Math $ stats -------------------------------------------------------------------------------------------------

sem <- function(x) {  # Calculates the standard error of the mean (SEM) for a numeric vector (it excludes NA-s by default)
	sd(x)/sqrt(length(x))}

cv <- function(x) { # Calculates the coefficient of variation (CV) for a numeric vector (it excludes NA-s by default)
	sd( x, na.rm=T)/mean(x, na.rm=T) }

fano <- function(x) { # Calculates the fano factor on a numeric vector (it excludes NA-s by default)
	var(x, na.rm=T)/mean(x, na.rm=T) }

modus <- function(x) { # Calculates the modus of a numeric vector (it excludes NA-s by default)
	x= unlist(na.exclude(x))
	ux <- unique(x)
	tab <- tabulate(match(x, ux));
	ux[tab == max(tab)]
}

gm_mean <- function(x, na.rm=TRUE){ # Calculates the geometric mean of a numeric vector (it excludes NA-s by default)
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) }

mean_of_log <- function(x, k=2, na.rm=TRUE){ # Calculates the mean of the log_k of a numeric vector (it excludes NA-s by default)
  negs = sum(x<0);  zeros = sum(x==0)
  if (negs | zeros) { any_print("The input vector has", negs, "negative values and", zeros, "zeros." )  }
  mean(log(x, base = k), na.rm = na.rm) }

movingAve <- function(x, oneSide) { # Calculates the moving / rolling average of a numeric vector.
	y = NULL
	for (i in oneSide:l(x)) {
		y[i] = mean( x[ (i-oneSide):(i+oneSide) ] )
	}; 	return (y)
}

movingSEM <- function(x, oneSide) { # Calculates the moving / rolling standard deviation of the mean (SEM) on a numeric vector.
	y = NULL
	for (i in oneSide:l(x)) {
		y[i] = sem( x[ (i-oneSide):(i+oneSide) ] )
	}; 	return (y)
}

imovingSEM <- function(x, oneSide = 5) { # Calculates the moving / rolling standard deviation of the mean (SEM). It calculates it to the edge of the vector with incrementally smaller window-size.
	y = NULL
	for (i in 1:l(x)) {
		oneSideDynamic = min(i-1,oneSide, l(x)-i); oneSideDynamic
		indexx = (i-oneSideDynamic):(i+oneSideDynamic);indexx
		y[i] = sem( x[ indexx ] )
	}; 	return (y)
}

## String operations  -------------------------------------------------------------------------------------------------
eval_parse_kollapse <- function( ... ){ # evaluate and parse (dyn_var_caller)
	substitute(eval(parse(text=kollapse( ... , print=F))))
}

substrRight <- function(x, n){ # Take the right substring of a string
	substr(x, nchar(x)-n+1, nchar(x))
}

lookup <- function(needle, haystack, exact =TRUE, report = FALSE) { # Awesome pattern matching for a set of values in another set of values. Returns a list with all kinds of results.
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
}

translate = replace_values <- function(vec, oldvalues, newvalues) { # Replaces a set of values in a vector with another set of values, it translates your vector. Oldvalues and newvalues have to be 1-to-1 corespoding vectors.
  Nr = l(oldvalues)
  if (Nr > l(newvalues) ) {
    if (l(newvalues) == 1) {
      newvalues =  rep(newvalues, l(oldvalues))
    } else if (l(newvalues) > 1) { any_print("PROVIDE ONE NEWVALUE, OR THE SAME NUMEBR OF NEWVALUES AS OLDVALUES.")}
  }
  tmp = vec
  for (i in 1:Nr) {
    oldval = oldvalues[i]
    tmp[vec==oldval] = newvalues[i]
  }
  return(tmp)
}
'chartr("a-cX", "D-Fw", x) does the same as above in theory, but it did not seem very robust regarding your input...'

capitalize_Firstletter <- function(s, strict = FALSE) { # Capitalize every first letter of a word
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

## Plotting and Graphics -----------------------------------------------------------------------------------------------------

HeatMapCol_BGR <- grDevices::colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
# HeatMapCol_BWR <- grDevices::colorRampPalette(c("blue", "white", "red"), bias=1)
HeatMapCol_RedBlackGreen <- grDevices::colorRampPalette(c("red", "black", "green"), bias=1)

lm_equation_formatter <- function(lm) { # Renders the lm() function's output into a human readable text. (e.g. for subtitles)
	eq = (lm$coefficients);
	kollapse ("Intercept:", eq[1], " Slope:", eq[2]);
}

## Read and write plotting functions READ -------------------------------------
# rw-funcitons cannot call the w-functions, there would be too much things lost: rw writes where the file comes from, w writes in Outdir & current dir

rwplot <- function(FnP, ..., w=7, h=7) { # read in a table, plot it and save it.
	print ('file without header, and file path between ""')
	variable=read.simple  (FnP);
	fname=gsub (".*/", "",FnP);  trunk = gsub ("\\..*", "",fname)
	plot (variable, ..., main=trunk )
	dev.copy2pdf (file=kollapse(FnP,".plot.pdf"), width=w, height=h )
}

rwscatterplot <- function(FnP, ..., w=7, h=7) { # read in a table, plot it as scatter plot and save it.
	print ('file without header, and file path between ""')
	variable=read.simple.table  (FnP);
	fname=gsub (".*/", "",FnP);  trunk = gsub ("\\..*", "",fname)
	line = lm(variable[,2]~variable[,1])
	plot (variable, ..., main=trunk, sub = lm_equation_formatter (line) )
	abline(line)
	dev.copy2pdf (file=kollapse(FnP,".scatter.pdf"), width=w, height=h )
}

rwboxplot <- function(FnP, col ="gold1", ..., w=7, h=7) { # read in a table, plot it and save it.
	print ('inputfile is a dataframe, with header')
	variable=read.simple.table (FnP);
	fname=gsub (".*/", "",FnP);  trunk = gsub ("\\..*", "",fname)
	boxplot (variable, ..., main=trunk, col =col, las=2)
	dev.copy2pdf (file=kollapse(FnP,".boxplot.pdf"), width=w, height=h )
}

rwhist <- function(FnP, col ="gold1", ..., w=7, h=7) { # read in a table, plot it and save it. File should be without header, and file path should be given between ""-s.
	print ('file without header, and file path between ""')
	variable=read.simple (FnP);
	fname=gsub (".*/", "",FnP); trunk = gsub ("\\..*", "",fname)
	if ( length (variable) > 0 ) {
		if ( !is.numeric(variable)) { variable = table (variable) ; wbarplot (variable); print ("FILENAME IS: VARIABLE.BARPLOT.PDF") }
		else { 	hist (variable, ..., plotname=trunk, col =col, las=2)
				dev.copy2pdf (file=kollapse(FnP,".hist.pdf"), width=w, height=h )
		} # if is.numeric
	} # if non empty
}

rwbarplot <- function(FnP, col ="gold1", ..., w=7, h=7) { # read in a vector, plot it and save it. File should be without header, and file path should be given between ""-s.
	print ('file without header, and file path between ""')
	variable=read.simple (FnP);
	fname=gsub (".*/", "",FnP);  trunk = gsub ("\\..*", "",fname)
	barplot (variable, ..., main=trunk, col =col, las=2)
	dev.copy2pdf (file=kollapse(FnP,".barplot.pdf"), width=w, height=h )
}

## Plots -------------------------------------------------------------------------------------------------


pdfA4plot_on <- function (pname = date(), ..., w = 8.27, h = 11.69, rows = 4, cols = 3, mdlink = FALSE,
                          title = paste0(basename(fname), " by ", if (exists("scriptname")) scriptname else "Rscript")) { # Print (multiple) plots to an (A4) pdf.
  assign("mfrow_default", par("mfrow"), fname, envir = .GlobalEnv)
  fname = kollapse(OutDir,"/" , pname, ".pdf")
  pdf(fname,width=w, height=h, title = title)
  par(mfrow = c(rows, cols))
  any_print(" ----  Don't forget to call the pair of this function to finish the pdf: pdfA4plot_off ()")
  if (mdlink) { MarkDown_Img_Logger_PDF_and_PNG(fname_wo_ext = plotname) }
}

pdfA4plot_off <- function () {
  if (exists("mfrow_default")) {
    x = mfrow_default
  } else { x =  c(1,1)}
  par(mfrow = x)
  try(dev.off()) # close pdf
  oo();
}

## Generic -------------------------------------------------------------------------------------------------


stopif <- function(condition, message ="") { if(condition) {any_print (message); stop()} } # Stop script if the condition is met

most_frequent_elements <- function(thingy, topN=10) { # Show the most frequent elements of a table
	tail(sort(table(thingy, useNA = "ifany")), topN)
}

what <- function(x, printme=0) { # A better version of is(). It can print the first "printme" elements.
	any_print (is (x),"; nr. of elements:", length (x))
	if (is.numeric (x) ) 		{ any_print ("min&max:", range(x) ) } else {print ("Not numeric")}
	if ( length(dim(x) ) > 0 ) 	{ any_print ("Dim:", dim (x) )	}
	if ( printme>0) 			{ any_print ("Elements:", x[0:printme] )	}
	head (x)
}

idim <- function(any_object) { # A dim() function that can handle if you pass on a vector: then, it gives the length.
	if (is.null(dim(any_object))) {print(length(any_object))}
	else {	print(dim(any_object))	}
}

idimnames <- function(any_object) { # A dimnames() function that can handle if you pass on a vector: it gives back the names.
	if (!is.null(dimnames(any_object))) 	{ print(dimnames(any_object)) }
	else if (!is.null(colnames(any_object))) { any_print("colnames:", colnames(any_object))	}
	else if (!is.null(rownames(any_object))) { any_print("rownames:", rownames(any_object))	}
	else if (!is.null(names(any_object))) { any_print("names:", names(any_object))	}
}

table_fixed_categories <- function(vector, categories_vec) { # generate a table() with a fixed set of categories. It fills up the table with missing categories, that are relevant when comparing to other vectors.
	if ( !is.vector(vector)) {print (is(vector[]))}
	table (factor(unlist(vector), levels = categories_vec))
}

percentile2value <- function(distribution, percentile = 0.95, FirstValOverPercentile =T) { # Calculate what is the actual value of the N-th percentile in a distribution or set of numbers. Useful for calculating cutoffs, and displaying them by whist()'s "vline" paramter.
	index = percentile * l(distribution)
	if (FirstValOverPercentile){ index = ceiling(index)
	} else {index = floor(index) }
	value = sort(distribution)[index]
	return (value)
}

top_indices <- function(x, n = 3, top = T){ # Returns the position / index of the n highest values. For equal values, it maintains the original order
	head( order(x, decreasing =top), n )
}

attach_w_rownames <- function(df_w_dimnames) { # Take a data frame (of e.g. metadata) from your memory space, split it into vectors so you can directly use them. E.g.: Instead of metadata$color[blabla] use color[blabla]
	if(!is.null(rownames(df_w_dimnames)) & !is.null(colnames(df_w_dimnames))) {
		namez= rownames(df_w_dimnames)
		any_print("Now directly available in the workspace:      ", colnames(df_w_dimnames))
		attach (df_w_dimnames)
		for (n in colnames(df_w_dimnames)) {
			x=get(n); names(x) = namez
			assign (n,x,envir =.GlobalEnv) } # for
	} else { print ("ERROR: the DF does not have some of the dimnames!")}
}

Color_Check <- function(..., incrBottMarginBy=0 ) { # Display the colors encoded by the numbers / color-ID-s you pass on to this function
	if (incrBottMarginBy) { .ParMarDefault <- par("mar"); 	par(mar=c(par("mar")[1]+incrBottMarginBy, par("mar")[2:4]) ) } 	# Tune the margin
	Numbers  = c(...)
	if (l(names(Numbers)) == l(Numbers)) {labelz = names(Numbers)} else {labelz = Numbers}
	barplot (rep(10,length(Numbers)), col = Numbers, names.arg = labelz, las=2 )
	if (incrBottMarginBy) { par("mar" = .ParMarDefault )}
}

## New additions -----------------------------------------------------------------------------------------------------


cormethod = "spearman"
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, method = cormethod) { # A function to display correlation values for pairs() function. Default is pearson correlation, that can be set to  "kendall" or "spearman".
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y, method = method)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}


any.duplicated <- function (vec, summarize=T){ # How many entries are duplicated
  y=sum(duplicated(vec))
  if(summarize & y){
    x = table(vec); x= x[x>1]-1;
    print("The following elements have >1 extra copies:")
    print(x) # table formatting requires a separate entry
  }
  return(y)
}


merge_numeric_df_by_rn <-function(x, y) { # Merge 2 numeric data frames by rownames
  merged =  merge(x ,y, by="row.names", all=TRUE)  # merge by row names (by=0 or by="row.names")
  rownames(merged) = merged$Row.names
  merged = merged[ ,-1] # remove row names
  merged[is.na(merged)] <- 0  
  return(merged)
}

# CodeAndRoll - A collection of custom R functions

[Get CodeAndRoll](https://github.com/vertesy/TheCorvinas/blob/master/R/CodeAndRoll.R)


### Install

1.) Download and source: https://github.com/vertesy/TheCorvinas/blob/master/R/CodeAndRoll.R


2.) Source from the web:  

`source("https://raw.githubusercontent.com/vertesy/TheCorvinas/master/R/CodeAndRoll.R")`

Read: https://github.com/vertesy/TheCorvinas/blob/master/R/CodeAndRoll.md



## Usage





### If you encounter a bug, or something doesn't work, please let me know by raising an issue on [The Corvinas](https://github.com/vertesy/TheCorvinas/issues/new?milestone=CodeAndRoll)




## CHAPTERS

- ### File handling, export, import [read & write]
	- ### Clipboard interaction (OS X)
	- ### Reading files in
	- ### Writing files out
- ### Vector operations
	- ### Vector filtering
- ### Matrix operations
- ### List operations
- ### Set operations
- ### Math $ stats
- ### String operations
- ### Plotting and Graphics
- ### Read and write plotting functions READ
- ### Generic
- ### Plots
- ### New additions



## Index

|  Function 	| Description |
|---|---|
| try.dev.off 	|  |
| topN.dfCol 	|  |
| bottomN.dfCol 	|  |
| colSums.barplot 	|  |
| as.factor.numeric 	|  |
| sstrsplit 	|  |
| coolor 	|  |
| unlapply 	| lapply, then unlist |
| oo 	| function () {toClipboard(OutDir); print("OutDir is copied to the Clipbiard")} |
| toClipboard 	| Copy an R-object to your clipboard on OS X. |
| fromClipboard 	| Paste data from your clipboard (e.g. a table from Excel) into R, parse it to a code-snippet defining an  |
| fromClipboard.as_vec 	| Paste a list of numbers from your clipboard (e.g. from Excel) into R, parse it to a code-snippet  |
| fromClipboard.as_num_vec 	| Paste a list of strings from your clipboard (e.g. from Excel) into R, parse it to a numeric  |
| fromClipboard.as_named_vec 	| Paste a list of strings from your clipboard (e.g. from Excel) into R, parse it to a  |
| inline_vec.char 	| Paste data into your code easily. Take a character vector, parse it to a code-snippet defining an R  |
| inline_vec.num 	| Paste data into your code easily. Take a numeric vector, parse it to a code-snippet defining an R character  |
| inline_named_vec 	| Paste data into your code easily. Take a numeric vector, parse it to a code-snippet defining an R character  |
| inline_list_char 	| Paste data into your code easily. Take a list of character vectors, parse it to a code-snippet defining an R  |
| inline_vec.char.from_Clipboard 	| Paste data into your code easily. Take a list of strings from your clipboard, parse it to a code-snippet  |
| inline_vec.num.from_Clipboard 	| Paste data into your code easily. Take a list of numbers from your clipboard, parse it to a code-snippet  |
| FnP_parser 	| Parses the full path from the filename & location of the file. |
| read.simple.vec 	| Read each line of a file to an element of a vector (read in new-line separated values, no header!). |
| read.simple 	| It is essentially read.table() with file/path parsing. |
| read.simple_char_list 	| Read in a file. |
| read.simple.table 	| Read in a file. default: header defines colnames, no rownames. For rownames give the col nr. with  |
| read.simple.tsv 	| Read in a file with excel style data: rownames in col1, headers SHIFTED. The header should start with a TAB / First  |
| read.simple.tsv.named.vector 	| Read in a file with excel style named vectors, names in col1, headers SHIFTED. The header should start  |
| write.simple 	| Write out a matrix-like R-object to a file with as tab separated  |
| write.simple.vec 	| Write out a vector-like R-object to a file with as newline  |
| write.simple.tsv 	| Write out a matrix-like R-object WITH ROW- AND COLUMN- NAMES to a  |
| write.simple.append 	| Append an R-object WITHOUT ROWNAMES, to an existing .tsv file  |
| as.named.vector 	| Convert a dataframe column or row into a vector, keeping the corresponding dimension name. |
| as.numeric.wNames 	| Converts any vector into a numeric vector, and puts the original character values into the names of the new vector,  |
| as.logical.wNames 	| Converts your input vector into a logical vector, and puts the original character values into the names of the new  |
| as.character.wNames 	| Converts your input vector into a character vector, and puts the original character values into the names of the  |
| rescale 	| linear transformation to a given range of values |
| flip_value2name 	| Flip the values and the names of a vector with names |
| # sortbyitsnames 	| Sort a vector by the alphanumeric order of its names (instead of its values). |
| sortbyitsnames 	| Sort a vector by the alphanumeric order of its names (instead of its values). |
| which_names 	| Return the names where the input vector is TRUE. The input vector is converted to logical. |
| na.omit.strip 	| Omit NA values from a vector and return a clean vector without any spam. |
| inf.omit 	| Omit infinite values from a vector. |
| zero.omit 	| Omit zero values from a vector. |
| pc_TRUE 	| Percentage of true values in a logical vector, parsed as text (useful for reports.) |
| pc_in_total_of_match 	| Percentage of a certain value within a vector or table. |
| filter_survival_length 	| Parse a sentence reporting the % of filter survival. |
| remove_outliers 	| Remove values that fall outside the trailing N % of the distribution. |
| simplify_categories 	| Replace every entry that is found in "replaceit", by a single value provided by "to" |
| sortEachColumn 	| Sort each column of a numeric matrix / data frame. |
| colMedians 	| Calculates the median of each column of a numeric matrix / data frame. |
| rowMedians 	| Calculates the median of each row of a numeric matrix / data frame. |
| colCV 	| Calculates the CV of each column of a numeric matrix / data frame. |
| rowMin 	| Calculates the minimum of each row of a numeric matrix / data frame. |
| rowMax 	| Calculates the maximum of each row of a numeric matrix / data frame. |
| colMin 	| Calculates the minimum of each column of a numeric matrix / data frame. |
| colMax 	| Calculates the maximum of each column of a numeric matrix / data frame. |
| sort.mat 	| Sort a matrix. ALTERNATIVE: dd[with(dd, order(-z, b)), ]. Source: https:// |
| rowNameMatrix 	| Create a copy of your matrix, where every entry is replaced by the corresponding row name. Useful if you  |
| colNameMatrix 	| Create a copy of your matrix, where every entry is replaced by the corresponding column name. Useful if you  |
| matrix_from_vector 	| Create a matrix from values in a vector repeated for each column / each row.  |
| matrix_from_dimnames  	| Set up an empty matrix from two vectors that will define row- and column-names. |
| colsplit 	| split a data frame by a factor corresponding to columns. |
| median_normalize 	| normalize each column to the median of the columns |
| list.wNames 	| create a list with names from ALL variables you pass on to the function |
| as.list.df.by.row 	| Split a dataframe into a list by its columns. omit.empty for the  |
| as.list.df.by.col 	| oSplit a dataframe into a list by its rows. omit.empty for the  |
| reorder.list 	| reorder elements of lists in your custom order of names / indices. |
| range.list 	| range of values in whole list |
| intermingle2lists 	| Combine 2 lists (of the same length) so that form every odd and every even element of a unified list. Useful for  |
| intermingle2vec 	| Combine 2 vectors (of the same length) so that form every odd and every even element of a unified vector. |
| as.listalike 	| convert a vector to a list with certain dimensions, taken from the list it wanna resemble |
| list2df_presence 	| Convert a list to a full dataframe, summarizing the presence or absence of  |
| list_to_fullDF 	| convert a list to a full numeric data matrix. Designed for occurence counting, think tof table() |
| splitbyitsnames 	| split a list by its names |
| splititsnames_byValues 	| split a list by its names |
| symdiff 	| Quasy symmetric difference of any number of vectors |
| sem 	| Calculates the standard error of the mean (SEM) for  |
| cv 	| Calculates the coefficient of variation (CV) for a numeric vector (it excludes NA-s  |
| fano 	| Calculates the fano factor on a numeric vector (it excludes NA-s by default) |
| modus 	| Calculates the modus of a numeric vector (it excludes NA-s by default) |
| gm_mean 	| Calculates the geometric mean of a numeric vector (it excludes NA-s by default) |
| mean_of_log 	| Calculates the mean of the log_k of a numeric vector (it excludes NA-s by default) |
| movingAve 	| Calculates the moving / rolling average of a numeric vector. |
| movingSEM 	| Calculates the moving / rolling standard deviation of the mean (SEM) on a numeric vector. |
| imovingSEM 	| Calculates the moving / rolling standard deviation of the mean (SEM). It calculates it to the edge of the  |
| eval_parse_kollapse 	| evaluate and parse (dyn_var_caller) |
| substrRight 	| Take the right substring of a string |
| lookup 	| Awesome pattern matching for a set of values in another set of values. Returns a  |
| translate = replace_values 	| Replaces a set of values in a vector with another set of values, it translates your  |
| capitalize_Firstletter 	| Capitalize every first letter of a word |
| cap 	| Renders the lm() function's output into a human readable text. (e.g. for subtitles) |
| lm_equation_formatter 	| read in a table, plot it and save it. |
| rwplot 	| read in a table, plot it as scatter plot and save it. |
| rwscatterplot 	| read in a table, plot it and save it. |
| rwboxplot 	| read in a table, plot it and save it. File should be without header, and file path should be  |
| rwhist 	| read in a vector, plot it and save it. File should be without header, and file path should be  |
| rwbarplot 	|  |
| pdfA4plot_on 	|  |
| pdfA4plot_off 	|  |
| stopif 	| Stop script if the condition is met |
| most_frequent_elements 	| Show the most frequent elements of a table |
| what 	| A better version of is(). It can print the first "printme" elements. |
| idim 	| A dim() function that can handle if you pass on a vector: then, it gives the length. |
| idimnames 	| A dimnames() function that can handle if you pass on a vector: it gives back the names. |
| table_fixed_categories 	| generate a table() with a fixed set of categories. It fills up the table with missing  |
| percentile2value 	| Calculate what is the actual value of the N-th percentile in  |
| top_indices 	| Returns the position / index of the n highest values. For equal values, it maintains the original order |
| attach_w_rownames 	| Take a data frame (of e.g. metadata) from your memory space, split it into  |
| Color_Check 	| Display the colors encoded by the numbers / color-ID-s you pass on to this function |
| panel.cor 	| A function to display correlation values for pairs() function. Default  |
| any.duplicated 	| How many entries are duplicated |
| plot_filtering_RaceID 	| function(sc, minexpr=p$minexpr, minnumber = p$minnumber) { |
| printEveryN 	| Report at every e.g. 1000 |
| icolor_categories 	| function (vec, rndize=F) {  x= table(vec);colvec = coolor(l(x)); if(rndize) colvec=sample(colvec); names(colvec) =names(x); return( |
| wlegend2 	| Add a legend, and save the plot  |
| # llwrite_list 	| function(yalist) { |
| irequire 	| install package if  |
| zigzagger 	| mix entries so that they differ |
| name2id 	|  |
| wvenn 	|  |


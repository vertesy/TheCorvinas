
# CodeAndRoll - A collection of custom R functions
[Get CodeAndRoll](https://github.com/vertesy/CodeAndRoll/blob/master/CodeAndRoll.R)

> **NOTE CodeAndRoll has been moved to its [own repository](https://github.com/vertesy/CodeAndRoll)!**

## Install

1.) [Download `CodeAndRoll.R`](https://github.com/vertesy/CodeAndRoll/blob/master/CodeAndRoll.R), save as local `.R` file, and `source()`: 

2.) Source from the web:

`source("https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R")`

### Troubleshooting

*If you encounter a **bug**, or something doesn't work, please let me know by raising an issue on [CodeAndRoll](https://github.com/vertesy/CodeAndRoll/issues)*



## Chapters 

The script is roughly organised in the following sections / categories.

1. ##### File handling, export, import [read & write]
  1. ##### Clipboard interaction (OS X)
  1. ##### Reading files in
  1. ##### Writing files out
1. ##### Vector operations
  1. ##### Vector filtering
1. ##### Matrix operations
1. ##### List operations
1. ##### Set operations
1. ##### Math and stats
1. ##### String operations
1. ##### Plotting and Graphics
1. ##### Read and write plotting functions READ
1. ##### Generic
1. ##### Plots
1. ##### New additions


## Index

| Function  	| Description |
|---|---|
| grepv  	|  grep returning the value |
| read.simple.vec  	|  Read each line of a file to an element of a vector (read in new-line separated values, no header!). |
| read.simple  	|  It is essentially read.table() with file/path parsing. |
| read.simple_char_list  	|  Read in a file. |
| read.simple.table  	|  Read in a file. default: header defines colnames, no rownames. For rownames give the col nr. with rownames, eg. 1 The header should start with a TAB / First column name should be empty. |
| read.simple.tsv  	|  Read in a file with excel style data: rownames in col1, headers SHIFTED. The header should start with a TAB / First column name should be empty. |
| read.simple.tsv.named.vector  	|  Read in a file with excel style named vectors, names in col1, headers SHIFTED. The header should start with a TAB / First column name should be empty. |
| read.simple.xls  	|  Read multi-sheet excel files. row_namePos = NULL for automatic names |
| write.simple  	|  Write out a matrix-like R-object to a file with as tab separated values (.tsv). Your output filename will be either the variable's name. The output file will be located in "OutDir" specified by you at the beginning of the script, or under your current working directory. You can pass the PATH and VARIABLE separately (in order), they will be concatenated to the filename. |
| write.simple.vec  	|  Write out a vector-like R-object to a file with as newline separated values (.vec). Your output filename will be either the variable's name. The output file will be located in "OutDir" specified by you at the beginning of the script, or under your current working directory. You can pass the PATH and VARIABLE separately (in order), they will be concatenated to the filename. |
| write.simple.tsv  	|  Write out a matrix-like R-object WITH ROW- AND COLUMN- NAMES to a file with as tab separated values (.tsv). Your output filename will be either the variable's name. The output file will be located in "OutDir" specified by you at the beginning of the script, or under your current working directory. You can pass the PATH and VARIABLE separately (in order), they will be concatenated to the filename. |
| write.simple.append  	|  Append an R-object WITHOUT ROWNAMES, to an existing .tsv file of the same number of columns. Your output filename will be either the variable's name. The output file will be located in "OutDir" specified by you at the beginning of the script, or under your current working directory. You can pass the PATH and VARIABLE separately (in order), they will be concatenated to the filename. |
| as.factor.numeric  	|  Turn any vector into numeric categories as.numeric(as.factor(vec)) |
| sstrsplit  	|  Alias for str_split_fixed in the stringr package |
| topN.dfCol  	|  Find the n highest values in a named vector |
| bottomN.dfCol  	|  Find the n lowest values in a named vector |
| as.named.vector  	|  Convert a dataframe column or row into a vector, keeping the corresponding dimension name. |
| col2named.vector  	|  Convert a dataframe column into a vector, keeping the corresponding dimension name. |
| row2named.vector  	|  Convert a dataframe row into a vector, keeping the corresponding dimension name. |
| as.numeric.wNames  	|  Converts any vector into a numeric vector, and puts the original character values into the names of the new vector, unless it already has names. Useful for coloring a plot by categories, name-tags, etc. |
| as.numeric.wNames.old  	|  Converts any vector into a numeric vector, and puts the original character values into the names of the new vector, unless it already has names. Useful for coloring a plot by categories, name-tags, etc. |
| as.logical.wNames  	|  Converts your input vector into a logical vector, and puts the original character values into the names of the new vector, unless it already has names. |
| as.character.wNames  	|  Converts your input vector into a character vector, and puts the original character values into the names of the new vector, unless it already has names. |
| rescale  	|  linear transformation to a given range of values |
| flip_value2name  	|  Flip the values and the names of a vector with names |
| # sortbyitsnames  	|  Sort a vector by the alphanumeric order of its names (instead of its values). |
| sortbyitsnames  	|  Sort a vector by the alphanumeric order of its names (instead of its values). |
| any.duplicated  	|  How many entries are duplicated |
| unique.wNames  	|  Get the unique elements from a vector keeping names (of the first elements of each category). |
| which_names  	|  Return the names where the input vector is TRUE. The input vector is converted to logical. |
| na.omit.strip  	|  Omit NA values from a vector and return a clean vector without any spam. |
| na.omit.mat  	|  Omit rows with NA values from a matrix. Rows with any, or full of NA-s |
| inf.omit  	|  Omit infinite values from a vector. |
| zero.omit  	|  Omit zero values from a vector. |
| pc_TRUE  	|  Percentage of true values in a logical vector, parsed as text (useful for reports.) |
| pc_in_total_of_match  	|  Percentage of a certain value within a vector or table. |
| filter_survival_length  	|  Parse a sentence reporting the % of filter survival. |
| remove_outliers  	|  Remove values that fall outside the trailing N % of the distribution. |
| rotate   	|  rotate a matrix 90 degrees. |
| sortEachColumn  	|  Sort each column of a numeric matrix / data frame. |
| rowMedians  	|  Calculates the median of each row of a numeric matrix / data frame. |
| colMedians  	|  Calculates the median of each column of a numeric matrix / data frame. |
| rowGeoMeans  	|  Calculates the median of each row of a numeric matrix / data frame. |
| colGeoMeans  	|  Calculates the median of each column of a numeric matrix / data frame. |
| colCV  	|  Calculates the CV of each column of a numeric matrix / data frame. |
| rowCV  	|  Calculates the CV of each column of a numeric matrix / data frame. |
| rowMin  	|  Calculates the minimum of each row of a numeric matrix / data frame. |
| colMin  	|  Calculates the minimum of each column of a numeric matrix / data frame. |
| rowMax  	|  Calculates the maximum of each row of a numeric matrix / data frame. |
| colMax  	|  Calculates the maximum of each column of a numeric matrix / data frame. |
| rowSEM  	|  Calculates the SEM of each row of a numeric matrix / data frame. |
| colSEM  	|  Calculates the SEM of each column of a numeric matrix / data frame. |
| colDivide  	|  divide by column |
| rowDivide  	|  divide by row |
| sort.mat  	|  Sort a matrix. ALTERNATIVE: dd[with(dd, order(-z, b)), ]. Source: https://stackoverflow.com/questions/1296646/how-to-sort-a-dataframe-by-columns-in-r |
| rowNameMatrix  	|  Create a copy of your matrix, where every entry is replaced by the corresponding row name. Useful if you want to color by row name in a plot (where you have different number of NA-values in each row). |
| colNameMatrix  	|  Create a copy of your matrix, where every entry is replaced by the corresponding column name. Useful if you want to color by column name in a plot (where you have different number of NA-values in each column). |
| matrix_from_vector  	|  Create a matrix from values in a vector repeated for each column / each row. Similar to rowNameMatrix and colNameMatrix. |
| matrix_from_dimnames   	|  Set up an empty matrix from two vectors that will define row- and column-names. |
| colsplit  	|  split a data frame by a factor corresponding to columns. |
| median_normalize  	|  normalize each column to the median of all the columns |
| mean_normalize  	|  normalize each column to the median of the columns |
| rownames.trimws  	|  trim whitespaces from the rownames |
| select.rows.and.columns  	|  Subset rows and columns. It checks if the selected dimension names exist and reports if any of those they aren't found. |
| getRows  	|  Get the subset of rows with existing rownames, report how much it could not find. |
| get.oddoreven  	|  Get odd or even columns or rows of a data frame |
| combine.matrices.intersect  	|  combine matrices by rownames intersect |
| merge_numeric_df_by_rn  	|  Merge 2 numeric data frames by rownames |
| attach_w_rownames  	|  Take a data frame (of e.g. metadata) from your memory space, split it into vectors so you can directly use them. E.g.: Instead of metadata$color[blabla] use color[blabla] |
| panel.cor.pearson  	|  A function to display correlation values for pairs() function. Default is pearson correlation, that can be set to  "kendall" or "spearman". |
| panel.cor.spearman  	|  A function to display correlation values for pairs() function. Default is pearson correlation, that can be set to  "kendall" or "spearman". |
| remove.na.rows  	|  cols have to be a vector of numbers corresponding to columns |
| intersect.ls  	|  Intersect any number of lists. |
| unlapply  	|  lapply, then unlist |
| list.wNames  	|  create a list with names from ALL variables you pass on to the function |
| as.list.df.by.row  	|  Split a dataframe into a list by its columns. omit.empty for the listelments; na.omit and zero.omit are applied on entries inside each list element. |
| as.list.df.by.col  	|  oSplit a dataframe into a list by its rows. omit.empty for the listelments; na.omit and zero.omit are applied on entries inside each list element. |
| reorder.list  	|  reorder elements of lists in your custom order of names / indices. |
| range.list  	|  range of values in whole list |
| intermingle2lists  	|  Combine 2 lists (of the same length) so that form every odd and every even element of a unified list. Useful for side-by-side comparisons, e.g. in wstripchart_list(). |
| list2df_presence  	|  Convert a list to a full dataframe, summarizing the presence or absence of elements |
| list2fullDF  	|  convert a list to a full numeric data matrix. Designed for occurence counting, think tof table() |
| splitbyitsnames  	|  split a list by its names |
| splititsnames_byValues  	|  split a list by its names |
| intermingle2vec  	|  Combine 2 vectors (of the same length) so that form every odd and every even element of a unified vector. |
| intermingle.cbind  	|  Combine 2 data frames (of the same length) so that form every odd and every even element of a unified list. Useful for side-by-side comparisons, e.g. in wstripchart_list(). |
| pad.na  	|  Fill up with a vector to a given length with NA-values at the end. |
| list2df_NA_padded  	|  When converting a list to a data frame, the list elements can have different lengths. This function fills up the data frame with NA values. |
| clip.values  	|  Signal clipping. Cut values above or below a threshold. |
| list2df  	|  Basic list-to-df functionality in R |
| symdiff  	|  Quasy symmetric difference of any number of vectors |
| sem  	|  Calculates the standard error of the mean (SEM) for a numeric vector (it excludes NA-s by default) |
| cv  	|  Calculates the coefficient of variation (CV) for a numeric vector (it excludes NA-s by default) |
| fano  	|  Calculates the fano factor on a numeric vector (it excludes NA-s by default) |
| modus  	|  Calculates the modus of a numeric vector (it excludes NA-s by default) |
| geomean  	|  Calculates the geometric mean of a numeric vector (it excludes NA-s by default) |
| mean_of_log  	|  Calculates the mean of the log_k of a numeric vector (it excludes NA-s by default) |
| movingAve  	|  Calculates the moving / rolling average of a numeric vector. |
| movingSEM  	|  Calculates the moving / rolling standard error of the mean (SEM) on a numeric vector. |
| imovingSEM  	|  Calculates the moving / rolling standard error of the mean (SEM). It calculates it to the edge of the vector with incrementally smaller window-size. |
| eval_parse_kollapse  	|  evaluate and parse (dyn_var_caller) |
| substrRight  	|  Take the right substring of a string |
| lookup  	|  Awesome pattern matching for a set of values in another set of values. Returns a list with all kinds of results. |
| translate = replace_values  	|  Replaces a set of values in a vector with another set of values, it translates your vector. Oldvalues and newvalues have to be 1-to-1 corespoding vectors. |
| richColors  	|  Alias for rich.colors in gplots |
| Color_Check  	|  Display the colors encoded by the numbers / color-ID-s you pass on to this function |
| colSums.barplot  	|  Draw a barplot from ColSums of a matrix. |
| lm_equation_formatter  	|  Renders the lm() function's output into a human readable text. (e.g. for subtitles) |
| hist.XbyY  	|  Split a one variable by another. Calculates equal bins in splitby, and returns a list of the corresponding values in toSplit. |
| nameiftrue  	|  returns the name if its value is true |
| flag.name_value  	|  returns the name if its value is true |
| quantile_breaks  	|  Quantile breakpoints in any data vector http://slowkow.com/notes/heatmap-tutorial/ |
| vec.fromNames  	|  create a vector from a vector of names |
| list.fromNames  	|  create list from a vector with the names of the elements |
| matrix.fromNames  	|  create a matrix from 2 vectors defining the row- and column names of the matrix |
| what  	|  A better version of is(). It can print the first "printme" elements. |
| idim  	|  A dim() function that can handle if you pass on a vector: then, it gives the length. |
| idimnames  	|  A dimnames() function that can handle if you pass on a vector: it gives back the names. |
| table_fixed_categories  	|  generate a table() with a fixed set of categories. It fills up the table with missing categories, that are relevant when comparing to other vectors. |
| stopif  	|  Stop script if the condition is met |
| stopif2  	|  Stop script if the condition is met. You can parse anything (e.g. variables) in the message |
| most_frequent_elements  	|  Show the most frequent elements of a table |
| top_indices  	|  Returns the position / index of the n highest values. For equal values, it maintains the original order |
| percentile2value  	|  Calculate what is the actual value of the N-th percentile in a distribution or set of numbers. Useful for calculating cutoffs, and displaying them by whist()'s "vline" paramter. |
| hclust.getOrder.row  	|  Extract ROW order from a pheatmap object. |
| hclust.getOrder.col  	|  Extract COLUMN order from a pheatmap object. |
| hclust.getClusterID.row  	|  Extract cluster ID's for ROWS of a pheatmap object. |
| hclust.getClusterID.col  	|  Extract cluster ID's for COLUMNS of a pheatmap object. |
| hclust.ClusterSeparatingLines.row  	|  Calculate the position of ROW separating lines between clusters in a pheatmap object. |
| hclust.ClusterSeparatingLines.col  	|  Calculate the position of COLUMN separating lines between clusters in a pheatmap object. |
| Gap.Postions.calc.pheatmap  	|  calculate gap positions for pheatmap, based a sorted annotation vector of categories |
| matlabColors.pheatmap  	|  Create a Matlab-like color gradient using "colorRamps". |
| annot_col.create.pheatmap.vec  	|  For VECTORS. Auxiliary function for pheatmap. Prepares the 2 variables needed for "annotation_col" and "annotation_colors" in pheatmap |
| annot_col.create.pheatmap.df  	|  For data frames. Auxiliary function for pheatmap. Prepares the 2 variables needed for "annotation_col" and "annotation_colors" in pheatmap |
| annot_col.fix.numeric  	|  fix class and color annotation in pheatmap annotation data frame's and lists. |
| wPairConnector  	|  Connect Pairs of datapoints with a line on a plot. |
| numerate  	|  numerate from x to y with additonal zeropadding |
| printEveryN  	|  Report at every e.g. 1000 |
| zigzagger  	|  mix entries so that they differ |
| irequire  	|  Load a package. If it does not exist, try to install it from CRAN. |
| NrAndPc  	|  Summary stat. text formatting for logical vectors (%, length) |

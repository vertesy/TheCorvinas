## _Old_Functions that are kept for reference!
# source("~/Github_repos/TheCorvinas/R/Old_Functions.r")


# ------------------------
# ------------------------
# ------------------------
# ------------------------
# ------------------------

rich.colors.vec <- function(vec, pre=0, post=0, randomize=F, seed=11) { # Generates a vector of colors from a vector of categories with rich.colors() for a numeric vector
  colz = gplots::rich.colors(pre+l(unique(vec))+post)
  colz = colz[(pre+1):( l(colz)-post  )] # subset
  if (randomize) {
    set.seed(seed)
    colz = sample(colz)
  }
  names(colz) = unique(vec)
  if (randomize) Color_Check(colz)# , ylab=p0("seed:", seed)
  colz[as.character(vec)] # convert to character because they are referred by name
}

# icolor_categories <- function (vec, rndize=F) {  x= table(vec);colvec = richColors(l(x)); if(rndize) colvec=sample(colvec); names(colvec) =names(x); return(colvec) } # create color categories
icolor_categories <- function (vec, rndize=F, trail=0, seed=354, plotit=F) {
  x= table(vec);
  colvec = richColors(l(x)+trail );
  colvec = colvec[ trail:l(x)+trail];
  if(rndize) {set.seed(seed); colvec=sample(colvec)};
  names(colvec) =names(x);
  if(plotit) Color_Check(colvec);
  return(colvec)
} # create color categories # ; colvec = colvec[ (1+trail):(l(x)+trail)]


# ------------------------

# This is a legacy function
#' wlegend.old
#' Add a legend, and save the plot immediately
#'
#' @param x location of legend
#' @param legend Labels displayed (Text)
#' @param fill Color of the boxes next to the text
#' @param bty Background of legend, transparent by default
#' @param OverwritePrevPDF Save the plot immediately with the same name the last wplot* function made (It is stored in plotnameLastPlot variable).
#' @param ... Pass any other parameter of the corresponding text function (most of them should work).
#' @examples wlegend(...)
#' @export


wlegend.old <- function(x=c("topleft", "topright", "bottomright", "bottomleft")[4],
                        legend, fill = NULL, ..., bty = "n", OverwritePrevPDF =T) { # Add a legend, and save the plot immediately
  legend(x=x, legend=legend, fill=fill, ..., bty=bty)
  if (OverwritePrevPDF) {   wplot_save_this(plotname = plotnameLastPlot)  }
}

# ------------------------

# old <- `[`
'If the rownames are retained, the subselected column and row remains the original type of the object, eg df. hist () does not take df-s for instance'
# `[` <- function(...) { old(..., drop=FALSE) }
'Convenience in the original [ function, objects are simplified df[,col] -> vector.'
# `[` <- old


capitalize_Firstletter <- function(s, strict = FALSE) { # Capitalize every first letter of a word
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


normalize_per_column <- function (m, factor) { # Normalize each column by its by a corresponding value. It is solely t(t(matrix)/factor) with a couple of error checks.
	any_print("Range of normalization factors:" , range(iround(factor)))
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



# 15 August 2017 (Tuesday) 13:05  -------------------------------------------------------

#' #' qlegend
#' #' # Quickly add a legend, and save the plot immediately
#' #'
#' #' @param NamedColorVec A vector defining the Color of the boxes, while its names define the labels.
#' #' @param poz 1:4 corresponding to "topleft","topright", "bottomright", "bottomleft"
#' #' @param ... Additional parameters for legend()
#' #' @param w_ Width of the saved pdf image, in inches.
#' #' @param h_ Height of the saved pdf image, in inches.
#' #' @param bty The type of box to be drawn around the legend. The allowed values are "o" (the default) and "n".
#' #' @param OverwritePrevPDF Save the plot immediately with the same name the last wplot* function made (It is stored in plotnameLastPlot variable).
#' #' @export
#' #'
#' #' @examples ccc = as.factor.numeric(c(6,7,8,7,7,6,6,8))  ; qlegend(ccc) # Uses as.factor.numeric() from Github / Vertesy / TheCorvinas
#'
#' qlegend <- function(NamedColorVec, poz=3, ..., w_=7, h_=w_, bty = "n", OverwritePrevPDF =T) {
#'   pozz = translate(poz, oldvalues = 1:4, newvalues = c("topleft","topright", "bottomright", "bottomleft"))
#'   fill_ = getCategories(NamedColorVec)
#'   legend(x=pozz, legend=names(fill_), fill=fill_, ..., bty=bty)
#'   if (OverwritePrevPDF) {   wplot_save_this(plotname = plotnameLastPlot, w= w_, h = h_)  }
#' }
#'


# cormethod = "spearman"
# panel.cor <- function(x, y, digits=2, prefix="", cex.cor, method = cormethod) { # A function to display correlation values for pairs() function. Default is pearson correlation, that can be set to  "kendall" or "spearman".
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   r <- abs(cor(x, y, method = method))
#   txt <- format(c(r, 0.123456789), digits=digits)[1]
#   txt <- paste(prefix, txt, sep="")
#   if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
#
#   test <- cor.test(x,y)
#   # borrowed from printCoefmat
#   Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
#                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
#                    symbols = c("***", "**", "*", ".", " "))
#
#   text(0.5, 0.5, txt, cex = cex * r)
#   text(.8, .8, Signif, cex=cex, col=2)
# }


# http://stackoverflow.com/questions/20127282/r-color-scatterplot-points-by-col-value-with-legend
# scatter_fill <- function (x, y, color, xlim=range(x), ylim=range(y), zlim=range(color),
#                           nlevels = 20, plot.title, plot.axes, pch=21, cex=1,
#                           key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
#                           axes = TRUE, frame.plot = axes, ...) {
#   mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
#   on.exit(par(par.orig))
#   w <- (3 + mar.orig[2L]) * par("csi") * 2.54
#   layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
#   par(las = las)
#   mar <- mar.orig
#   mar[4L] <- mar[2L]
#   mar[2L] <- 1
#   par(mar = mar)
#
#   # choose colors to interpolate
#   levels <- seq(zlim[1], zlim[2], length.out = nlevels)
#   col <- colorRampPalette(c("red", "yellow", "dark green"))(nlevels)
#   colz <- col[cut(color, nlevels)]
#
#   plot.new()
#   plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
#
#   rect(0, levels[-length(levels)], 1, levels[-1L], col=col, border=col)
#   if (missing(key.axes)) {if (axes){axis(4)}}
#   else key.axes
#   box()
#   if (!missing(key.title))
#     key.title
#   mar <- mar.orig
#   mar[4L] <- 1
#   par(mar = mar)
#
#   # points
#   plot(x, y, type = "n", xaxt='n', yaxt='n', xlab="", ylab="", xlim=xlim, ylim=ylim, bty="n")
#   points(x, y, bg = colz, xaxt='n', yaxt='n', xlab="", ylab="", bty="n", pch=pch,...)
#
#   ## options to make mapping more customizable
#   if (missing(plot.axes)) {
#     if (axes) {
#       title(main = "", xlab = "", ylab = "")
#       Axis(x, side = 1)
#       Axis(y, side = 2)
#     }
#   }
#   else plot.axes
#   if (frame.plot)
#     box()
#   if (missing(plot.title))
#     title(...)
#   else plot.title
#   invisible()
# }


# annot_col.create <- function(df, annot_vec, annot_names=NULL) { # Auxiliary function for pheatmap. Prepares the 2 variables needed for "annotation_col" and "annotation_colors" in pheatmap
#   stopifnot( l(annot_vec) == dim(df)[2] )
#   print(substitute(annot_vec))
#   df = as.data.frame(annot_vec)
#   if (!is.null(annot_names)) {  stopifnot(length(annot_vec) == length(annot_names));
#     colnames(df) = annot_names
#   } else {    colnames(df) =  substitute(annot_vec)    }
#
#   df[,1] = as.character(df[,1])
#   assign(x = "annot", value = df, envir = .GlobalEnv)
#
#   xx = list(annot_vec = val2col(annot_vec[!duplicated(annot_vec)]))
#   if (!is.null(annot_names)) {    stopifnot(length(annot_col) == length(annot_names));
#     names(xx) = annot_names  }
#   else {    names(xx) = substitute(annot_vec)  }
#
#   assign(x = "annot_col", value = xx, envir = .GlobalEnv)
#   print("annot and annot_col variables are created. Use: pheatmap(..., annotation_col = annot, annotation_colors = annot_col)")
# }
# # annot_col.create

# JUNKYARD --------------------------------------------------------------------------------

# select.rows <- function(df, RowIndices ) {
#   true_rownames = intersect(rownames(df), RowIndices)
#   NotFound = setdiff(RowIndices, rownames(df))
#   if (l(NotFound)) { any_print(l(NotFound), "Row Indices Not Found:", head(NotFound), "...     Rows found:", l(true_rownames))  } #if
#   return(df[ true_rownames, ])
# }


# require("corrr")
# require("dplyr")
# require(ggplot2)
# corrr::correlate(mtcars[, 1:4], method = "pearson") %>%
#   corrr::rearrange(absolute=F) %>% #  rearrange cols and rows
#   corrr::shave()%>% # shave upper triangle
#   corrr::rplot(print_cor = T,legend = T)# dot plot
#
# corrr::correlate(mtcars[, 1:4], method = "pearson") %>% network_plot()


# ## A largish data set
# n <- 10000
# x1  <- matrix(rnorm(n), ncol = 2)
# x2  <- matrix(rnorm(n, mean = 3, sd = 1.5), ncol = 2)
# x   <- rbind(x1, x2)
# oldpar <- par(mfrow = c(2, 2))
# smoothScatter(x, nrpoints = 0)
# dim(x)
# l(x)


# # http://stackoverflow.com/questions/14271584/r-legend-for-color-density-scatterplot-produced-using-smoothscatter
# install.packages("fields")
# fudgeit <- function(){
#   xm <- get('xm', envir = parent.frame(1))
#   ym <- get('ym', envir = parent.frame(1))
#   z  <- get('dens', envir = parent.frame(1))
#   colramp <- get('colramp', parent.frame(1))
#   fields::image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F)
# }
# #
# par(mar = c(5,4,4,5) + .1)
# smoothScatter(x, nrpoints = 0, postPlotHook = fudgeit)


# wvenn <- function (yalist, imagetype = "png", alpha = .5, fill = 1:length(yalist), ..., w = 7, h = 7, mdlink = F, plotname = substitute(yalist)) {
#   if (!require("VennDiagram")) { print("Please install VennDiagram: install.packages('VennDiagram')") }
#   fname = kollapse(plotname, ".", imagetype, print = F)
#   LsLen = length(yalist)
#   if(length(names(yalist)) < LsLen) { names(yalist) =1:LsLen; print("List elements had no names.") }
#   print(names(yalist))
#
#   filename = kollapse(OutDir,"/", fname, print = F)
#   subt = kollapse("Total = ", length(unique(unlist(yalist))), " elements in total.", print = F)
#   venn.diagram(x = yalist, imagetype = imagetype, filename = filename, main = plotname, ... ,
#                sub = subt, fill = fill, alpha = alpha, sub.cex = .75, main.cex = 2)
#   if (mdlink) {
#     llogit(MarkDown_ImgLink_formatter(fname))
#     if (exists("png4Github") & png4Github == T) { llogit(MarkDown_ImgLink_formatter(paste0("Reports/", fname) ) )	}
#   }
# }

# wLinRegression <- function(DF, coeff = c("pearson", "spearman", "r2")[3], textlocation = "topleft", savefile =T, ...) { # Add linear regression, and descriptors to line to your scatter plot. Provide the same dataframe as you provided to wplot() before you called this function
#   print(coeff)
#   regression <- lm(DF[,2] ~ DF[,1])
#   abline(regression, ...)
#   legendText = NULL
#   if ( "pearson" %in% coeff) {    dispCoeff = iround(cor(DF[,2], DF[,1], method = "pearson"))
#   legendText  =  c(legendText, paste0("Pearson c.c.: ", dispCoeff))  }
#   if ("spearman" %in% coeff) {    dispCoeff = iround(cor(DF[,2], DF[,1], method = "spearman"))
#   legendText = c(legendText, paste0("Spearman c.c.: ", dispCoeff))  }
#   if ("r2" %in% coeff) {          r2 = iround(summary(regression)$r.squared)
#   legendText = c(legendText, paste0("R^2: ", r2))  }
#   print(legendText)
#   if (length(coeff)==1 & "r2" == coeff[1]) {  legend(textlocation, legend = superscript_in_plots(prefix = "R", sup = "2",suffix = paste0(": ", r2)) , bty="n")
#   } else {                                    legend(textlocation, legend = legendText , bty="n") }
#   if(savefile){   wplot_save_this(plotname = plotnameLastPlot) }
# }




#' #' wpie
#' #'
#' #' Create and save pie charts as .pdf, in "OutDir". If mdlink =T, it inserts a .pdf and a .png link in the markdown report, set by "path_of_report". The .png version is not created, only the link is put in place, not to overwrite previous versions.
#' #' @param variable The variable to plot.
#' #' @param ... Pass any other parameter of the corresponding plotting function (most of them should work).
#' #' @param percentage Display percentage instead of counts. TRUE by default.
#' #' @param both_pc_and_value Report both percentage AND number.
#' #' @param plotname Title of the plot (main parameter) and also the name of the file.
#' #' @param col Fill color. Defined by rich colours by default
#' #' @param savefile Save plot as pdf in OutDir, TRUE by default.
#' #' @param w Width of the saved pdf image, in inches.
#' #' @param h Height of the saved pdf image, in inches.
#' #' @param mdlink Insert a .pdf and a .png image link in the markdown report, set by "path_of_report".
#' #' @examples wpie (variable =  , ... =  , percentage = TRUE, plotname = substitute(variable), w = 7, h = 7, mdlink = F)
#' #' @export
#'
#' wpie <-function (variable, ..., percentage = TRUE, both_pc_and_value=F, plotname = substitute(variable), col = gplots::rich.colors(length(variable)), savefile = T, w = 7, h = 7, mdlink = F) {
#'   if (!require("gplots")) { print("Please install gplots: install.packages('gplots')") }
#'   fname = kollapse(plotname, ".pie")
#'   subt = kollapse("Total = ", sum(variable), print = F)
#'   if (percentage) {	labs <- paste("(", names(variable), ")", "\n", percentage_formatter(variable/sum(variable)), sep = "")
#'   if (both_pc_and_value) { labs <- paste("(", names(variable), ")", "\n", percentage_formatter(variable/sum(variable)),"\n", variable , sep = "")}
#'   } else {	labs <- paste("(", names(variable), ")", "\n", variable, sep = "")	}
#'   pie(variable, ..., main = plotname, sub = subt, clockwise = T, labels = labs, col = col )
#'   if (savefile) { dev.copy2pdf(file = FnP_parser(fname, "pdf"), width = w, height = h, title = ttl_field()) }
#'   if (mdlink) { MarkDown_Img_Logger_PDF_and_PNG(fname_wo_ext = fname) }
#' }
#'


# wvioplot_list <-function (yalist, ..., coll = c(2:(length(yalist)+1)),
#                           plotname = as.character(substitute(yalist)), sub = NULL, xlb = names(yalist), ylb = "", ylimm=F,
#                           incrBottMarginBy = 0, tilted_text = F, yoffset=0, savefile = T, w = 7, h = 7, mdlink = F) {
#   if (!require("vioplot")) { print("Please install vioplot: install.packages('vioplot')") }
#   if (incrBottMarginBy) { .ParMarDefault <- par("mar"); 	par(mar=c(par("mar")[1]+incrBottMarginBy, par("mar")[2:4]) ) } 	# Tune the margin
#   l_list = length(yalist)
#   fname = kollapse(plotname, ".vioplot")
#   if (length(coll) < l_list) { coll = rep(coll, l_list) }
#   if (tilted_text) {	xlb = NA } else { xlb = names(yalist) }
#   if (! (is.numeric(ylimm) & length(ylimm)==2)) { ylimm = range(unlist(yalist),na.rm = T)}
#   plot(0, 0, type = "n", xlim = c(0.5, (l_list + 0.5)), ylim = ylimm, xaxt = "n", xlab = "",
#        ylab = ylb, main = plotname, sub = sub)
#   for (i in 1:l_list) {
#     if( l(na.omit.strip(yalist[[i]])) ){
#       vioplot::vioplot(na.omit(yalist[[i]]), ..., at = i, add = T, col = coll[i])
#     }
#   }
#   axis(side = 1, at = 1:l_list, labels = xlb, las = 2)
#   if (tilted_text) {
#     text(x = 1:length(yalist), y = min(unlist(yalist))+yoffset, labels = names(yalist), xpd = TRUE, srt = 45)
#   }
#   if (savefile) { dev.copy2pdf(file = FnP_parser(fname, "pdf"), width = w, height = h, title = paste0(basename(fname), " by ", if (exists("scriptname")) scriptname else "Rscript")) }
#   if (incrBottMarginBy) { par("mar" = .ParMarDefault )}
#   assign("plotnameLastPlot", fname, envir = .GlobalEnv)
#   if (mdlink) { MarkDown_Img_Logger_PDF_and_PNG(fname_wo_ext = fname) }
# }


# create_set_SubDir <-function (..., makeOutDirOrig=T, setDir=T) {
#   NewOutDir = kollapse(OutDir,"/", ..., print = F)
#   any_print("All files will be saved under 'NewOutDir': ", NewOutDir)
#   if (!exists(NewOutDir)) {	dir.create(NewOutDir)	}
#   if (setDir) {	setwd(NewOutDir)}
#   if (makeOutDirOrig) {
#     if (exists("OutDirOrig")) any_print("OutDirOrig was defined as:",OutDirOrig)
#     any_print("OutDirOrig will be:", OutDir)
#     assign("OutDirOrig", OutDir, envir = .GlobalEnv)
#   } #if
#
#   assign("OutDir", NewOutDir, envir = .GlobalEnv)
# }


# wvenn <- function (yalist, imagetype = "png", alpha = .5, fill = 1:length(yalist), ..., w = 7, h = 7, mdlink = F) {
#   if (!require("VennDiagram")) { print("Please install VennDiagram: install.packages('VennDiagram')") }
#   fname = kollapse(substitute(yalist), ".", imagetype, print = F)
#   filename = kollapse(OutDir,"/", fname, print = F)
#   subt = kollapse("Total = ", length(unique(unlist(yalist))), " elements in total.", print = F)
#   venn.diagram(x = yalist, imagetype = imagetype, filename = filename, main = substitute(yalist), ... ,
#                sub = subt, fill = fill, alpha = alpha, sub.cex = .75, main.cex = 2)
#   if (mdlink) {
#     llogit(MarkDown_ImgLink_formatter(fname))
#     if (exists("png4Github") & png4Github == T) { llogit(MarkDown_ImgLink_formatter(paste0("Reports/", fname) ) )	}
#   }
# }


# ttl_field <- function (flname = basename(fname) ) { paste0(flname, " by ", if (exists("scriptname")) scriptname else "Rscript") }

### THIS IS DUPLICATE OF THE ONE IN MD REPORTS
# llwrite_list <- function(yalist) {
#   for (e in 1:l(yalist)) {
#     if (is.null( names(yalist) )) { llprint("#####",names(yalist)[e]) } else { llprint("#####", e)}
#     print(yalist[e]); llogit("`", yalist[e], "`")
#   }
# }


#  ------------------------------------
"NEW: 11 September 2017 (Monday) 15:40 "

# color.vec.pal <- function(vector=Size, set = "Set1", ReturnCategoriesToo=F) {
#   NrCol = l(unique(vector))
#   COLZ = RColorBrewer::brewer.pal(NrCol, name = set)[as.factor.numeric(vector)]
#   # if (l(names(vector))) names(COLZ) = names(vector)
#   names(COLZ) =vector
#   CATEG = unique.wNames(COLZ)
#   if (ReturnCategoriesToo) {COLZ = list("vec" = COLZ, "categ" = CATEG)}
#   COLZ
# }
#
# color.vec.base <- function(vector=Plate, set = c(F, "heat.colors", "terrain.colors", "topo.colors", "rainbow")[1], ReturnCategoriesToo=F) {
#   NrCol = l(unique(vector))
#   COLZ = as.factor.numeric(vector) # if basic
#   if(set == "rainbow") {rainbow(NrCol)[COLZ]} else if
#   (set == "heat.colors") {heat.colors(NrCol)[COLZ]} else if
#   (set == "terrain.colors") {terrain.colors(NrCol)[COLZ]} else if
#   (set == "topo.colors") {topo.colors(NrCol)[COLZ]}
#   names(COLZ) = vector
#   CATEG = unique.wNames(COLZ)
#   if (ReturnCategoriesToo) {COLZ = list("vec" = COLZ, "categ" = CATEG)}
#   COLZ
# }
#


# JUNK

# pdfA4plot_on <- function (pname = date(), ..., w = 8.27, h = 11.69, rows = 4, cols = rows-1, mdlink = FALSE,
#                           title = ttl_field(pname)) { # Print (multiple) plots to an (A4) pdf.
#   try.dev.off()
#   assign("mfrow_default", par("mfrow"), fname, envir = .GlobalEnv)
#   fname = FnP_parser(pname, "pdf")
#   pdf(fname, width=w, height=h, title = title)
#   par(mfrow = c(rows, cols))
#   iprint(" ----  Don't forget to call the pair of this function to finish plotting in the A4 pdf.: pdfA4plot_off ()")
#   if (mdlink) { MarkDown_Img_Logger_PDF_and_PNG(fname_wo_ext = pname) }
# }


# wscatter.fill <- function (df2col = cbind("A"=rnorm(100), "B"=rnorm(100)), ..., color, xlim=range(df2col[, 1]), ylim=range(df2col[, 2]), zlim=range(color), nlevels = 20, pch=21, cex=1,
#                            plotname = substitute(df2col), plot.title = plotname,
#                            plot.axes, key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
#                            axes = TRUE, frame.plot = axes, xlb, ylb,
#                            savefile = T, w = 7, h = w, incrBottMarginBy = 0, mdlink = F ) {
#   x = df2col[, 1]
#   y = df2col[, 2]
#   CNN = colnames(df2col)
#   xlb = if(length(CNN) & missing(xlb)) CNN[1]
#   ylb = if(length(CNN) & missing(ylb)) CNN[2]
#
#   fname = kollapse(plotname, ".barplot")
#   if (incrBottMarginBy) { .ParMarDefault <- par("mar"); 	par(mar=c(par("mar")[1]+incrBottMarginBy, par("mar")[2:4]) ) } 	# Tune the margin
#
#   mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
#   on.exit(par(par.orig))
#   WID <- (3 + mar.orig[2L]) * par("csi") * 2.54
#   layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(WID)))
#   par(las = las)
#   mar <- mar.orig
#   mar[4L] <- mar[2L]
#   mar[2L] <- 1
#   par(mar = mar)
#
#   # choose colors to interpolate
#   levels <- seq(zlim[1], zlim[2], length.out = nlevels)
#   col <- colorRampPalette(c("red", "yellow", "dark green"))(nlevels)
#   colz <- col[cut(color, nlevels)]
#
#   plot.new()
#   plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
#
#   rect(0, levels[-length(levels)], 1, levels[-1L], col=col, border=col)
#   if (missing(key.axes)) { if (axes){axis(4)} }
#   else key.axes
#   box()
#   if (!missing(key.title)) key.title
#   mar <- mar.orig
#   mar[4L] <- 1
#   par(mar = mar)
#
#   # points
#   plot(x, y, main =plot.title, type = "n", xaxt='n', yaxt='n', ..., xlim=xlim, ylim=ylim, bty="n", xlab=xlb, ylab=ylb)
#   points(x, y, bg = colz, xaxt='n', yaxt='n', xlab="", ylab="", bty="n", pch=pch, ...)
#
#   ## options to make mapping more customizable
#   if (missing(plot.axes)) {
#     if (axes) {
#       title(main = "", xlab = "", ylab = "")
#       Axis(x, side = 1)
#       Axis(y, side = 2)
#     }
#   }
#   else plot.axes
#   if (frame.plot) box()
#   if (missing(plot.title)) title(...)
#   else plot.title
#   invisible()
#
#   if (savefile) { dev.copy2pdf(file = FnP_parser(fname, "pdf"), width = w, height = h, title = ttl_field(fname)) }
#   if (incrBottMarginBy) { par("mar" = .ParMarDefault )}
#   assign("plotnameLastPlot", fname, envir = .GlobalEnv)
#   if (mdlink) { MarkDown_Img_Logger_PDF_and_PNG(fname_wo_ext = fname)	}
# }


# wvenn <- function (yalist, imagetype = "png", alpha = .5, fill = 1:length(yalist), subt, ..., w = 7, h = w, mdlink = F, plotname = substitute(yalist)) {
#   if (!require("VennDiagram")) { print("Please install VennDiagram: install.packages('VennDiagram')") }
#   fname = kollapse(plotname, ".", imagetype, print = F)
#   LsLen = length(yalist)
#   if(length(names(yalist)) < LsLen) { names(yalist) =1:LsLen; print("List elements had no names.") }
#   print(names(yalist))
#
#   filename = kollapse(OutDir, "/", fname, print = F)
#   if (missing(subt)) { subt = kollapse("Total = ", length(unique(unlist(yalist))), " elements in total.", print = F)  } #if
#   venn.diagram(x = yalist, imagetype = imagetype, filename = filename, main = plotname, ... ,
#                sub = subt, fill = fill, alpha = alpha, sub.cex = .75, main.cex = 2)
#   if (mdlink) {
#     llogit(MarkDown_ImgLink_formatter(fname))
#     if (exists("png4Github") & png4Github == T) { llogit(MarkDown_ImgLink_formatter(paste0("Reports/", fname) ) )	}
#   }
# }


# filter_LP <- function(numeric_vector, threshold, passequal = F, prepend ="", return_survival_ratio=F, na_rm = T) { # Filter values that fall below the low-pass threshold (X <).
#   survivors <- if (passequal) { numeric_vector <= threshold } else { numeric_vector < threshold }
#   pc = percentage_formatter(sum(survivors, na.rm = na_rm)/length(survivors))
#   conclusion = kollapse(prepend, pc, " or ", sum(survivors, na.rm = na_rm), " of ", length(numeric_vector), " entries in ", substitute (numeric_vector), " fall below a threshold value of: ", iround(threshold))
#   if (file.exists(path_of_report) ) {	llogit (conclusion)	} else { print  ("NOT LOGGED") }
#   if (return_survival_ratio) {return (sum(survivors, na.rm = na_rm)/length(survivors))} else if (!return_survival_ratio) { return (survivors) }
# }

#
# md.LogSettingsFromList <-function (parameterlist=p, maxlen =20) {
#   LZ = unlapply(parameterlist, l) # collapse paramters with multiple entires
#   LNG = names(which(LZ>1))
#   for (i in LNG ) {
#     if (l(parameterlist[[LNG]]) > maxlen) parameterlist[[LNG]] = parameterlist[[LNG]][1:maxlen]
#     parameterlist[[LNG]] = paste(parameterlist[[LNG]], collapse = ", ")
#   } #for
#   DF = t(as.data.frame(parameterlist))
#   colnames(DF) = "Value"
#   MarkDown_Table_writer_DF_RowColNames(DF, title_of_table = "Script Parameters and Settings")
# }

# wLinRegression <- function(DF, coeff = c("pearson", "spearman", "r2")[3], textlocation = "topleft", savefile =T, cexx =1, ...) { # Add linear regression, and descriptors to line to your scatter plot. Provide the same dataframe as you provided to wplot() before you called this function
#   # print(coeff)
#   regression <- lm(DF[, 2] ~ DF[, 1])
#   abline(regression, ...)
#   legendText = NULL
#   if ( "pearson" %in% coeff) {    dispCoeff = iround(cor(DF[, 2], DF[, 1], method = "pearson"))
#   legendText  =  c(legendText, paste0("Pears.: ", dispCoeff))  }
#   if ("spearman" %in% coeff) {    dispCoeff = iround(cor(DF[, 2], DF[, 1], method = "spearman"))
#   legendText = c(legendText, paste0("Spear.: ", dispCoeff))  }
#   if ("r2" %in% coeff) {          r2 = iround(summary(regression)$r.squared)
#   legendText = c(legendText, paste0("R^2: ", r2))  }
#   # print(legendText)
#   if (length(coeff)==1 & "r2" == coeff[1]) {  legend(textlocation, legend = superscript_in_plots(prefix = "R", sup = "2", suffix = paste0(": ", r2)) , bty="n", cex = cexx)
#   } else {                                    legend(textlocation, legend = legendText , bty="n", cex = cexx) }
#   if(savefile){   wplot_save_this(plotname = plotnameLastPlot) }
# }


# wlegend <- function(fill_ = "NULL", poz=4, legend, bty = "n", ..., w_=7, h_=w_, OverwritePrevPDF =T) { # Add a legend, and save the plot immediately
#   stopif(is.null(fill_))
#   fNames = names(fill_)
#   if( !is.null(fNames ) ) legend = fNames
#   check_ =(  is.null(fNames) && missing(legend) )
#   stopif( check_, message = "The color vector (fill_) has no name, and the variable 'legend' is not provided.")
#   stopif( ( length(fill_)  != length(legend)), message = "fill and legend are not equally long.")
#
#   pozz = translate(poz, oldvalues = 1:4, newvalues = c("topleft", "topright", "bottomright", "bottomleft"))
#   legend(x=pozz, legend=legend, fill=fill_, ..., bty=bty)
#   if (OverwritePrevPDF) {   wplot_save_this(plotname = plotnameLastPlot, w= w_, h = h_)  }
# }
#
#
# llwrite_list <- function(yalist, printName="self") {
#   if (printName == "self")  llprint("####", substitute(yalist))  else if (printName == F) { ""} else { llprint("####", printName) }  #  else do not print
#   for (e in 1:l(yalist)) {
#     if (is.null( names(yalist) )) { llprint("#####", names(yalist)[e]) } else { llprint("#####", e)}
#     print(yalist[e]); llogit("`", yalist[e], "`")
#   }
# }



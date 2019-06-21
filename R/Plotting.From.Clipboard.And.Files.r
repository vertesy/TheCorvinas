############################################################
## Plotting.From.Clipboard.And.Files.r
############################################################
# source("~/GitHub/TheCorvinas/R/Plotting.From.Clipboard.And.Files.r")


### Plot from clipboard directly -------------------------------------------------------------------------------------------------
# require(MarkdownReports) # See: https://vertesy.github.io/MarkdownReports/

clplot.scatter <-function(..., header = F, col = 1) { # Draw a scatterplot from a 2-column data pasted from clipboard. Works on OS X only.
  DF = fromClipboard(header = header)
  stopifnot(NCOL(DF)==2)
  wplot(DF, savefile = F, col=col)
}

clhist <-function(..., breakz = 20, col = "gold1", xlb = "-") { # Draw a histogram from data pasted from clipboard. Works on OS X only.
  whist(fromClipboard.as_num_vec(), breaks = breakz, savefile = F)
}

clpie <-function(..., percentage_ = TRUE, both_pc_and_value = F, plotname = "Distribution" ) { #  Draw a pie chart from data pasted from clipboard.  Works on OS X only.
  wpie(fromClipboard.as_num_vec(), percentage = percentage_, both_pc_and_value = both_pc_and_value, savefile = F)
}

clbarplot <-function( ..., col_ = "gold1", sub = F) { #  Draw a barplot from data pasted from clipboard.  Works on OS X only.
  wbarplot(fromClipboard.as_num_vec(), col =col_, savefile = F)
}


## Read and write plotting functions READ -------------------------------------
# rw-funcitons cannot call the w-functions, there would be too much things lost: rw writes where the file comes from, w writes in Outdir & current dir


rwplot <- function(FnP, ..., w=7, h=7) { # read in a table, plot it and save it.
  print ('file without header, and file path between ""')
  variable=read.simple  (FnP);
  fname=gsub (".*/", "", FnP);  trunk = gsub ("\\..*", "", fname)
  plot (variable, ..., main=trunk )
  dev.copy2pdf (file=kollapse(FnP, ".plot.pdf"), width=w, height=h )
}

rwscatterplot <- function(FnP, ..., w=7, h=7) { # read in a table, plot it as scatter plot and save it.
  print ('file without header, and file path between ""')
  variable=read.simple.table  (FnP);
  fname=gsub (".*/", "", FnP);  trunk = gsub ("\\..*", "", fname)
  line = lm(variable[, 2]~variable[, 1])
  plot (variable, ..., main=trunk, sub = lm_equation_formatter (line) )
  abline(line)
  dev.copy2pdf (file=kollapse(FnP, ".scatter.pdf"), width=w, height=h )
}

rwboxplot <- function(FnP, col ="gold1", ..., w=7, h=7) { # read in a table, plot it and save it.
  print ('inputfile is a dataframe, with header')
  variable=read.simple.table (FnP);
  fname=gsub (".*/", "", FnP);  trunk = gsub ("\\..*", "", fname)
  boxplot (variable, ..., main=trunk, col =col, las=2)
  dev.copy2pdf (file=kollapse(FnP, ".boxplot.pdf"), width=w, height=h )
}

rwhist <- function(FnP, col ="gold1", ..., w=7, h=7) { # read in a table, plot it and save it. File should be without header, and file path should be given between ""-s.
  print ('file without header, and file path between ""')
  variable=read.simple (FnP);
  fname=gsub (".*/", "", FnP); trunk = gsub ("\\..*", "", fname)
  if ( length (variable) > 0 ) {
    if ( !is.numeric(variable)) { variable = table (variable) ; wbarplot (variable); print ("FILENAME IS: VARIABLE.BARPLOT.PDF") }
    else { 	hist (variable, ..., plotname=trunk, col =col, las=2)
      dev.copy2pdf (file=kollapse(FnP, ".hist.pdf"), width=w, height=h )
    } # if is.numeric
  } # if non empty
}

rwbarplot <- function(FnP, col ="gold1", ..., w=7, h=7) { # read in a vector, plot it and save it. File should be without header, and file path should be given between ""-s.
  print ('file without header, and file path between ""')
  variable=read.simple (FnP);
  fname=gsub (".*/", "", FnP);  trunk = gsub ("\\..*", "", fname)
  barplot (variable, ..., main=trunk, col =col, las=2)
  dev.copy2pdf (file=kollapse(FnP, ".barplot.pdf"), width=w, height=h )
}


# RNA-seq specific R functions -----------------------------------------------------------------------------------------------------
# source("~/Github_repos/TheCorvinas/R/RNA_seq_specific_functions.r")




# Differential Gene Expression ------------------------------------------------------------------------------------------------------------------------------

# source: http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html
wplot_Volcano <- function (DEseqResults, thr_log2fc_ = thr_log2fc, thr_padj_ =thr_padj, showNames = F, saveit =T, pname_ =F, highlight =T) {
  if (pname_ == FALSE) { pname_ = substitute(DEseqResults) }
  Columns =c("Gene", "log2FoldChange", "padj")
  DE = as.data.frame(DEseqResults)
  if (!is.null(rownames(DE)) & !"Gene" %in% colnames(DE)) {    DE = cbind( "Gene" = rownames(DE), DE)  }
  if (sum(! (Columns %in% colnames(DE)))) { any_print("A dataframe with 3 columns needed:", Columns )}

  # Make a basic volcano plot
  subb = paste0("Red if padj<",thr_padj_,", orange of log2FC>",thr_log2fc_,", green if both.")
  with(DE, plot(log2FoldChange, main=pname_, sub=subb, -log10(padj), pch=20, cex=.5, col = rgb(0,0,0,.25),xlim=range(DE$"log2FoldChange") ))

  if(highlight) { # Add colored points:
    with(subset(DE, padj< thr_padj_ ),      points(log2FoldChange, -log10(padj), pch=20, cex=.5, col="red"))
    with(subset(DE, abs(log2FoldChange)>=thr_log2fc_), points(log2FoldChange, -log10(padj), pch=20, cex=.5, col="orange"))
    with(subset(DE, (padj< thr_padj_) & (abs(log2FoldChange)>=thr_log2fc_)), points(log2FoldChange, -log10(padj), pch=20, cex=.75, col="green"))
  }
  # Label points with the textxy function from the calibrate plot
  if (showNames) {    with(subset(DE, padj<thr_padj_ & abs(log2FoldChange)>thr_log2fc_), calibrate::textxy(log2FoldChange, -log10(padj), labs=Gene, cex=.8))  }
  if (saveit) { wplot_save_this(plotname = paste0(pname_,".volcano") )  }
}

filter_DESeq <- function(DESeq_results, thr_log2fc_ =thr_log2fc, thr_padj_=thr_padj) {
  DE = as.data.frame(DESeq_results)
  llprint("#### ", substitute(DESeq_results))
  index_isSign = DE$padj < thr_padj_
  llprint(sum(index_isSign, na.rm = T), "or", pc_TRUE(index_isSign), "of the results is significant at p=",thr_padj_)
  index_FoldChange = (DE$log2FoldChange < -thr_log2fc_ | DE$log2FoldChange >  thr_log2fc_)
  llprint(sum(index_FoldChange, na.rm = T), "or", pc_TRUE(index_FoldChange), "of the results has a fold change more extreme than (+/-)", 2^thr_log2fc)
  index_Hits = index_isSign & index_FoldChange
  llprint(sum(index_Hits, na.rm = T), "or", pc_TRUE(index_Hits), "of the results meet both criteria.")
  DE_hits = iround(DE[ which(index_Hits), ])
  return(DE_hits)
}

prepare4plotMA <- function(DESeq_results, thr_padj_=thr_padj, thr_log2fc_ =F) { # highlight results using 2 thresholds
  DE = as.data.frame(DESeq_results)[, c("baseMean", "log2FoldChange", "padj")]
  index_isSign = DE$"padj" < thr_padj_
  if (thr_log2fc_ != F) {
    index_FoldChange = (DE$log2FoldChange < -thr_log2fc_ | DE$log2FoldChange >  thr_log2fc_)
    DE$"padj" = (index_isSign & index_FoldChange)
  } else { DE$"padj" = index_isSign }
  return(DE)
}


# Correlation plots ------------------------------------------------------------------------------------------------------------------------------

panel.cor.pearson <- function(x, y, digits=2, prefix="", method = "pearson", cex.cor, ...) {  # A function to display correlation values for pairs() function. Default is pearson correlation, that can be set to  "kendall" or "spearman".
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  excludeNAs = which(is.na(x) | is.na(y))
  r <- abs(cor(x[-excludeNAs], y[-excludeNAs], method = method))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.cor.spearman <- function(x, y, digits=2, prefix="", method = "spearman", cex.cor, ...) {  # A function to display correlation values for pairs() function. Default is pearson correlation, that can be set to  "kendall" or "spearman".
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  excludeNAs = which(is.na(x) | is.na(y))
  r <- abs(cor(x[-excludeNAs], y[-excludeNAs], method = method))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, method = cormethod) {
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

# RaceID -----------------------------------------------------------------------------------------------------

id2name <- function(x) sub("\\_\\_chr\\w+","",x) # From RaceID

name2id <- function(x,id=rownames(sc@expdata)) id[sub("\\_\\_chr\\w+","",id) %in% x] # From RaceID



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
	for (i in 1:DisLen) { PCtiles[i] = sum(distribution < distribution[i]) / DisLen	}
	return(PCtiles)
}



#' BaseFrequencies
#'
#' @param mygene Gene of interest. You need either the ID or Gene symbol
#' @param genome mm10 of hg19 for now
#' @export
#'
#' @examples BaseFrequencies()

BaseFrequencies <- function(mygene="Rn45s", genome="mm10", silent=F){ # Gives you the base distribution of a gene of interest, and how extreme it is compared to all transcripts
  MetaDdir = "~/Github_repos/TheCorvinas/Mapping/Reference_Stats/"

  if ( !exists("BaseFrequencies_")) {
    if (genome=="mm10") {        BaseFrequencies_ = read.simple.tsv(MetaDdir, "mm10/BaseFrequencies.mm10.tsv")  }
    else if (genome=="hg19") {   BaseFrequencies_ = read.simple.tsv(MetaDdir, "hg19/BaseFrequencies.hg19.tsv")  }
    assign("BaseFrequencies_", BaseFrequencies_, envir = .GlobalEnv)
  }
  mygene = grep(mygene, rownames(BaseFrequencies_), value = T)
  stopif(condition = (l(mygene)==0),message =  "Gene not found in BaseFrequencies.mm10.tsv or in BaseFrequencies.hg19.tsv")
  frz = BaseFrequencies_[mygene, 1:4]

  x =NULL
  Bases = c("A","C","G","T" )
  if (!silent) {
    for (L in 1:l(Bases)) {    x[L]=ecdf(BaseFrequencies_[ ,Bases[L]])(frz[L])  }
    print("",quote = F)
    print("Position in the distribution of base frequencies across all genes")
    print (percentage_formatter(x))
  }
  return(frz)
}

# BaseFrequencies()


# Transcriptome / Genome Stats -----------------------------------------------------------------------------------------------------

TrLength <- function(mygene="Rn45s", genome="mm10", silent=T){ # Gives you the transctipt length of a gene of interest, and how extreme it is compared to all transcripts
  MetaDdir = "~/Github_repos/TheCorvinas/Mapping/Reference_Stats/"
  if ( !exists("TrLength_")) {
    if (genome=="mm10") {        TrLength_ = read.simple.tsv.named.vector(MetaDdir, "mm10/TranscriptLength.mm10.tsv")  }
    else if (genome=="hg19") {   TrLength_ = read.simple.tsv(MetaDdir, "hg19/TranscriptLength.hg19.tsv")  }
    assign("TrLength_", TrLength_, envir = .GlobalEnv)
  }
  mygene = grep(mygene, names(TrLength_), value = T)
  stopif(condition = (l(mygene)==0),message =  "Gene not found in TranscriptLength.hmm10.tsv or in TranscriptLength.hg19.tsv")
  Len_MyGene = TrLength_[mygene]

  if (!silent) {
    for (L in 1:l(Bases)) {    x[L]=ecdf(BaseFrequencies_)(Len_MyGene)  }
    print("",quote = F)
    print("Position in the distribution of lengths across all genes")
    print (percentage_formatter(x))
  }
  return(Len_MyGene)
}

#  -----------------------------------------------------------------------------------------------------
#  -----------------------------------------------------------------------------------------------------
#  -----------------------------------------------------------------------------------------------------
#  -----------------------------------------------------------------------------------------------------


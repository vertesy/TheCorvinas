
# RNA-seq specific R functions -----------------------------------------------------------------------------------------------------
# source("~/Github_repos/TheCorvinas/R/RNA_seq_specific_functions.r")




# Differential Gene Expression ------------------------------------------------------------------------------------------------------------------------------

# source: http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html
wplot_Volcano <- function (DEseqResults, thr_log2fc_ = thr_log2fc, thr_padj_ =thr_padj, showNames = F, saveit =T, pname_ =F, highlight =T, ...) {
  if (pname_ == FALSE) { pname_ = substitute(DEseqResults) }
  Columns =c("Gene", "log2FoldChange", "padj")
  DE = as.data.frame(DEseqResults)
  if (!is.null(rownames(DE)) & !"Gene" %in% colnames(DE)) {    DE = cbind( "Gene" = rownames(DE), DE)  }
  if (sum(! (Columns %in% colnames(DE)))) { any_print("A dataframe with 3 columns needed:", Columns )}
  DEseqResults = DEseqResults[!is.na(rowSums(DEseqResults)), ]

  # Make a basic volcano plot
  subb = paste0("Red if padj<",thr_padj_,", orange of log2FC>",thr_log2fc_,", green if both.")
  logFC.4plot = DEseqResults$"log2FoldChange"
  padj.4plot = DEseqResults$"padj"
  with(DE, plot(logFC.4plot, main=pname_, sub=subb, -log10(padj.4plot), pch=20, cex=.5, col = rgb(0,0,0,.25), xlim=range(logFC.4plot) ))

  if(highlight) { # Add colored points:
    with(subset(DE, padj< thr_padj_ ),      points(log2FoldChange, -log10(padj), pch=20, cex=.5, col="red"))
    with(subset(DE, abs(log2FoldChange)>=thr_log2fc_), points(log2FoldChange, -log10(padj), pch=20, cex=.5, col="orange"))
    with(subset(DE, (padj< thr_padj_) & (abs(log2FoldChange)>=thr_log2fc_)), points(log2FoldChange, -log10(padj), pch=20, cex=.75, col="green"))
  }
  # Label points with the textxy function from the calibrate plot
  if (showNames) {    with(subset(DE, padj<thr_padj_ & abs(log2FoldChange)>thr_log2fc_), calibrate::textxy(log2FoldChange, -log10(padj), labs=Gene, cex=.8))  }
  if (saveit) { wplot_save_this(plotname = paste0(pname_,".volcano") )  }
}

filter_DESeq <- function(DESeq_results, thr_log2fc_ =thr_log2fc, thr_padj_=thr_padj, usepAdj=T,foldChange_GeoMean) {
  DE = as.data.frame(DESeq_results)
  llprint("#### ", substitute(DESeq_results))

  condition = F
  index_isSign = if (usepAdj) { DE$"padj" <= thr_padj_  } else {DE$"pval" <= thr_padj_}
  llprint(sum(index_isSign, na.rm = T), "or", pc_TRUE(index_isSign), "of the results is significant at p=",thr_padj_)

  index_FoldChange = (DE$"log2FoldChange" <= -thr_log2fc_ | DE$"log2FoldChange" >=  thr_log2fc_)
  llprint(sum(index_FoldChange, na.rm = T), "or", pc_TRUE(index_FoldChange), "of the results has a fold change more extreme than (+/-)", 2^thr_log2fc_)

  if (!missing(foldChange_GeoMean)) {
    index_foldChange_GeoMean = (DE$"foldChange_GeoMean" <= 1/foldChange_GeoMean | DE$"foldChange_GeoMean" >=  foldChange_GeoMean)
    llprint(sum(index_foldChange_GeoMean, na.rm = T), "or", pc_TRUE(index_foldChange_GeoMean), "of the results has a Geometric Mean fold change more extreme than (+/-)", foldChange_GeoMean)
  } #if

  index_Hits = index_isSign & index_FoldChange
  llprint(sum(index_Hits, na.rm = T), "or", pc_TRUE(index_Hits), "of the results meet both criteria.")
  DE_hits = iround(DE[ which(index_Hits), ])
  return(DE_hits)
}


prepare4plotMA <- function(DESeq_results, thr_padj_=thr_padj, thr_log2fc_ =F) { # highlight results using 2 thresholds
  DE = as.data.frame(DESeq_results)[, c("baseMean", "log2FoldChange", "padj")]
  index_notNA = !is.na(DE$"padj") & !is.na(DE$"log2FoldChange")
  index_isSign = (DE$"padj" <= thr_padj_)
  if (thr_log2fc_ != F) {
    index_FoldChange = (na.omit.strip(DE$log2FoldChange <= -thr_log2fc_ | DE$log2FoldChange >=  thr_log2fc_))
    DE$"padj" = (index_isSign & index_FoldChange & index_notNA)
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

id2name <- function(x) sub("\\_\\_chr\\w+","",x ) # From RaceID

name2id <- function(x,id=rownames(sc@expdata)) {
  found = id[sub("\\_\\_chr\\w+","",id) %in% x]
  not_found = setdiff(x, found)
  any_print("NOT FOUND: ", not_found, "or", inline_vec.char(not_found))
  return(found)
} # From RaceID



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

wellname2index <- function(wellnames_vec, wells =384, ZeroPaddedIndices = T, ZeroPaddedWellNames = F) {
  if (ZeroPaddedWellNames) {  wellnames=  paste0(sort(rep(LETTERS[1:16],24)), stringr::str_pad(1:24, 2, pad = "0"))  }
  else {                      wellnames=  paste0(sort(rep(LETTERS[1:16],24)), 1:24) }

  if(ZeroPaddedIndices) { wellindices = stringr::str_pad(1:wells, width = nchar(wells), pad = "0") }
  else {                  wellindices = 1:wells }
  names(wellindices) = wellnames
  assign(x = "wellindices",value =  wellindices, envir = .GlobalEnv)
  print("wellindices vector is in now the memory.")
  wellindices[wellnames_vec]
}

index2wellname <- function(numeric_vec, wells =384, ZeroPaddedIndices = T, ZeroPaddedWellNames = F) {
  if (ZeroPaddedWellNames) {  wellnames=  paste0(sort(rep(LETTERS[1:16],24)), stringr::str_pad(1:24, 2, pad = "0"))  }
  else {                      wellnames=  paste0(sort(rep(LETTERS[1:16],24)), 1:24) }

  if(ZeroPaddedIndices) { names(wellnames) = stringr::str_pad(1:wells, width = nchar(wells), pad = "0") }
  else {                  names(wellnames) = 1:wells }
  assign(x = "wellnames", value = wellnames, envir = .GlobalEnv)
  print("wellnames variable is in now the memory.")
  wellnames[numeric_vec]
}


find.Gene <- function(PartialSymbol, model=c("human", "mouse", "worm")[2], IgnoreCase =F, ...) {
  SubDir = if (model == "human") { "hg19" } else if (model == "mouse") { "mm10" } else if (model == "worm") { "C_elegans" } else {print ("model has to be either of: human, mouse, worm.")}
  MetaDir = p0("~/Github_repos/TheCorvinas/Biology/Sequencing/", SubDir); stopifnot(dir.exists(MetaDir))
  fnp = p0(MetaDir,"/TranscriptLength.",SubDir,".tsv"); stopifnot(file.exists(fnp))
  TL = read.simple.tsv.named.vector(fnp)
  Hits = grep(pattern = PartialSymbol, x=names(TL), ignore.case = IgnoreCase, value = T, ...)
  iprint("Median transcript length in",SubDir,"is:", median(TL))
  print(cbind("Length[nt]" = TL[Hits]))
  return(Hits)
}




stat.Gene.sequence <- function(GeneSymbols=find.Gene("Ssx"), model=c("human", "mouse", "worm")[2], IgnoreCase =F, ...) {
  SubDir = if (model == "human") { "hg19" } else if (model == "mouse") { "mm10" } else if (model == "worm") { "C_elegans" } else {print ("model has to be either of: human, mouse, worm.")}
  MetaDir = p0("~/Github_repos/TheCorvinas/Biology/Sequencing/", SubDir); stopifnot(dir.exists(MetaDir))
  fnpTL = p0(MetaDir,"/TranscriptLength.",SubDir,".tsv"); stopifnot(file.exists(fnpTL))
  TL = read.simple.tsv.named.vector(fnpTL)
  fnpBD = p0(MetaDir,"/BaseFrequencies.",SubDir,".tsv"); stopifnot(file.exists(fnpBD))
  BD = read.simple.tsv(fnpBD)

  NotFound = setdiff(GeneSymbols, names(TL))
  if (length(NotFound)) { iprint("GeneSymbols not found:", NotFound)  } #if
  Found = intersect(GeneSymbols, names(TL))

  median.TL = median(TL)
  median.BD = colMedians(BD, na.rm = T)[1:4]

  RelativeBaseComposition = round(100*(BD[Found  ,1:4]/median.BD))
  colnames(RelativeBaseComposition)  = p0(colnames(RelativeBaseComposition),"%")

  Statz = cbind(  "Length" = TL[Found],
                  "Median Length (% of)" = p0(percentage_formatter(TL[Found]/ median.TL),"    "),
                  RelativeBaseComposition)
  return(Statz)
}



MergeCS1s <- function( InputNames4Libs,PlotIt = T) {
  if (length(InputNames4Libs) !=4 ) { any_print("Too many files provided:", l(InputNames4Libs)); stop()} #if
  Ls_dfs = list.fromNames(p0("Lib",1:4))
  i=1
  for (i in 1:4 ) {
    Ls_dfs[[i]] = read.simple.tsv(p0(InputDir, InputNames4Libs[i]))
  } #for
  l(Ls_dfs)
  RowNZ = lapply(Ls_dfs, rownames)
  SharedGenes = Reduce(intersect, RowNZ)
  colSums96 = list2df(lapply(Ls_dfs, colSums))

  PcOfGenesInMerged = 100*l(SharedGenes)/ unlapply(Ls_dfs, NROW)
  Ls_dfs.filt = lapply(Ls_dfs, select.rows.and.columns, RowIDs = SharedGenes)  # subset down to the shared genes


  L12 = intermingle.cbind(Ls_dfs.filt$"Lib1", Ls_dfs.filt$"Lib2") # intermingle by columns
  L34 = intermingle.cbind(Ls_dfs.filt$"Lib3", Ls_dfs.filt$"Lib4")

  rowStart384 = seq(1, ncol(L12), by = 24);l(rowStart384)
  rowEnd384 = seq(24, ncol(L12), by = 24);l(rowEnd384)
  MergeIDXs = cbind(rowStart384, rowEnd384)

  # stringr::str_pad()
  MTX =matrix.fromNames(rowname_vec = SharedGenes, colname_vec = 1:384)
  i=3
  for (i in 1:8 ) {
    idx1 = MergeIDXs[i,1]
    idx2 = MergeIDXs[i,2]
    colz384 = ((i-1)*48+1): (i*48)
    MTX[ , colz384] = cbind(L12[ , (idx1:idx2)], L34[ , (idx1:idx2)])
  } #for


  if (PlotIt) {

    pdfA4plot_on.layout("CS1.Lib.Merger.Stats")
    wbarplot(colSums(MTX), main = "mRNA per well",ylab="mRNA", xlab="cells", col = richColors(5)[-1],savefile=F)
    wbarplot(PcOfGenesInMerged, ylab="% of genes", tilted_text = T, savefile = F)
    barplot_label(PcOfGenesInMerged, labels = percentage_formatter(PcOfGenesInMerged/100), bottom = T, OverwritePrevPDF = F)
    pdfA4plot_off()

    Combinations = t(combn(1:4, 2))
    pdfA4plot_on("CS1.Primer.Correlations.in.4.Libs", rows = 3)
    for (i in 1:NROW(Combinations) ) {
      COLZ = Combinations[i,]
      DAT = log10(colSums96[, COLZ]+1)
      plot(DAT, pch="", main="log10 library correlation")
      text(DAT, labels = 1:96)
    } #for
    pairs(log10(colSums96+1), pch=1:25)
    pdfA4plot_off(); try.dev.off()

    mRNA.log10 = log10(colSums96[ ,1:2]+1)
    wscatter.fill(mRNA.log10, color = 1:96, pch=21:25, plotname = "mRNA - 1vs2", savefile = T); try.dev.off()

  } #if Plotit
  return(MTX)
}


#  -----------------------------------------------------------------------------------------------------
#  -----------------------------------------------------------------------------------------------------
#  -----------------------------------------------------------------------------------------------------


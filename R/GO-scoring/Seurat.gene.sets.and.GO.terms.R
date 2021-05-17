######################################################################
# Seurat.gene.sets.and.GO.terms.R
######################################################################
# source('~/GitHub/TheCorvinas/R/GO-scoring/Seurat.gene.sets.and.GO.terms.R')

# require(MarkdownReports)
# source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R')
# Setup ------------------------------------------------------------
# BiocManager::install('grimbough/biomaRt')
library(biomaRt); # package.version('biomaRt')
ensembl = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl"
                     , host = "https://grch37.ensembl.org"
                     # , ensemblRedirect = FALSE
                       ) #uses human ensembl annotations


# - PlotGoTermScores
# - IntersectWithExpressed
# - GetGOTerms
# - AddGOGeneList.manual
# - fix.metad.Colname.rm.trailing.1
# - AddGOScore
# - FeaturePlotSaveGO
# - AddCustomScore
# - FeaturePlotSaveCustomScore
# - PasteUniqueGeneList
# - CalcTranscriptomePercentage
# - CalcTranscriptomePercentageGO
# - ww.convert.GO_term.2.score
# - ww.convert.score.2.GO_term
# - GetAllGOTerms


# ------------------------------------------------------------------------
PlotGoTermScores <- function(obj = combined.obj, only.draw.plot = F # Automate retrieving, processing and plotting GO term based gene scores.
                             , openBrowser = F, plot.each.gene = F, verbose = T
                             , GO = "GO:0061621", desc = "canonical.glycolysis", ...) {
  GO.wDot <- make.names(GO)
  ScoreName <- paste0("Score.", GO.wDot)
  print(ScoreName)
  if (DefaultAssay(obj) != 'RNA') { print("DefaultAssay set to RNA"); DefaultAssay(obj) <-  'RNA'}
  if (only.draw.plot) {
    stopifnot(ScoreName %in% colnames(obj@meta.data))
  } else {
    obj <- GetGOTerms(obj = obj, GO = GO, web.open = openBrowser);
    GO.genes <- obj@misc$GO[[ GO.wDot ]]
    if (verbose) iprint(desc, head(GO.genes))
    obj <- AddGOScore(obj = obj, GO = GO);
  }

  plot <- FeaturePlotSaveGO(obj = obj, GO = ScoreName, name_desc = desc)
  if (plot.each.gene) multiFeaturePlot.A4(obj = obj, list.of.genes = GO.genes, foldername = ppp(GO.wDot, desc, 'UMAPs'), ...)
  if (only.draw.plot) return(plot) else return(obj)
}
# "GO:0061621"  "canonical.glycolysis"
# PlotGoTermScores(GO = "GO:0061621", desc = "canonical.glycolysis")


# ------------------------------------------------------------------------
IntersectWithExpressed <- function(genes, obj=combined.obj, genes.shown = 10) { # Intersect a set of genes with genes in the Seurat object.
  print('IntersectWithExpressed()')
  # print(head(genes, n=15))
  diff = setdiff(genes, rownames(obj))
  iprint(length(diff),"genes (of",length(genes), ") are MISSING from the Seurat object:",head(diff, genes.shown))
  return(intersect(rownames(obj), genes))
}
# GO.0010941.regulation.of.cell.death <- IntersectWithExpressed(GO.0010941.regulation.of.cell.death)

# ------------------------------------------------------------------------
GetGOTerms <- function(obj = combined.obj, GO = 'GO:0034976', web.open = T, genes.shown = 10) { # Get GO terms via Biomart package
  print('GetGOTerms()')
  genes <- getBM(attributes = c('hgnc_symbol'), #  'ensembl_transcript_id', 'go_id'
                 filters = "go_parent_term",  uniqueRows = TRUE,
                 values = GO, mart = ensembl)[,1]

  (GO.wDot <- make.names(GO))
  iprint(length(genes), "Gene symbols downloaded:", head(genes, n = genes.shown))
  genes <- IntersectWithExpressed(obj = obj, genes = genes)

  if (is.null(obj@misc$GO)) obj@misc$GO <- list()
  obj@misc$GO[[ GO.wDot ]] <- genes
  iprint("Genes in", GO, "are saved under obj@misc$GO$", GO.wDot)
  if (web.open) system(paste0("open https://www.ebi.ac.uk/QuickGO/search/", GO))
  return(obj)
}
# combined.obj <- GetGOTerms(obj = combined.obj, GO = 'GO:0034976'); combined.obj@misc$GO$GO.0034976

# ------------------------------------------------------------------------
AddGOGeneList.manual <- function(obj = combined.obj, GO = 'GO:0034976', web.open=F  # Add GO terms via Biomart package.
                                 , genes =  c("A0A140VKG3", "ARX", "CNTN2", "DRD1", "DRD2", "FEZF2", "LHX6")) {
  print(head(genes, n = 15))
  genes <- IntersectWithExpressed(obj = obj, genes = genes)

  if (is.null(obj@misc$GO)) obj@misc$GO <- list()
  obj@misc$GO[[make.names(GO)]] <- genes
  iprint("Genes in", GO, "are saved under obj@misc$GO$", make.names(GO))
  if (web.open) system(paste0("open https://www.ebi.ac.uk/QuickGO/search/", GO))
  return(obj)
}
# combined.obj <- AddGOGeneList.manual(obj = combined.obj, GO = 'GO:1904936'
#       , genes =  c("A0A140VKG3", "ARX", "CNTN2", "DRD1", "DRD2", "FEZF2", "LHX6")); combined.obj@misc$GO$GO.0034976


# clean.duplicate.scorenames ------------------------------------------------------------------------------------
clean.duplicate.scorenames <- function(obj = obj) { # Helper. When AddGOScore(), a '1' is added to the end of the column name. It is hereby removed.
  obj <- combined.obj
  cn <- colnames(obj@meta.data)

  nonGO <- sort(grepv(x = cn, pattern = paste0('^Score.GO.[0-9]{7}.*'), invert = TRUE))
  clean <- grepv(x = cn, pattern = paste0('^Score.GO.[0-9]{7}$'))

  appended <- grepv(x = cn, pattern = paste0('^Score.GO.[0-9]{7}\\.[0-9]$'))
  fixed <- gsub(x = appended, pattern = paste0('\\.[0-9]$'), replacement = "")
  fixed.keep <- which( !(fixed %in% clean))
  uniqueGO <- c(clean, fixed.keep)
  obj@meta.data <- obj@meta.data[ , c(nonGO, uniqueGO)]

  iprint(l(clean), "GO's are clean, ", l(appended), "GO's are suffixed by .1 etc, of which"
         , l(fixed.keep), "GO's had no clean counterpart. All",l(uniqueGO), "scores, are cleaned, fixed and unique now.")
  iprint("Metadata column order re-organized alphabetically, and GO-scores at the end.")
  return(obj)
}
# combined.obj <- clean.duplicate.scorenames(obj = combined.obj)


# fix.metad.Colname.rm.trailing.1 ------------------------------------------------------------------------------------
fix.metad.Colname.rm.trailing.1 <- function(obj = obj, colname=ScoreName) { # Helper. When AddGOScore(), a '1' is added to the end of the column name. It is hereby removed.
  colnames(obj@meta.data) <-
    gsub(x = colnames(obj@meta.data)
         , pattern = paste0(colname,1)
         , replacement = colname
    )
  iprint("Trailing '1' in metadata column name is removed. Column name:", colname)
  return(obj)
}
# obj <- fix.metad.Colname.rm.trailing.1(obj = obj, colname=ScoreName)


# ------------------------------------------------------------------------
AddGOScore <- function(obj = combined.obj, GO = "GO:0034976", FixName = TRUE ) { # Call after GetGOTerms. Calculates Score for gene set. Fixes name.
  print("AddGOScore()")
  GO.wDot <- make.names(GO)
  (genes.GO = list(obj@misc$GO[[GO.wDot]]))
  # print(genes.GO)
  (ScoreName = paste0("Score.", make.names(GO)))
  if (!is.list(genes.GO)) genes.GO <-  list(genes.GO) # idk why this structure is not consistent...
  obj <- AddModuleScore(object = obj, features = genes.GO, name = ScoreName)

  if (FixName) obj <- fix.metad.Colname.rm.trailing.1(obj = obj, colname = ScoreName)
  return(obj)
}
# combined.obj <- AddGOScore(obj = combined.obj, GO = "GO:0034976", FixName = TRUE)
# combined.obj$Score.GO.0034976


# ------------------------------------------------------------------------
# name_desc="esponse to endoplasmic reticulum stress"
FeaturePlotSaveGO <- function(obj = combined.obj, GO.score = "Score.GO.0034976", name_desc=NULL
                              , title_ = paste(GO.score, name_desc)
                              , h=7, PNG = T, ...) { # Plot and save a FeaturePlot, e.g. showing gene set scores.
  print("FeaturePlotSaveGO")
  proper.GO <- paste(sstrsplit(GO.score, pattern = "\\.", n = 3)[2:3], collapse = ":")
  (genes.GO = obj@misc$GO[[make.names(proper.GO)]])

  ggplot.obj <-
    FeaturePlot(obj, features = GO.score, min.cutoff = "q05", max.cutoff = "q95", reduction = 'umap', ...) +
    labs(title = title_, caption = paste("Score calc. from",length(genes.GO), "expr. genes from BioMart.", paste0("https://www.ebi.ac.uk/QuickGO/search/", proper.GO)))
  pname = paste0("FeaturePlot.",(GO.score))
  fname = ww.FnP_parser(kpp(pname,name_desc), if (PNG) "png" else "pdf")
  save_plot(filename = fname, plot = ggplot.obj, base_height = h)
  ggplot.obj
}
# FeaturePlotSaveGO()



# AddCustomScore ------------------------------------------------------------------------
AddCustomScore <- function(obj = combined.obj, genes=ALLEN.FRONTAL.found, FixName = TRUE ) { # Call after GetGOTerms. Calculates Score for gene set. Fixes name.
  ls.genes = list(genes)
  if (!is.list(ls.genes)) ls.genes <- list(ls.genes) # idk why this structure is not consistent...
  (ScoreName = ppp("Score", substitute(genes)) )
  obj <- AddModuleScore(object = obj, features = ls.genes, name = ScoreName)

  if (FixName) obj <- fix.metad.Colname.rm.trailing.1(obj = obj, colname = ScoreName)
  return(obj)
}
# combined.obj <- AddCustomScore(obj = combined.obj, genes=ALLEN.FRONTAL.found, FixName = F); colnames(combined.obj@meta.data)



# FeaturePlotSaveCustomScore ------------------------------------------------------------------------
FeaturePlotSaveCustomScore <- function(obj = combined.obj, genes =ALLEN.FRONTAL.found, name_desc=NULL, h=7, PNG =T, ...) { # Plot and save a FeaturePlot, e.g. showing gene set scores.
  ScoreName <- p0('Score.',substitute(genes))

  ggplot.obj <-
    FeaturePlot(obj, features = ScoreName, min.cutoff = "q05", max.cutoff = "q95", reduction = 'umap', ...) +
    labs(title = paste(ScoreName, name_desc), caption = paste("Score calc. from",length(genes), "expr. genes ."))
  pname = paste0("FeaturePlot.",ScoreName)
  fname = ww.FnP_parser(kpp(pname,name_desc), if (PNG) "png" else "pdf")
  save_plot(filename = fname,  plot = ggplot.obj, base_height = h)
  ggplot.obj
}
# FeaturePlotSaveCustomScore()




# ------------------------------------------------------------------------
PasteUniqueGeneList <- function() {
  dput(sort(unique(clipr::read_clip())))
}


# ------------------------------------------------------------------------
CalcTranscriptomePercentage <- function(obj = combined.obj, genes = genes.GO.0061621.can.glyc) {
  total_expr = Matrix::colSums(GetAssayData(object = obj))
  Matrix::colSums(obj[ genes, ]) / total_expr
}

# genes.GO.0061621.canonical.glycolysis <- c("ENO1", "PKLR", "HK2", "HK3", "FOXK1", "PGAM2", "GCK", "BPGM",
#   "PGK1", "ALDOB", "PGM2L1", "PFKP", "HK1", "PGAM1", "GAPDH", "TPI1",
#   "ENO2", "PFKM", "PKM", "ADPGK", "ALDOA", "ENO3", "ALDOC", "FOXK2",
#   "GPI", "GAPDHS", "PFKL")
# Percentile.GO.0061621 <- CalcTranscriptomePercentage(genes = genes.GO.0061621.canonical.glycolysis)
# vioplot::vioplot(Percentile.GO.0061621)

# ------------------------------------------------------------------------
CalcTranscriptomePercentageGO <- function(obj = combined.obj, GO.score = "GO.0061621") {
  total_expr = Matrix::colSums(GetAssayData(object = obj))
  Matrix::colSums(obj[ obj@misc$GO[[GO.score]], ]) / total_expr
}
# Percentile.GO.0061621 <- CalcTranscriptomePercentageGO(GO.score = "GO.0061621")
# vioplot::vioplot(Percentile.GO.0061621)



# ------------------------------------------------------------------------
# extract a named vector of all terms
ww.convert.GO_term.2.score <- function(GO_term = colnames(cor.GO)) {
  gsub(x = GO_term, pattern = ".*?GO:", replacement = "Score.GO.")
}
# ww.convert.GO_term.2.score("GO:0001666")
# ww.convert.GO_term.2.score("lalaGO:0001666")

# ------------------------------------------------------------------------
# extract a named vector of all terms
ww.convert.score.2.GO_term <- function(ScoreNames = colnames(cor.GO)) {
  gsub(x = ScoreNames, pattern = "Score.GO.", replacement = "GO:")
}

# ------------------------------------------------------------------------
GetAllGOTerms <- function(obj=combined.obj, return.obj = T) {
  x <- obj@meta.data
  GOz <- grepv(pattern = "^Score.GO", x = colnames(x), perl = T)
  GOx <- sort.natural(unique(grepv(pattern = "\\.[0-9]$", x = GOz, perl = T, invert = T)))
  GO.names <- Term(object = ww.convert.score.2.GO_term(GOx))

  if (return.obj) {
    obj@misc$"GO.Lookup" <- GO.names
    iprint("GO IDs present in @meta.data are now saved in misc$GO.Lookup")
    cat(head(GO.names), "...")
    return(obj)
  } else {
    return(GO.names)
  }
}

# ------------------------------------------------------------------------
clUMAP.thresholding <- function(q.meta.col = 'Score.GO.0034976', c.meta.col =  "integrated_snn_res.30"
                                , quantile = .95, absolute.cutoff = NULL
                                , obj = combined.obj, plot.barplot = F
                                , subt ="response to ER stress", plotUMAP =T, ... ) {
  clusters.keep <- calc.cluster.averages(col_name = q.meta.col, split_by = c.meta.col, absolute.thr = absolute.cutoff
                                         , quantile.thr = quantile, plotit = plot.barplot)
  cells.2.granules <- as.character(obj[[c.meta.col]][ ,1])


  cluster.nrs.discard <- clusters.keep[which(!clusters.keep)]
  filtered.out <- gsub(pattern = "^cl\\.", replacement = "", x = names(cluster.nrs.discard))
  is.2.highlight <- cells.2.granules %in% filtered.out
  idx.highlight.these <- which(is.2.highlight)
  cells.highlight <- colnames(obj)[idx.highlight.these]

  ttl <- paste('Cells above', 'quantile', quantile, 'in', q.meta.col)
  if (plotUMAP) {
    clUMAP(obj = obj, cells.highlight = cells.highlight
           , plotname = ppp("UMAP.thresholding_q", q.meta.col, quantile, c.meta.col)
           , title = ttl, sub = c.meta.col, raster = F, label.cex = 5, ...) # , shape.by = c.meta.col
  } else {
    # "Not sure about this"
    # x <- calc.cluster.averages(col_name = q.meta.col, split_by = c.meta.col, quantile.thr = quantile, plotit = plot.barplot, simplify = F)
    is.2.highlight
  }
}


# FilterStressedCells ------------------------------------------------------------------------
FilterStressedCells <- function(obj = combined.obj
                                , res = "integrated_snn_res.30"
                                , quantile.thr = 0.9
                                , GOterms = c('glycolytic process' = 'GO:0006096', 'response to endoplasmic reticulum stress' = 'GO:0034976')
                                , direction = "above"
                                , saveRDS = T, saveRDS.Removed = F
                                , PlotExclusionByEachScore = T
                                , PlotSingleCellExclusion = T
                                , GranuleExclusionScatterPlot = T
) {

  # Check arguments ------------------------------------------------------------
  Meta <- obj@meta.data
  MetaVars <- colnames(Meta)
  if( ! res %in% MetaVars ) { iprint('res',res,'is not found in the object.')}

  ScoreNames <- ww.convert.GO_term.2.score(GOterms)
  if( ! all(ScoreNames %in% MetaVars) ) { iprint('Some of the GO-term scores were not found in the object:', ScoreNames, 'Please call first: PlotGoTermScores(), or the actual GetGOTerms()')}

  cells.2.granules <- combined.obj[[res]][ ,1]

  # Exclude granules ------------------------------------------------------------
  mScoresFiltPass <- mScores <-  matrix.fromNames(fill = NaN, rowname_vec = sort(unique(Meta[,res])), colname_vec = make.names(GOterms))
  for (i in 1:length(GOterms) ) {
    scoreX <- as.character(ScoreNames[i]); print(scoreX)
    scoreNameX <- names(ScoreNames)[i]
    mScoresFiltPass[,i] <- calc.cluster.averages(col_name = scoreX, split_by = res, quantile.thr = quantile.thr
                                                 , histogram = T, subtitle = names(ScoreNames)[i]
                                                 , filter =  direction, obj = obj)

    if (GranuleExclusionScatterPlot) {
      mScores[,i] <- calc.cluster.averages(col_name = scoreX, split_by = res, quantile.thr = quantile.thr
                                           , plot.UMAP.too = F, plotit = F
                                           , filter = F)
    }

    if (PlotExclusionByEachScore) {
      # CellsExcludedByScoreX <- which(cells.2.granules %in% which_names(mScoresFiltPass[,i]))
      clUMAP(ident = res, highlight.clusters = which_names(mScoresFiltPass[,i])
             , label = F, title = p0("Stressed cells removed by ", scoreX), plotname =  p0("Stressed cells removed by ", scoreX), sub = scoreNameX
             , sizes.highlight = .5, raster = F)

    }

  }

  # PlotSingleCellExclusion ------------------------------------------------------------
  if (PlotSingleCellExclusion) {
    mSC_ScoresFiltPass <- matrix.fromNames(fill = NaN, rowname_vec = rownames(Meta), colname_vec = make.names(GOterms))
    for (i in 1:length(GOterms) ) {
      scoreX <- as.character(ScoreNames[i]); print(scoreX)
      scoreNameX <- names(ScoreNames)[i]
      ScoreValuesSC <- Meta[,scoreX]
      mSC_ScoresFiltPass[,i] <- filter_HP(numeric_vector = ScoreValuesSC, threshold = quantile(x = ScoreValuesSC, quantile.thr), breaks = 100)
    }
    cells.remove.SC <- rowSums(mSC_ScoresFiltPass)>0
    table(cells.remove.SC)
    sc.filtered <- which_names(cells.remove.SC)
    sc.kept <- which_names(!cells.remove.SC)

    clUMAP(cells.highlight = sc.filtered, label = F
           , title = "Single-cell filtering is not robust to remove stressed cells"
           , plotname = "Single-cell filtering is not robust to remove stressed cells"
           , suffix = "Stress.Filtering", sizes.highlight = .5, raster = F)

    {
      #"TEST"
      MeanZ.ER.stress <- c(
        'mean.stressed' = mean(Meta[sc.filtered, scoreX]),
        'mean.kept' = mean(Meta[sc.kept, scoreX])
      )
      qbarplot(MeanZ.ER.stress, ylab= scoreNameX, xlab.angle = 45, xlab = '', col = as.logical(0:1))
    }

    # saveRDS.Removed
    if (T) {
      obj.Removed.by.SingleCells <- subset(x = obj, cells = sc.filtered) # remove stressed
      isave.RDS(obj.Removed.by.SingleCells, inOutDir = T, suffix = "Removed")
    }



  }

  # Plot excluded cells ------------------------------------------------------------
  PASS <- (rowSums(mScoresFiltPass) == 0)
  granules.excluded <- which_names(!PASS)
  clUMAP(ident = res, highlight.clusters = granules.excluded, label = F, title = "Stressed cells removed", suffix = "Stress.Filtering", sizes.highlight = .5, raster = F)

  if (GranuleExclusionScatterPlot) {
    Av.GO.Scores <- tibble("Average ER-stress score" = mScores[,2]
                           ,"Average Glycolysis score" = mScores[,1]
                           , 'Stressed' = !PASS
                           , 'name' = rownames(mScores)
    )

    if (F) colnames(Av.GO.Scores)[1:2] <- names(GOterms)

    qscatter(Av.GO.Scores, cols = 'Stressed', w = 6, h = 6
             , vline = quantile(mScores[,2], quantile.thr)
             , hline = quantile(mScores[,1], quantile.thr)
             , title = "Groups of Stressed cells have high scores."
             , subtitle = "Thresholded at 90th percentile"
             , label = 'name'
             , repel = T
             , label.rectangle = T
    )

  }


  # Exclude cells & subset ------------------------------------------------------------

  cells.discard <- which(cells.2.granules %in% granules.excluded)
  cells.keep <- which(cells.2.granules %!in% granules.excluded)
  llprint(percentage_formatter(l(cells.keep) / l(cells.2.granules)), "cells kept.")
  # clUMAP(highlight.clusters = filtered.out, ident = rgx, title = "Stressed cells removed", suffix = "Stress.Filtering.cl", sizes.highlight = .5, raster = F
  #        , MaxCategThrHP = 500, label = F, legend = F)

  obj.noStress <- subset(x = obj, cells = cells.keep) # remove stressed
  if (saveRDS) isave.RDS(obj.noStress, inOutDir = T, suffix = "Cleaned")
  if (saveRDS.Removed) {
    obj.RemovedCells <- subset(x = obj, cells = cells.discard) # remove stressed
    isave.RDS(obj.RemovedCells, inOutDir = T, suffix = "Removed")
  }



  return(obj.noStress)
}




# ------------------------------------------------------------------------


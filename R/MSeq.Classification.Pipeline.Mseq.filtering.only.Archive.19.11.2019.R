######################################################################
# MSeq.Pipeline.R
######################################################################
# Based on https://github.com/chris-mcginnis-ucsf/MULTI-seq
# source ('~/GitHub/MULTI-seq.TSC2/MULTI-seq.sample.BCs/MSeq.Pipeline.101147.R')
rm(list=ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
try (source ('~/GitHub/CodeAndRoll/CodeAndRoll.R'),silent= F)
# https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R

require('MarkdownReportsDev') # https://github.com/vertesy/MarkdownReportsDev

# source ('~/Github/TheCorvinas/R/DatabaseLinke.r')
# devtools::install_github('chris-mcginnis-ucsf/MULTI-seq')
library(deMULTIplex)
require(ggplot2);try(ggplot2::theme_set( theme_bw()), silent = TRUE)
require(magrittr)
require(dplyr)
require(cowplot)
# theme_set(theme_cowplot())

BarTableSweepList <- function(min=0.01, max=0.99, step=0.02, bar_table =bar.table) {
  bar.table_sweep.list <- list()
  n <- 0
  Quantiles = seq(from = min, to = max, by=step)
  for (n in 1:length(Quantiles)) { # print(q)
    bar.table_sweep.list[[n]] <- classifyCells(bar_table, q=Quantiles[n])
    names(bar.table_sweep.list)[n] <- paste("q=",Quantiles[n], sep="")
  }
  return(bar.table_sweep.list)
}



# Parameters ------------------------
p=NULL
p$'nr.bc' = 4
p$'LP.bc.ratio' = .25
p$'LP.bc.ratio.2x' = .5
p$'HP.bc' = 250
p$'filterby' ='Abel'

subdir = '101142/'
# subdir = '101143/'

if (subdir == '101142/') {
  MSBC = "bc.used.101142.csv"
  CBC = "barcodes.filtered.101146.tsv"
  r1 = "101442_S9_ME_L001_R1_001.fastq.gz"
  r2 = "101442_S9_ME_L001_R2_001.fastq.gz"
} else {
  MSBC = "bc.used.101143.csv"
  CBC = "barcodes.filtered.101147.tsv"
  r1 = "101443_S10_ME_L001_R1_001.fastq.gz"
  r2 = "101443_S10_ME_L001_R1_001.fastq.gz"
}
seq.library = substr(subdir,start = 1,stop = nchar(subdir)-1)


# Setup ------------------------
OutDirOrig = OutDir = kollapse("~/Dropbox/Abel.IMBA/AnalysisD/Oli/MSeq/",  subdir) # flag.nameiftrue(test)),"/"
setup_MarkdownReports(OutDir = OutDir, scriptname = "MSeq.Pipeline.R")
md.LogSettingsFromList(p)

# 101142/101442_S9_ME_L001_R2_001.fastq.gz
# (inDir = p0("/Volumes/abel/Data/raw.fastq/Oli.set2.MultiSeq/HGHF2DRXX_all/HGHF2DRXX/",subdir))
(inDir = p0("~/Dropbox/Abel.IMBA/AnalysisD/Oli/HGHF2DRXX.local.copy/",subdir))
# inDirMRT = "~/Dropbox/merrits_multiseq/data/"
# Metadata ------------------------

# Read In ------------------------
ReadIn = T
if (ReadIn) {
  CellIDs = read.csv(file = p0(inDir, CBC)) # Vector of cellIDs which the user wants to align.
  CellIDs = substr(CellIDs[,1], start = 0, stop=16)

  # ------------------------
  bar.ref <- read.csv(p0(inDir,MSBC), header =  F)[,1] #V ector of reference MULTI-seq sample barcode sequences.

  # ------------------------
  readTable <- MULTIseq.preProcess(R1 = p0(inDir,r1),
                                   R2 = p0(inDir,r2 ),
                                   cellIDs = CellIDs ,
                                   cell = c(1,16), umi = c(17,26), tag = c(1,8)
  )
  ssaveRDS(object = readTable, filename = ppp("MSeq.readTable", seq.library, idate()))

  bar.table <- MULTIseq.align(readTable, CellIDs, bar.ref )
  ssaveRDS(object = bar.table, filename = ppp("MSeq.bar.table.orig", seq.library, idate(), ".Rds"))
  if(!exists('bar.tableX')) bar.tableX =bar.table

  str(bar.table)
  view.head2(bar.table)

} else {
  load("~/Documents/Rdata.files/",seq.library,"/2019.10.31_08h.Rdata.gz")
  memory.biggest.objects()
}


# bar.tablez=bar.table

"You might not need this cutoff"
topBCrc = apply(bar.table[,1:4], 1, max)
idx.BC.detected = which(topBCrc >= 10)

Top.BC.Index.per.cell = apply(bar.table, 1,which.max )
names.BCs = colnames(bar.table)[1:p$'nr.bc']
Top.BC = names.BCs[Top.BC.Index.per.cell]



ratio.max.per.secondmax = apply(bar.table[idx.BC.detected,1:4], 1, MaxN) / topBCrc[idx.BC.detected]
ratio.vs.top = cbind(ratio.max.per.secondmax,"log2.topBCrc" = log2(topBCrc[idx.BC.detected]))
wplot(ratio.vs.top, abline = 'h', a = log2(10) )

bar.table.HQ = bar.table[idx.BC.detected, ]
View(bar.table.HQ)

# -------

ratio.vs.top = tibble(ratio.max.per.secondmax, "log2.topBCrc" = log2(topBCrc[idx.BC.detected]))

pass = ratio.max.per.secondmax < p$LP.bc.ratio
sum(pass)

# ------------------------------------------------------------------------

UMI.distr = c(
  'UMIs.in.cells' = sum(bar.table.HQ[ pass, 'nUMI_total' ]) ,
  'UMIs.in.bg' = sum(bar.table.HQ[ !pass, 'nUMI_total' ]   )
)
wpie(UMI.distr)
fraction.in.cells = UMI.distr['UMIs.in.cells'] / sum(UMI.distr)
(expected.reads2 = 800e6 * 0.05 /20000 * fraction.in.cells)

fname = ppp("MSeq.readTable", seq.library, idate(), ".Rds")
ssaveRDS(object = readTable, filename = fname)

# ------------------------------------------------------------------------
maxYpoz = floor(max(ratio.vs.top$log2.topBCrc))

A =
  ggplot(data = ratio.vs.top, aes(x=ratio.max.per.secondmax, y = log2.topBCrc)) +
  # ggplot(ratio.vs.top %>% sample_frac(.1), aes(x=ratio.max.per.secondmax, y = log2.topBCrc)) +
  ggtitle(paste("Cells above", p$HP.bc,"BCs' and below", p$LP.bc.ratio,
                "BC-ratio are selected (", pc_TRUE(pass), ")")) +
  geom_vline(xintercept = p$LP.bc.ratio) +
  geom_hline(yintercept = log2(p$HP.bc)) +
  geom_hline(yintercept = log2(expected.reads2), lty=2, col='grey') +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  # geom_jitter(alpha = 0.5, size =.5,show.legend = FALSE, aes(color = pass),width = 0.05, height = 0.05  ) +
  geom_point(alpha = 0.25, size =.2,show.legend = FALSE, aes(color = pass)) +
  geom_text(label = "Expected UMI count", aes(x = .95, y = log2(expected.reads2*1.1) ) ) +
  geom_text(label = "Unlabeled Cells", aes(x = .95, y = log2(p$HP.bc)*0.9, size =3 ) ) +
  geom_text(label = "Multiplets", aes(x = .95, y = maxYpoz, size =3 ) ) +
  geom_text(label = "Labeled Singlets", aes(x = 0.05, y = maxYpoz, size =3 ) ) +
  # scale_y_log10() +
  theme(legend.position = "none") +
  annotation_logticks(); A

save_plot(filename = "filter.pdf", plot = A, base_height=12, ncol=1, nrow=1) #Figure 2


# -----

plotdat = cbind(bar.table, "Top.BC" =Top.BC, "Max.n.UMI" = rowMax(bar.table))# [names(idx.BC.detected), ]
wpie(table(Top.BC))

ggplot(plotdat, aes(Max.n.UMI, fill = Top.BC)) +
  geom_histogram(alpha = 0.5, bins = 60) +
  scale_x_log10()

ggplot(plotdat, aes(x = Top.BC, y = Max.n.UMI, fill = Top.BC)) +
  geom_violin(draw_quantiles = T) +
  geom_boxplot(width=0.1, outlier.size = 1) +
  # stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") +
  scale_y_log10()

MedianRC.per.label = unlapply(split(rowMax(bar.table), Top.BC), median)
wbarplot(MedianRC.per.label, hline = c(500, 600, 700), lcol = "grey66", filtercol = F, tilted_text = T)

ratio.above.BC.thr = ratio.max.per.secondmax[topBCrc[idx.BC.detected] > 10]
whist(log2(topBCrc+1), vline = log2(p$HP.bc+1), filtercol = T, breaks = 50)
whist(ratio.above.BC.thr, vline = p$LP.bc.ratio, filtercol = -1, breaks = 50)

if (!exists("bar.table.original")) bar.table.original = bar.table

if (p$'filterby' == "Abel") bar.table = bar.table[ pass,] # filter by my standards


# QC ------------------------

# ------------------------
# ------------------------

# Step 2: Visually inspect ------------------------

apply(bar.table[idx.BC.detected,1:4], 2, cor)

## Visualize barcode space
bar.tsne <- barTSNE(bar.table[,1:4])

pdf("bc.all.check.pdf")
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none")
  print(g)
}
dev.off()



## Step 3: Sample Classification------------------------
## Perform Quantile Sweep -----------------------------------------------------------------------------------------------------
# bar.table = bar.tableX
(bar.table <- bar.table.full <- bar.table[,1:p$'nr.bc'])
idim(bar.table)

all.neg.cells = NULL
plotpie=T
r = neg.cells = 1
for (r in 1:99) { iprint("classification round", r)

  ## Quantile Sweep List
  bar.table_sweep.list <- BarTableSweepList()

  ## Identify ideal inter-maxima quantile to set barcode-specific thresholds
  threshold.results1 <- findThresh(call.list=bar.table_sweep.list)

  # Plot IQM in the first round
  if (r==1) {
    IMQ = ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
      geom_vline(xintercept=threshold.results1$extrema, lty=2) +
      ggtitle("Inter-Maxima Quantile") +
      scale_color_manual(values=c("red","black","blue")); IMQ
    save_plot(plot = IMQ, filename = "IMQ.pdf")
  }

  ## Finalize round 1 classifications,
  calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
  if (plotpie) wpie(table(calls), plotname = p0('Calls, round ', r))
  table(calls)

  ## remove negative cells
  neg.cells <- which_names(calls == "Negative"); iprint("Negative:", length(neg.cells))
  poz.cells <- which_names(calls != "Negative"); iprint("Positive:", length(poz.cells))

  ## subset to poz cells only
  bar.table <- bar.table[poz.cells, ]
  all.neg.cells <- c(all.neg.cells, neg.cells); iprint("all.neg.cells: ", length(all.neg.cells))

  if (r > 10) {print ("Classification takes too long"); break}
  if (length(neg.cells) == 1) {print ("There are no - cells"); break}
}
"Loop to identify false positives from the previous rounds proper singlets"
## Repeat until all no negative cells remain (usually 3 rounds)...

#  ------------------------
final.calls <- c(calls, vec.fromNames(all.neg.cells, fill = "Negative"))
wpie(table(final.calls));table(final.calls)

Compare.w.Abel=T
if (Compare.w.Abel) {
  pass2 = final.calls[names(pass)]

  A2 =
    ggplot(data = ratio.vs.top, aes(x=ratio.max.per.secondmax, y = log2.topBCrc)) +
    # ggplot(ratio.vs.top %>% sample_frac(.1), aes(x=ratio.max.per.secondmax, y = log2.topBCrc)) +
    ggtitle(paste("Cells above", p$HP.bc,"BCs' and below", p$LP.bc.ratio,
                  "BC-ratio are selected (", pc_TRUE(pass), ")")) +
    geom_vline(xintercept = p$LP.bc.ratio) +
    geom_hline(yintercept = log2(p$HP.bc)) +
    geom_hline(yintercept = log2(expected.reads2), lty=2, col='grey') +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
    # geom_jitter(alpha = 0.5, size =.5,show.legend = FALSE, aes(color = pass),width = 0.05, height = 0.05  ) +
    geom_point(alpha = 1, size =1.5,show.legend = FALSE, aes(color = pass2)) +
    geom_text(label = "Expected UMI count", aes(x = .95, y = log2(expected.reads2*1.1) ) ) +
    geom_text(label = "Unlabeled Cells", aes(x = .95, y = log2(p$HP.bc)*0.9, size =3 ) ) +
    geom_text(label = "Multiplets", aes(x = .95, y = maxYpoz, size =3 ) ) +
    geom_text(label = "Labeled Singlets", aes(x = 0.05, y = maxYpoz, size =3 ) ) +
    # scale_y_log10() +
    annotation_logticks()

  save_plot(filename = "filter.MSeq.final.wo.rescue", plot = A2, base_height=12, ncol=1, nrow=1); A2

}



## Step 4 (optional): Semi-Supervised Negative Cell Reclassification ------------
p$'reclassify' = TRUE
p$'thr.ClassSt' = 16 # def

if (p$'reclassify') {

  ## Perform semi-supervised negative cell reclassification
  reclass.cells <- findReclassCells(bar.table.full, all.neg.cells)
  reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)

  ## Visualize Results
  ClassStability =
  ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) +
    geom_point() + xlim(c(nrow(reclass.res)-1,1)) +
    ylim(c(0,1.05))  +
    geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
    geom_vline(xintercept = p$'thr.ClassSt', color="blue",lty=2)

  save_plot(plot = ClassStability, filename = "ClassStability.pdf"); ClassStability

  ## Finalize negative cell rescue results
  final.calls.rescued <- final.calls
  rescue.ind <- which(reclass.cells$ClassStability >= p$'thr.ClassSt') ## Note: Value will be dataset-specific
  final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]
  wpie(table(final.calls.rescued))
  table(final.calls.rescued)

}

Compare.w.Abel=T
if (Compare.w.Abel) {
  pass2 = final.calls.rescued[names(pass)]

  A3 =
    ggplot(data = ratio.vs.top, aes(x=ratio.max.per.secondmax, y = log2.topBCrc)) +
    # ggplot(ratio.vs.top %>% sample_frac(.1), aes(x=ratio.max.per.secondmax, y = log2.topBCrc)) +
    ggtitle(paste("Cells above", p$HP.bc,"BCs' and below", p$LP.bc.ratio,
                  "BC-ratio are selected (", pc_TRUE(pass), ")")) +
    geom_vline(xintercept = p$'LP.bc.ratio.2x') +
    geom_hline(yintercept = log2(p$'HP.bc')) +
    geom_hline(yintercept = log2(expected.reads2), lty=2, col='grey') +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
    # geom_jitter(alpha = 0.5, size =.5,show.legend = FALSE, aes(color = pass),width = 0.05, height = 0.05  ) +
    geom_point(alpha = 1, size =1.5,show.legend = FALSE, aes(color = pass2)) +
    geom_text(label = "Expected UMI count", aes(x = .95, y = log2(expected.reads2*1.1) ) ) +
    geom_text(label = "Unlabeled Cells", aes(x = .95, y = log2(p$HP.bc)*0.9, size =3 ) ) +
    geom_text(label = "Multiplets", aes(x = .95, y = maxYpoz, size =3 ) ) +
    geom_text(label = "Labeled Singlets", aes(x = 0.05, y = maxYpoz, size =3 ) ) +
    # scale_y_log10() +
    annotation_logticks()
  save_plot(filename = "filter.MSeq.reclassified.pdf", plot = A2, base_height=12, ncol=1, nrow=1); A3
}

"recalculate r2ndmax coz no cutoff"

HardCutOff= F
if (HardCutOff) {
  topBCrc = apply(bar.tableX[,1:4], 1, max)

  ratio.max.per.secondmax = apply(bar.tableX[1:4], 1, MaxN) / topBCrc
  symdiff(names(final.calls.rescued) , names(ratio.max.per.secondmax))
  ratio.max.per.secondmax = ratio.max.per.secondmax[names(final.calls.rescued)]
  # ? any mismatched name?
  sum(names(final.calls.rescued) != names(ratio.max.per.secondmax))

  range(apply(bar.tableX[,1:4], 1, max))

  idx.LowConf =
  (ratio.max.per.secondmax < p$'LP.bc.ratio.2x') &  # less than 2x enriched
    (!final.calls.rescued %in% "Doublet")           # not called doublet
  sum(idx.LowConf)
  range(idx.LowConf)

  final.calls.rescued[idx.LowConf] <- "Doublet"
  table(final.calls.rescued)
}

FixLowCounts= F
if (FixLowCounts) {
  "Does not work, somehow I get too many under 10, does not correspond to the plot"
  # table(final.calls.rescued[which_names(less.then10)])
  #
  # Doublet Negative
  # 105      120
  less.then10 = topBCrc[names(final.calls.rescued)] < p$HP.bc
  sum(less.then10)
  table(final.calls.rescued[which_names(less.then10)])
  final.calls.rescued[which_names(less.then10)] <- "Negative"
}


write.simple.tsv(t(t(final.calls.rescued)), ManualName = ppp('final.calls.rescued', seq.library, idate(),'tsv'))

isave()
oo()

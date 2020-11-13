######################################################################
# Plot.UMI.BC.freuqency.R
######################################################################
# source('~/GitHub/TheCorvinas/R/Plot.UMI.BC.freuqency.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
try (source('~/GitHub/CodeAndRoll/CodeAndRoll.R'),silent= F)
require('MarkdownReportsDev')
require(magrittr)
require(cowplot)
# source('~/Github/TheCorvinas/R/DatabaseLinke.r')

# Setup ------------------------
OutDir = "~/Dropbox/Abel.IMBA/AnalysisD/UMI.and.Cell.Barcode.Frequency/"
InFile= 'reads.per.UMI.1e6.ssv'

setup_MarkdownReports(OutDir = OutDir, scriptname = "Plot.UMI.BC.freuqency.R")
OutDirOrig = OutDir

# Metadata ------------------------

# Parameters ------------------------


# Read In ------------------------
reads.per.UMI <- read.simple.ssv(p0(OutDir, InFile), wRownames = F, colnames =  F); nrow(reads.per.UMI);head(reads.per.UMI)
colnames(reads.per.UMI) <- c("RC", "UMI")
RC<- as.numeric(trimws(reads.per.UMI[,1]))
reads.per.UMI$RC <- RC

# subset ---
reads.per.UMI <- reads.per.UMI[RC>1,]; head(reads.per.UMI); nrow(reads.per.UMI)
UMIs<- reads.per.UMI[,2]

# QC ------------------------

GC_content <- function(string, len=8, pattern = c("G","C")) { # GC-content of a string (frequency of G and C letters among all letters).
  char.list <- stringr::str_split_fixed(string, pattern = "", n = nchar(string))
  tbl = table (factor(unlist(char.list), levels = c("A", "T", "G", "C")))
  sum(tbl[  pattern ]) / sum(tbl)
}

# Calculate ------------------------

UMI.GC_content <- sapply(UMIs, FUN = GC_content); system("say GC Ready")
UMI.A_content <- sapply(UMIs, FUN = GC_content, pattern = c("A")); system("say A Ready")
UMI.T_content <- sapply(UMIs, FUN = GC_content, pattern = c("T")); system("say T Ready")
UMI.C_content <- sapply(UMIs, FUN = GC_content, pattern = c("G")); system("say C Ready")
UMI.G_content <- sapply(UMIs, FUN = GC_content, pattern = c("C")); system("say G Ready")

# ------------------------
reads.per.UMI.stats <- cbind(
                          reads.per.UMI,
                          "GC" = UMI.GC_content,
                          "A" = UMI.A_content,
                          "T" = UMI.T_content,
                          "C" = UMI.C_content,
                          "G" = UMI.G_content
)

write.simple.tsv(reads.per.UMI.stats); system("say Fully Ready")
head(reads.per.UMI.stats)

head(reads.per.UMI.stats)

# Plot ------------------------
pl.UMI <- NULL
{
  pl.UMI$"GC" <-
    reads.per.UMI.stats[1:1000, ] %>%
    ggplot(aes(x = GC, y = RC)) + geom_jitter(size=2) + theme_bw()  +
    scale_y_continuous(trans = 'log10') + annotation_logticks() +
    xlab("Base content - GC (%)") + ylab("Read count") + ggtitle("UMI frequency vs. composition")

  pl.UMI$"A" <-
    reads.per.UMI.stats[1:1000, ] %>%
    ggplot(aes(x = A, y = RC, col=GC) ) + geom_jitter(size=2) + theme_bw() + theme(legend.position = "none") +
    scale_y_continuous(trans = 'log10') + annotation_logticks() +
    xlab("Base content - A (%)") + ylab("Read count") + ggtitle("UMI frequency vs. composition")

  pl.UMI$"T" <-
    reads.per.UMI.stats[1:1000, ] %>%
    ggplot(aes(x = T, y = RC, col=GC) ) + geom_jitter(size=2) + theme_bw() + theme(legend.position = "none") +
    scale_y_continuous(trans = 'log10') + annotation_logticks() +
    xlab("Base content - T (%)") + ylab("Read count") + ggtitle("UMI frequency vs. composition")

  pl.UMI$"G" <-
    reads.per.UMI.stats[1:1000, ] %>%
    ggplot(aes(x = G, y = RC, col=GC) ) + geom_jitter(size=2) + theme_bw() + theme(legend.position = "none") +
    scale_y_continuous(trans = 'log10') + annotation_logticks() +
    xlab("Base content - G (%)") + ylab("Read count") + ggtitle("UMI frequency vs. composition")

  pl.UMI$"C" <-
    reads.per.UMI.stats[1:1000, ] %>%
    ggplot(aes(x = C, y = RC, col=GC) ) + geom_jitter(size=2) + theme_bw()  +
    scale_y_continuous(trans = 'log10') + annotation_logticks() +
    xlab("Base content - C (%)") + ylab("Read count") + ggtitle("UMI frequency vs. composition")
}

plgr <- plot_grid(plotlist = pl.UMI, nrow = 3, ncol = 2)
save_plot(filename ="UMI frequency vs. composition.png", plot = plgr, base_height = hA4, base_width = wA4)




# ------------------------
# ------------------------
Top25 <- iround(head(reads.per.UMI.stats, n = 25))
md.tableWriter.DF.w.dimnames()

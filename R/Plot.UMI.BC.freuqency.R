######################################################################
# Plot.CBC.BC.freuqency.R
######################################################################
# source('~/GitHub/TheCorvinas/R/Plot.CBC.BC.freuqency.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
try (source('~/GitHub/CodeAndRoll/CodeAndRoll.R'),silent= F)
require('MarkdownReportsDev')
require(magrittr)
require(cowplot)
# source('~/Github/TheCorvinas/R/DatabaseLinke.r')

# Setup ------------------------
OutDir = "~/Dropbox/Abel.IMBA/AnalysisD/UMI.and.Cell.Barcode.Frequency/"
# InFile= 'reads.per.CBC.1e6.ssv'
InFile= 'reads.per.CBC.1e6.ssv'

setup_MarkdownReports(OutDir = OutDir, scriptname = "Plot.CBC.BC.freuqency.R")
OutDirOrig = OutDir

# Metadata ------------------------

# Parameters ------------------------


# Read In ------------------------
reads.per.CBC <- read.simple.ssv(p0(OutDir, InFile), wRownames = F, colnames =  F); nrow(reads.per.CBC)
colnames(reads.per.CBC) <- c("RC", "CBC");head(reads.per.CBC)
RC<- as.numeric(trimws(reads.per.CBC[,1]))
reads.per.CBC$RC <- RC
whist(log10(RC))
# subset ---
# reads.per.CBC <- reads.per.CBC[RC>1,]; head(reads.per.CBC); nrow(reads.per.CBC)
# CBCs<- reads.per.CBC[,2]

# QC ------------------------

GC_content <- function(string, len=nchar(string), pattern = c("G","C")) { # GC-content of a string (frequency of G and C letters among all letters).
  char.list <- stringr::str_split_fixed(string, pattern = "", n = nchar(string))
  tbl = table (factor(unlist(char.list), levels = c("A", "T", "G", "C")))
  sum(tbl[  pattern ]) / sum(tbl)
}

# Calculate ------------------------

CBC.GC_content <- sapply(CBCs, FUN = GC_content); system("say GC Ready")
CBC.A_content <- sapply(CBCs, FUN = GC_content, pattern = c("A")); system("say A Ready")
CBC.T_content <- sapply(CBCs, FUN = GC_content, pattern = c("T")); system("say T Ready")
CBC.C_content <- sapply(CBCs, FUN = GC_content, pattern = c("G")); system("say C Ready")
CBC.G_content <- sapply(CBCs, FUN = GC_content, pattern = c("C")); system("say G Ready")

# ------------------------
reads.per.CBC.stats <- cbind(
                          reads.per.CBC,
                          "GC" = CBC.GC_content,
                          "A" = CBC.A_content,
                          "T" = CBC.T_content,
                          "C" = CBC.C_content,
                          "G" = CBC.G_content
)

write.simple.tsv(reads.per.CBC.stats); system("say Fully Ready")
head(reads.per.CBC.stats)

# Plot ------------------------
pl.CBC <- NULL
{
  pl.CBC$"GC" <-
    reads.per.CBC.stats[1:1000, ] %>%
    ggplot(aes(x = GC, y = RC)) + geom_jitter(size=2) + theme_bw()  +
    scale_y_continuous(trans = 'log10') + annotation_logticks() +
    xlab("Base content - GC (%)") + ylab("Read count") + ggtitle("CBC frequency vs. composition")

  pl.CBC$"A" <-
    reads.per.CBC.stats[1:1000, ] %>%
    ggplot(aes(x = A, y = RC, col=GC) ) + geom_jitter(size=2) + theme_bw() + theme(legend.position = "none") +
    scale_y_continuous(trans = 'log10') + annotation_logticks() +
    xlab("Base content - A (%)") + ylab("Read count") + ggtitle("CBC frequency vs. composition")

  pl.CBC$"T" <-
    reads.per.CBC.stats[1:1000, ] %>%
    ggplot(aes(x = T, y = RC, col=GC) ) + geom_jitter(size=2) + theme_bw() + theme(legend.position = "none") +
    scale_y_continuous(trans = 'log10') + annotation_logticks() +
    xlab("Base content - T (%)") + ylab("Read count") + ggtitle("CBC frequency vs. composition")

  pl.CBC$"G" <-
    reads.per.CBC.stats[1:1000, ] %>%
    ggplot(aes(x = G, y = RC, col=GC) ) + geom_jitter(size=2) + theme_bw() + theme(legend.position = "none") +
    scale_y_continuous(trans = 'log10') + annotation_logticks() +
    xlab("Base content - G (%)") + ylab("Read count") + ggtitle("CBC frequency vs. composition")

  pl.CBC$"C" <-
    reads.per.CBC.stats[1:1000, ] %>%
    ggplot(aes(x = C, y = RC, col=GC) ) + geom_jitter(size=2) + theme_bw()  +
    scale_y_continuous(trans = 'log10') + annotation_logticks() +
    xlab("Base content - C (%)") + ylab("Read count") + ggtitle("CBC frequency vs. composition")
}

plgr <- plot_grid(plotlist = pl.CBC, nrow = 3, ncol = 2)
save_plot(filename ="CBC frequency vs. composition.png", plot = plgr, base_height = hA4, base_width = wA4)




# ------------------------
# ------------------------
Top25 <- iround(head(reads.per.CBC.stats, n = 25))
md.tableWriter.DF.w.dimnames()
system("say Script Ready")
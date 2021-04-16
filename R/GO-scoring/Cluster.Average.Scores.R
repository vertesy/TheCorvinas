######################################################################
# Cluster.Average.Scores.R
######################################################################
# source("~/GitHub/Projects/SEO/GO-scoring/Cluster.Average.Scores.R")
# try(dev.off(), silent = T)


# Parameters ------------------------

# Setup ------------------------

# ------------------------

# ------------------------

# ------------------------

# ------------------------

# ------------------------

# ------------------------





{
  "No good"
  BNIP3.STRING.cl.genes <-  c("BCL2", "BCL2L1", "BECN1", "BNIP3", "GABARAP", "GABARAPL1",
                              "GABARAPL2", "HIF1A", "MAP1LC3B", "RHEB", "PINK1")
  multiFeaturePlot.A4(BNIP3.STRING.cl.genes)

  "Seem good, but low"
  Hypoxic.genes.in.cl9.DEG.50 <- c('BNIP3', 'BNIP3L', 'PDK1', 'PGK1', 'FAM162A', 'HILPDA', 'DDIT4', 'P4HB', 'PKM', 'LDHA', 'CA9', 'EGLN3', 'VEGFA', 'ALKBH5')
  multiFeaturePlot.A4(Hypoxic.genes.in.cl9.DEG.50)

  tic(); combined.obj <- RunALRA(object = combined.obj); toc()

  {
    features.plot <- Hypoxic.genes.in.cl9.DEG.50
    DefaultAssay(combined.obj) <- "RNA"
    plot1 <- FeaturePlot(object = combined.obj, features = features.plot, ncol = 1)

    DefaultAssay(combined.obj) <- "alra"
    plot2 <- FeaturePlot(object = combined.obj, features = features.plot, ncol = 1) # , cols = c("lightgrey", "red")

    plotlist <- CombinePlots(list(plot1, plot2), ncol = 2)
    ggsave(plotlist, filename = "umaps.ALRA.vs.RNA.png", width = 10, height = 40)
  }

}




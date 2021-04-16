######################################################################
# Hypoxia.gene.lists.R
######################################################################
# source("~/GitHub/Projects/SEO/GO-scoring/Hypoxia.gene.lists.R")
# try(dev.off(), silent = T)


# Parameters ------------------------

# Setup ------------------------
# Setup ------------------------
OutDir = "~/Dropbox/Abel.IMBA/AnalysisD/SEO/Stress.Scoring/Hypoxia.gene.lists"
setup_MarkdownReports(OutDir = OutDir, scriptname = "Hypoxia.gene.lists.R")
OutDirOrig = OutDir

# ------------------------

# ------------------------

# ------------------------

# ------------------------

# ------------------------

# ------------------------




{
  # combined.obj.bac <- RunALRA(combined.obj.bac)
  # isave.RDS(combined.obj.bac)
  Hypoxia.Genes.NAR.top50 = c('TMEM45A', 'AKAP12', 'ALDOC', 'JMJD1A', 'DDIT4', 'P4HA1', 'SEC24A', 'ANKRD37', 'RSBN1', 'STC2', 'PGK1', 'GOPC',
                              'SAMD12', 'CRKL', 'EDEM3', 'GAPDH', 'TRIM9', 'GOSR2', 'MIF', 'ASPH', 'WDR33', 'DHX40', 'HIG2', 'KLF10', 'PPP1R3C',
                              'R3HDM1', 'RARA', 'LOC162073', 'PGRMC2', 'ZWILCH', 'TPCN1', 'WSB1', 'SPAG4', 'GYS1', 'RRP9', 'SLC25A28', 'NTRK2',
                              'NARF', 'ASCC1', 'UFM1', 'TXNIP', 'MGAT2', 'VDAC1', 'JMJD2B', 'SEC61G', 'SRP19', 'ERO1L', 'PBEF1', 'JMJD2C', 'LOX')
  Hypoxia.Genes.NAR.top100 <- c("TMEM45A", "AKAP12", "ALDOC", "JMJD1A", "DDIT4", "P4HA1", "SEC24A",
                                "ANKRD37", "RSBN1", "STC2", "PGK1", "GOPC", "SAMD12", "CRKL",
                                "EDEM3", "GAPDH", "TRIM9", "GOSR2", "MIF", "ASPH", "WDR33", "DHX40",
                                "HIG2", "KLF10", "PPP1R3C", "R3HDM1", "RARA", "LOC162073", "PGRMC2",
                                "ZWILCH", "TPCN1", "WSB1", "SPAG4", "GYS1", "RRP9", "SLC25A28",
                                "NTRK2", "NARF", "ASCC1", "UFM1", "TXNIP", "MGAT2", "VDAC1",
                                "JMJD2B", "SEC61G", "SRP19", "ERO1L", "PBEF1", "JMJD2C", "LOX",
                                "SNRPD1", "RASSF4", "P4HA2", "PHLDA1", "BHLHB3", "C19orf7", "PPP1R10",
                                "GOLGA1", "INSIG2", "MRPS12", "PIM1", "SOX6", "PJA2", "RRAGD",
                                "GRK6", "C3orf28", "ATF7IP", "MXI1", "XPNPEP1", "C1orf51", "TMF1",
                                "CLINT1", "ILF3", "ZNF710", "CNOT8", "CCNB1", "FGFR1OP2", "SLC7A6",
                                "UGCGL1", "NDRG1", "HES1", "CA12", "MED14", "PPP1R13L", "C1orf119",
                                "GANAB", "EDN1", "STIP1", "TUBG1", "C12orf11", "CXCR4", "AVPI1",
                                "DDIT3", "TTC7A", "SEPHS1", "BCL2L11", "ANKHD1", "NR3C1", "ALDH4A1",
                                "CDCA3")
  dput(Hypoxia.Genes.NAR.top100)
  Overlap.of.Hypox.annotations <- list(
    'GO:0071456' = combined.obj@misc$GO$GO.0071456,
    'NAR Benita 2009 top 100' = Hypoxia.Genes.NAR.top100
  )
  wvenn(Overlap.of.Hypox.annotations)

  multiFeaturePlot.A4(Hypoxia.Genes.NAR.top50)
}

{
  MarkerGenes.7and.9.lfc1 = c("MT1X", "BNIP3", "MT2A", "GAPDH", "FTH1", "VIM", "FTL", "ENO1",
                              "DDIT4", "FAM162A", "DDIT3", "ANGPTL4", "HILPDA", "PGK1", "ADM",
                              "NUPR1", "NEAT1", "PTN", "SLC3A2", "CD9", "PKM", "EIF1", "LGALS3",
                              "SLC2A1", "CHPF", "VEGFA", "SQSTM1", "CD63", "ERO1A", "HERPUD1",
                              "P4HA1", "CDKN1A", "GBE1", "PGAM1", "TRIB3", "HES1", "TPT1",
                              "ANXA5", "PHGDH", "HSPA9", "PSAT1", "P4HB", "HSPA5", "CRYAB",
                              "SCD", "MT3", "SAT1", "LDHA", "IGFBP2", "IGFBP5", "MT1E", "S100A6",
                              "S100A11", "PDK1", "NRN1", "EYS")
  multiFeaturePlot.A4(MarkerGenes.7and.9.lfc1)
  MarkerGenes.7and.9.lfc0.75 = c("CEBPB", "INSIG2", "TIMP1", "CLU", "TPI1", "ATF3", "RPL41",
                                 "SHMT2", "GPX3", "ZFP36L1", "SLC16A3", "STC2", "RCN1", "GADD45A",
                                 "SERPINH1", "PPP1R3C", "PLOD2", "CEBPD", "SERF2", "ZFAS1", "RPS19",
                                 "SEC61G", "BLVRB", "SLC16A1", "EIF2S2", "HLA-C", "HSPB1", "RPL36",
                                 "C4orf3", "NDUFA4L2", "XBP1", "RPS27L", "ANKRD37", "ATF4", "AARS",
                                 "ANXA2", "SDCBP", "HLA-B", "TRMT112", "PSAP", "B2M", "CANX",
                                 "SELENOK", "HLA-A", "BCAN", "ENO2", "C6orf48", "DNAJB9", "AKAP12",
                                 "BNIP3L", "HOPX", "EIF4EBP1", "SERPINE2", "ARF4", "APOE", "MAP1LC3B",
                                 "S100A10", "RPL21", "PKM", "KCTD16", "WSB1", "PGAM1", "EIF1",
                                 "P4HB", "HERPUD1", "FSTL5")
  multiFeaturePlot.A4(MarkerGenes.7and.9.lfc0.75)
}

Idents(combined.obj) <- 'cl.names.KnownMarkers.0.5'

combined.obj@misc$'AverageExpression' <- AverageExpression(object = combined.obj
                                                         , slot = 'data') # , assays = 'RNA'

combined.obj@misc$'Gini'$'cl.names.KnownMarkers.0.5' <- iround(edgeR::gini(t(combined.obj@misc$AverageExpression$RNA)))


Expression.and.Specificity <- tibble(
  "Expression (q90)" = combined.obj@misc$'expr.q90'[MarkerGenes.7and.9.lfc1],
  "Gini Coefficient" = combined.obj@misc$'Gini'$'cl.names.KnownMarkers.0.5'[MarkerGenes.7and.9.lfc1],
  "Genes" = MarkerGenes.7and.9.lfc1
)

qscatter(tbl_X_Y_Col_etc = Expression.and.Specificity, cols = 1,  label = "Genes", repel = TRUE
         , label.select = list(criteria = "`Expression (q90)` >= 0 & `Gini Coefficient` > 0.3"))
#q90 (or 90th quantile) expression is often 0 for LE genes

is(combined.obj@misc$'AverageExpression'$RNA)
g <- "S100A11"
gene.expr.X <- sort(unlist(combined.obj@misc$'AverageExpression'$RNA[g, ]), decreasing = T)
qbarplot(vec = gene.expr.X, title = paste("Average expression", g))

g <- "GAPDH"
gene.expr.X <- sort(unlist(combined.obj@misc$'AverageExpression'$RNA[g, ]), decreasing = T)
qbarplot(vec = gene.expr.X, title = paste("Average expression", g))


x <- GetAssay(combined.obj, assay = "RNA", slot="data")
geomean

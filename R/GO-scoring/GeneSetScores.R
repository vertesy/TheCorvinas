######################################################################
# GeneSetScores.R
######################################################################
# source("~/GitHub/Projects/SEO/GO-scoring/GeneSetScores.R")
# try(dev.off(), silent = T)

if (F) {
  combined.obj <- RecallReduction(obj = combined.obj, dim = 2)
  qUMAP()
}

# Parameters ------------------------
plot.GO.NonNeural <- FALSE
plot.GO.Weird <- FALSE
plotGenes <- FALSE
Go.Term.Correlation  <- FALSE

# Setup ------------------------
source('~/GitHub/Projects/SEO/GO-scoring/Seurat.gene.sets.and.GO.terms.R')
create_set_OutDir(OutDirOrig, "GeneSetScores")

DefaultAssay(combined.obj) <- "RNA"
iprint("DefaultAssay:", DefaultAssay(combined.obj))


ensembl = useEnsembl("ensembl", dataset="hsapiens_gene_ensembl"
                     , host ="https://grch37.ensembl.org"
                     # , ensemblRedirect = FALSE
) #uses human ensembl annotations


# Basic linages  -----------------------------------------------------------------------------------------------
if (F) {

  create_set_SubDir("GO-Basic.linages")
  # GO:0007492 endoderm development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0007492", desc = "endoderm development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0001706 endoderm formation
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001706", desc = "endoderm formation", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0001707 mesoderm formation
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001707", desc = "mesoderm formation", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0007498 mesoderm development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0007498", desc = "mesoderm development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0001705 ectoderm formation
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001705", desc = "ectoderm formation", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0007398 ectoderm development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0007398", desc = "ectoderm development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0001712 ectodermal cell fate commitment
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001712", desc = "ectodermal cell fate commitment", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0001710 mesodermal cell fate commitment
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001710", desc = "mesodermal cell fate commitment", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0001711 endodermal cell fate commitment
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001711", desc = "endodermal cell fate commitment", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  create_set_OutDir(ParentDir)
}

# Non-neural lineages  -----------------------------------------------------------------------------------------------
if (plot.GO.NonNeural) {
  create_set_SubDir("GO-Non.neural.lineages")
  # GO:0051216 cartilage development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0051216", desc = "cartilage development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0048565 digestive tract development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0048565", desc = "digestive tract development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0060348 bone development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0060348", desc = "bone development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0031012 extracellular matrix
  try(combined.obj <- PlotGoTermScores(GO = "GO:0031012", desc = "extracellular matrix", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0062023 collagen-containing extracellular matrix
  try(combined.obj <- PlotGoTermScores(GO = "GO:0062023", desc = "collagen-containing extracellular matrix", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0005615 extracellular space
  try(combined.obj <- PlotGoTermScores(GO = "GO:0005615", desc = "extracellular space", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0043062 extracellular structure organization
  try(combined.obj <- PlotGoTermScores(GO = "GO:0043062", desc = "extracellular structure organization", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0050840 extracellular matrix binding
  try(combined.obj <- PlotGoTermScores(GO = "GO:0050840", desc = "extracellular matrix binding", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  create_set_OutDir(ParentDir)
}

# Neural lineage  -----------------------------------------------------------------------------------------------

  create_set_SubDir("GO-Neural lineage")
  # GO:0042551 neuron maturation
  try(combined.obj <- PlotGoTermScores(GO = "GO:0042551", desc = "neuron maturation", obj = combined.obj, plot.each.gene = plotGenes), silent = F)

  # GO:0030182 neuron differentiation
  try(combined.obj <- PlotGoTermScores(GO = "GO:0030182", desc = "neuron differentiation", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0043005 neuron projection
  try(combined.obj <- PlotGoTermScores(GO = "GO:0043005", desc = "neuron projection", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0048666 neuron development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0048666", desc = "neuron development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  #  GO:0007399 nervous system development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0007399", desc = "nervous system development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  #   GO:0022008 neurogenesis
  try(combined.obj <- PlotGoTermScores(GO = "GO:0022008", desc = "neurogenesis", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0098794 postsynapse
  try(combined.obj <- PlotGoTermScores(GO = "GO:0098794", desc = "postsynapse", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0001764 neuron migration
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001764", desc = "neuron migration", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0048699 generation of neurons
  try(combined.obj <- PlotGoTermScores(GO = "GO:0048699", desc = "generation of neurons", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0030182 neuron differentiation
  try(combined.obj <- PlotGoTermScores(GO = "GO:0030182", desc = "neuron differentiation", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0008038 neuron recognition
  try(combined.obj <- PlotGoTermScores(GO = "GO:0008038", desc = "neuron recognition", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  create_set_OutDir(ParentDir)




# Glia  -----------------------------------------------------------------------------------------------

  create_set_SubDir("GO-Glia")
  # GO:0014002 astrocyte development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0014002", desc = "astrocyte development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0010001 glial cell differentiation
  try(combined.obj <- PlotGoTermScores(GO = "GO:0010001", desc = "glial cell differentiation", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0042063 gliogenesis
  try(combined.obj <- PlotGoTermScores(GO = "GO:0042063", desc = "gliogenesis", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  create_set_OutDir(ParentDir)




# Stress (non HGA) -----------------------------------------------------------------------------------------------

  create_set_SubDir("GO-Stress.general")
  # GO:0006950 response to stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0006950", desc = "response to stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0033554 cellular response to stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0033554", desc = "cellular response to stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0009628 response to abiotic stimulus
  try(combined.obj <- PlotGoTermScores(GO = "GO:0009628", desc = "response to abiotic stimulus", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0033554 cellular response to stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0033554", desc = "cellular response to stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0006970 response to osmotic stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0006970", desc = "response to osmotic stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0035900 response to isolation stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0035900", desc = "response to isolation stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0035902 response to immobilization stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0035902", desc = "response to immobilization stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:1990911 response to psychosocial stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:1990911", desc = "response to psychosocial stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0009651 response to salt stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0009651", desc = "response to salt stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0034059 response to anoxia
  try(combined.obj <- PlotGoTermScores(GO = "GO:0034059", desc = "response to anoxia", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0034616 response to laminar fluid shear stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0034616", desc = "response to laminar fluid shear stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0034405 response to fluid shear stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0034405", desc = "response to fluid shear stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0034976 response to endoplasmic reticulum stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0034976", desc = "response to endoplasmic reticulum stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0071470 cellular response to osmotic stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0071470", desc = "cellular response to osmotic stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0070417 cellular response to cold
  try(combined.obj <- PlotGoTermScores(GO = "GO:0070417", desc = "cellular response to cold", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:1902882 regulation of response to oxidative stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:1902882", desc = "regulation of response to oxidative stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:1905897 regulation of response to endoplasmic reticulum stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:1905897", desc = "regulation of response to endoplasmic reticulum stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # # GO:0000423 mitophagy # Too FEW genes
  # try(combined.obj <- PlotoTermScores(GO = "GO:0000423", desc = "mitophagy", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0000422 autophagy of mitochondrion
  try(combined.obj <- PlotGoTermScores(GO = "GO:0000422", desc = "autophagy of mitochondrion", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  create_set_OutDir(ParentDir)




# HGA | Hypoxia - Glycolysis - Apoptosis -----------------------------------------------------------------------------------------------

  create_set_SubDir("GO-HGA")
  # "GO:0061621"  "canonical.glycolysis"
  try(combined.obj <- PlotGoTermScores(GO = "GO:0061621", desc = "canonical.glycolysis", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0006096 glycolytic.process
  try(combined.obj <- PlotGoTermScores(GO = "GO:0006096", desc = "glycolytic.process", obj = combined.obj, plot.each.gene = plotGenes, openBrowser = T), silent = T)

  # "GO:0043065"  "positive.regulation.of.apoptotic.process"
  try(combined.obj <- PlotGoTermScores(GO = "GO:0043065", desc = "positive.regulation.of.apoptotic.process", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # "GO:0071456"  "cellular.response.to.hypoxia"
  try(combined.obj <- PlotGoTermScores(GO = "GO:0071456", desc = "cellular.response.to.hypoxia", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # "GO:0010942"  "positive.regulation.of.cell.death"
  try(combined.obj <- PlotGoTermScores(GO = "GO:0010942", desc = "positive.regulation.of.cell.death", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # # GO:0034976 response to endoplasmic reticulum stress
  # try(combined.obj <- PlotGoTermScores(GO = "GO:0034976", desc = "response to endoplasmic reticulum stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  create_set_OutDir(NewOutDir = ParentDir)


# cell death-----------------------------------------------------------------------------------------------

  create_set_SubDir("GO-Cell.death")
  # GO:0008219 cell death
  try(combined.obj <- PlotGoTermScores(GO = "GO:0008219", desc = "cell death", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0010942 positive regulation of cell death
  try(combined.obj <- PlotGoTermScores(GO = "GO:0010942", desc = "positive regulation of cell death", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0010941 regulation of cell death
  try(combined.obj <- PlotGoTermScores(GO = "GO:0010941", desc = "regulation of cell death", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0036473 cell death in response to oxidative stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0036473", desc = "cell death in response to oxidative stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0070997 neuron death
  try(combined.obj <- PlotGoTermScores(GO = "GO:0070997", desc = "neuron death", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0070265 necrotic cell death
  try(combined.obj <- PlotGoTermScores(GO = "GO:0070265", desc = "necrotic cell death", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0060548 negative regulation of cell death
  try(combined.obj <- PlotGoTermScores(GO = "GO:0060548", desc = "negative regulation of cell death", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0012501 programmed cell death
  try(combined.obj <- PlotGoTermScores(GO = "GO:0012501", desc = "programmed cell death", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0006915 apoptotic process
  try(combined.obj <- PlotGoTermScores(GO = "GO:0006915", desc = "apoptotic process", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0006914 autophagy
  try(combined.obj <- PlotGoTermScores(GO = "GO:0006914", desc = "autophagy ", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0090398 cellular senescence
  try(combined.obj <- PlotGoTermScores(GO = "GO:0090398", desc = "cellular senescence", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  create_set_OutDir(ParentDir)



# reactive oxygen species ROS  -----------------------------------------------------------------------------------------------

  create_set_SubDir("GO-reactive.oxygen.species.ROS")
  # GO:0000302 response to reactive oxygen species
  try(combined.obj <- PlotGoTermScores(GO = "GO:0000302", desc = "response to reactive oxygen species", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0034614 cellular response to reactive oxygen species
  try(combined.obj <- PlotGoTermScores(GO = "GO:0034614", desc = "cellular response to reactive oxygen species", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:1903426 regulation of reactive oxygen species biosynthetic process
  try(combined.obj <- PlotGoTermScores(GO = "GO:1903426", desc = "regulation of reactive oxygen species biosynthetic process", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:1901031 regulation of response to reactive oxygen species
  try(combined.obj <- PlotGoTermScores(GO = "GO:1901031", desc = "regulation of response to reactive oxygen species", obj = combined.obj, plot.each.gene = plotGenes), silent = T)


  create_set_OutDir(ParentDir)


# dopamine -----------------------------------------------------------------------------------------------

  create_set_SubDir("GO-dopamine")
  # GO:0042417 dopamine metabolic process
  try(combined.obj <- PlotGoTermScores(GO = "GO:0042417", desc = "dopamine metabolic process", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:1903350 response to dopamine
  try(combined.obj <- PlotGoTermScores(GO = "GO:1903350", desc = "response to dopamine", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:1903351 cellular response to dopamine
  try(combined.obj <- PlotGoTermScores(GO = "GO:1903351", desc = "cellular response to dopamine", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0004952 dopamine neurotransmitter receptor activity
  try(combined.obj <- PlotGoTermScores(GO = "GO:0004952", desc = "dopamine neurotransmitter receptor activity", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0007212 dopamine receptor signaling pathway
  try(combined.obj <- PlotGoTermScores(GO = "GO:0007212", desc = "dopamine receptor signaling pathway", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0042416 dopamine biosynthetic process
  try(combined.obj <- PlotGoTermScores(GO = "GO:0042416", desc = "dopamine biosynthetic process", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0004500 dopamine beta-monooxygenase activity
  try(combined.obj <- PlotGoTermScores(GO = "GO:0004500", desc = "dopamine beta-monooxygenase activity", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0001963 synaptic transmission, dopaminergic
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001963", desc = "synaptic transmission, dopaminergic", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  create_set_OutDir(ParentDir)



# Retina.and.eye  -----------------------------------------------------------------------------------------------

  create_set_SubDir("GO-Retina.and.eye")
  # GO:0001895 retina homeostasis
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001895", desc = "retina homeostasis", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0003407 neural retina development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0003407", desc = "neural retina development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0010842 retina layer formation
  try(combined.obj <- PlotGoTermScores(GO = "GO:0010842", desc = "retina layer formation", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0001745 compound eye morphogenesis
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001745", desc = "compound eye morphogenesis", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0048749 compound eye development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0048749", desc = "compound eye development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  create_set_OutDir(ParentDir)



# vasculature  -----------------------------------------------------------------------------------------------

  create_set_SubDir("GO-vasculature")
  # GO:0072359 circulatory system development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0072359", desc = "circulatory system development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0003018 vascular process in circulatory system
  try(combined.obj <- PlotGoTermScores(GO = "GO:0003018", desc = "vascular process in circulatory system", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0001944 vasculature development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001944", desc = "vasculature development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0001568 blood vessel development
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001568", desc = "blood vessel development", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  create_set_OutDir(ParentDir)



# Gene plots -----------------------------------------------------------------------------------------------

genes.GO.0061621.canonical.glycolysis <- c("ENO1", "PKLR", "HK2", "HK3", "FOXK1", "PGAM2", "GCK", "BPGM",
  "PGK1", "ALDOB", "PGM2L1", "PFKP", "HK1", "PGAM1", "GAPDH", "TPI1",
  "ENO2", "PFKM", "PKM", "ADPGK", "ALDOA", "ENO3", "ALDOC", "FOXK2",
  "GPI", "GAPDHS", "PFKL")
Percentile.GO.0061621 <- CalcTranscriptomePercentage(genes = genes.GO.0061621.canonical.glycolysis)
vioplot::vioplot(Percentile.GO.0061621)

genes.GO.0071456.cellular.response.to.hypoxia <- c("GNB1", "MTOR", "PINK1", "HP1BP3", "SLC9A1", "OPRD1", "MPL",
                                                   "AK4", "CLCA1", "LMNA", "PTGS2", "MDM4", "KCNK2", "FMN2", "DNMT3A",
                                                   "KCNK3", "SLC8A1", "PRKCE", "EPAS1", "FABP1", "SCN2A", "NFE2L2",
                                                   "CFLAR", "HIGD1A", "MST1", "CASR", "FAM162A", "OPA1", "CPEB2",
                                                   "PPARGC1A", "NDNF", "CCNA2", "MGARP", "TERT", "PRKAA1", "CCNB1",
                                                   "KCNMB1", "STC2", "EDN1", "PPARD", "VEGFA", "SLC29A1", "FOXO3",
                                                   "TWIST1", "AQP1", "TBL2", "GNGT1", "KCND2", "HILPDA", "PTN",
                                                   "HIPK2", "CYBB", "SUV39H1", "ATP7A", "PGK1", "IRAK1", "NKX3-1",
                                                   "STC1", "BNIP3L", "EIF4EBP1", "SFRP1", "MYC", "NDRG1", "AQP3",
                                                   "UBQLN1", "ENDOG", "BAD", "TRPC6", "CARD16", "HYOU1", "B3GAT1",
                                                   "UCN3", "SUV39H2", "SIRT1", "PTEN", "ANKRD1", "BNIP3", "ADAM8",
                                                   "TIGAR", "PHB2", "MDM2", "SIRT4", "RGCC", "ERO1A", "HIF1A", "ZFP36L1",
                                                   "SLC8A3", "AKT1", "RORA", "CPEB1", "STUB1", "VASN", "MT3", "SLC2A4",
                                                   "TP53", "NPEPPS", "BRIP1", "P4HB", "GATA6", "ACAA2", "PMAIP1",
                                                   "BCL2", "ANGPT4", "COX4I2", "E2F1", "SRC", "PTGIS", "BMP7", "PCK1",
                                                   "ICAM1", "SIRT2", "BBC3", "HMOX1", "MIEF1", "S100B")

Percentile.genes.GO.0071456 <- CalcTranscriptomePercentage(genes = genes.GO.0071456.cellular.response.to.hypoxia)



genes.GO.0010942.positive.regulation.of.cell.death <- c(
                      "MTOR", "ZC3H12A", "FAF1", "PTGS2", "KCNK2", "DNMT3A", "HTRA2",
                      "BCL2L11", "ACOX2", "EIF4G1", "APC", "PRR7", "RIPK1", "CDKN1A",
                      "NOD1", "RAMP3", "CD36", "CALCR", "AIFM1", "NKX3-1", "CRH", "RIPK2",
                      "ACER2", "LPAR1", "HBB", "UCP2", "AKR1C3", "DKK1", "CDKN1B",
                      "LRP1", "DUSP6", "PDX1", "TRIM13", "CIDEB", "BMP4", "HBA2", "HBA1",
                      "MT3", "KATNB1", "HP", "HPR", "ALOX12", "SLFN11", "PHB", "AXIN2",
                      "SAP30BP", "PRNP", "BMP7", "PRODH")
Percentile.GO.0010942 <- CalcTranscriptomePercentage(genes = genes.GO.0010942.positive.regulation.of.cell.death)
vioplot::vioplot(Percentile.GO.0010942*100, ylab="% UMIs")


# irrelevant or very few genes ------------------------------------------------------------------------
if (plot.GO.Weird) {
  create_set_Original_OutDir()
  # GO:0097501 stress response to metal ion
  try(combined.obj <- PlotGoTermScores(GO = "GO:0097501", desc = "stress response to metal ion", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0071472 cellular response to salt stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0071472", desc = "cellular response to salt stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:1990451 cellular stress response to acidic pH
  try(combined.obj <- PlotGoTermScores(GO = "GO:1990451", desc = "cellular stress response to acidic pH", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # did not even run ------------------------------------------------------------------------------------------
  # GO:1990170 stress response to cadmium ion
  try(combined.obj <- PlotGoTermScores(GO = "GO:1990170", desc = "stress response to cadmium ion", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:1990169 stress response to copper ion
  try(combined.obj <- PlotGoTermScores(GO = "GO:1990169", desc = "stress response to copper ion", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # # GO:1990359 stress response to zinc ion
  # try(combined.obj <- PlotGoTermScores(GO = "GO:1990359", desc = "stress response to zinc ion", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # # GO:0097532 stress response to acid chemical
  # try(combined.obj <- PlotGoTermScores(GO = "GO:0097532", desc = "stress response to acid chemical", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # # GO:0010335 response to non-ionic osmotic stress
  # try(combined.obj <- PlotGoTermScores(GO = "GO:0010335", desc = "response to non-ionic osmotic stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0034059 response to anoxia
  try(combined.obj <- PlotGoTermScores(GO = "GO:0034059", desc = "response to anoxia", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0072342 response to anion stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0072342", desc = "response to anion stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:0006979 response to oxidative stress
  try(combined.obj <- PlotGoTermScores(GO = "GO:0006979", desc = "response to oxidative stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  # GO:1903201 regulation of oxidative stress-induced cell death
  try(combined.obj <- PlotGoTermScores(GO = "GO:1903201", desc = "regulation of oxidative stress-induced cell death", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

}


if (F) {
  # GO:0070483 detection of hypoxia
  # GO:0001666 response to hypoxia
  # GO:1900037 regulation of cellular response to hypoxia
  # GO:1900038 negative regulation of cellular response to hypoxia
  # GO:1900039 positive regulation of cellular response to hypoxia
  # GO:1990144 intrinsic apoptotic signaling pathway in response to hypoxia
  # GO:0061418 regulation of transcription from RNA polymerase II promoter in response to hypoxia

  plotGenes =F
  try(combined.obj <- PlotGoTermScores(GO = "GO:0070483", desc = "detection of hypoxia", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  try(combined.obj <- PlotGoTermScores(GO = "GO:0001666", desc = "response to hypoxia", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  try(combined.obj <- PlotGoTermScores(GO = "GO:1900037", desc = "regulation of cellular response to hypoxia", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  try(combined.obj <- PlotGoTermScores(GO = "GO:1900038", desc = "negative regulation of cellular response to hypoxia", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  try(combined.obj <- PlotGoTermScores(GO = "GO:1900039", desc = "positive regulation of cellular response to hypoxia", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  try(combined.obj <- PlotGoTermScores(GO = "GO:1990144", desc = "intrinsic apoptotic signaling pathway in response to hypoxia", obj = combined.obj, plot.each.gene = plotGenes), silent = T)
  try(combined.obj <- PlotGoTermScores(GO = "GO:0061418", desc = "regulation of transcription from RNA polymerase II promoter in response to hypoxia", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

  try(combined.obj <- PlotGoTermScores(GO = "GO:0006950", desc = "response to stress", obj = combined.obj, plot.each.gene = plotGenes), silent = T)

}


# Go.Term.Correlation ------------------------------------------------------------------------------
if (Go.Term.Correlation) source('~/GitHub/Projects/SEO/standalone/Go.Term.Correlation.Analysis.R')


system("say GeneSetScores script ready")
create_set_Original_OutDir()


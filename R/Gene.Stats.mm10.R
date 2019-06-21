
# source("~/GitHub/TheCorvinas/R/Gene.Stats.mm10.R")

"mm10"
BaseFrequencies.mm10 = read.simple.tsv("~/GitHub/TheCorvinas/Biology/Sequencing/mm10/BaseFrequencies.mm10.tsv.gz")
transcripts.length.CDS.mm10 = read.simple.tsv.named.vector("~/GitHub/TheCorvinas/Biology/Sequencing/mm10/TranscriptLength.mm10.tsv.gz")
# Essentiality.OGEE.mm10 = read.simple.tsv("~/GitHub/TheCorvinas/Biology/Sequencing/mm10/Essentiality.OGEE/Mus_musculus_consolidated.tsv.gz")


#  USAGE
# Escapees.or.Noisy = c('Slitrk2__chrX', 'Armcx4__chrX', 'Jpx__chrX', 'Ptchd1__chrX', 'Cask__chrX')
# pdfA4plot_on(pname = "Seq.Check.Escapees.or.Noisy", cols = 2)
# for (i in 1:length(Escapees.or.Noisy[-1]) ) {
#   g=Escapees.or.Noisy[i]
#   check.base.distr(ID = g)
#   check.transcripts.length(ID = g)
# } #for
# pdfA4plot_off()
iprint("BETTER USE WITH EXPRESSED GENES IN YOUR DATASET")


# FUNCTIONS ------------------------------------------------------------------------------------------

cdf.where <- function(value=3, distribution=1:100, verbose=F) {
  CDF = sum(value >= distribution , na.rm = T) / l(distribution)
  if (verbose) { iprint("Your value",value,"sits at", percentage_formatter(CDF), " of the distribution") } #if
  return(CDF)
}
# cdf.where()


check.GC.content <- function(ID="Gnas__chr2", GC_cont=rowSums(BaseFrequencies.mm10[,c("G","C")])) {
  GC.Content.Genes = GC_cont[ID]
  MN = p0("GC-content of ",ID[1], " within all transcripts")
  # SB = kollapse(ID, ": ", LEN.id, " bp. at ", percentage_formatter(CDF), " of the distr." )
  whist(GC_cont, vline = GC.Content.Genes, breaks = 75, main = MN, ylab="Transcripts", xlb="%GC")
  return(GC.Content.Genes)
}






check.transcripts.length <- function(ID="Gnas__chr2", lengths=transcripts.length.CDS.mm10) {
  tr.names = names(lengths)
  if (! ID %in% tr.names) { iprint(ID, "not found.")  } else {
    LEN.id = lengths[ID]
    MN = p0(ID, " and all transcripts' length")
    CDF = cdf.where(LEN.id, lengths)
    SB = kollapse(ID, ": ", LEN.id, " bp. at ", percentage_formatter(CDF), " of the distr." )
    iprint(ID, LEN.id,percentage_formatter(CDF))
    whist(log10(lengths), vline = log10(LEN.id), breaks = 100, main = MN, sub=SB, ylab="Transcripts", xlb="Length [log10(bp)]")
  }
}
# check.transcripts.length()


check.base.distr <- function(ID="Gnas__chr2", Frequencies=100*BaseFrequencies.mm10) {
  tr.names = rownames(Frequencies)
  bases = colnames(Frequencies)
  if (! ID %in% tr.names) { iprint(ID, "not found.")  } else {
    Freqz.ID = Frequencies[ID,]
    MN = p0(ID, " and all base frequencies")
    CDF.freqz =Freqz.ID
    for (i in 1:5 ) {
      CDF.freqz[i]=cdf.where(value = Freqz.ID[,i], distribution = LS.Frequencies[[i]], verbose = F)
    } #for
    LL = paste(names(Freqz.ID),":",Freqz.ID,"% |", percentage_formatter(CDF.freqz))
    EXTREME = c(bases[which(CDF.freqz[1:4] >.9)],bases[which(CDF.freqz[1:4] <.1)])
    if (l(EXTREME)) { SB = iprint("Extreme values: ", paste(EXTREME, percentage_formatter(CDF.freqz[EXTREME]),sep = ":")) } else { SB= ""}
    LS.Frequencies = colsplit(Frequencies)
    wvioplot_list(LS.Frequencies, plotname =  MN, sub = SB, ylb = "Base Frequency [%]")
    FQ = t(rbind(Freqz.ID,1:5))[,2:1]
    points(FQ, lwd=2, bg="gold", pch=21, cex=1.5)
    wlegend(fill_ = 2:5, legend = LL, poz = 2, cex=.5, title = "Base distr. | CDF")
  }
}
# check.base.distr()


compare.length.MWW <- function(IDZ= c('Slitrk2__chrX', 'Armcx4__chrX', 'Jpx__chrX'), lengths=transcripts.length.CDS.mm10) {
  iprint("Mann-Whitney test of genes versus the background distribution.")
  tr.names = names(lengths)
  IDX = !IDZ %in% tr.names
  if ( any(IDX)) { iprint("Not all gene symbols are found: ", IDZ[IDX]) }
  LEN.id = lengths[IDZ]
  # https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
  MWW = suppressWarnings(wilcox.test( LEN.id, lengths))
  SB=p0('p-value @ MWW: ', iround(as.numeric(MWW[3]),2) )
  iprint (SB)
  MN = p0(l(IDZ), " selected transcripts within all transcripts' length")
  CDF = cdf.where(LEN.id, lengths)
  iprint(IDZ, LEN.id,percentage_formatter(CDF))
  whist(log10(lengths), vline = log10(LEN.id), breaks = 100, main = MN, sub=SB, ylab="Transcripts", xlb="Length [log10(bp)]", col="dodgerblue", lwd=2, lcol = 2)
  LL =sort(LEN.id)
  wlegend(fill_ = 2, legend = paste(LL,id2name(names(LL)), sep =": "), poz = 1, title = "Length [bp]", cex=.75)
  wlegend(fill_ = 'dodgerblue', legend = paste(median(lengths),"median", sep =": "), poz = 2, cex=.75)
  return(sapply(LEN.id, cdf.where, distribution = lengths))
}



compare.base.distr.MWW <- function(IDZ= c('Slitrk2__chrX', 'Armcx4__chrX', 'Jpx__chrX'), Frequencies=100*BaseFrequencies.mm10) {
  iprint("Mann-Whitney test of genes versus the background distribution.")
  # tr.names = rownames(Frequencies)
  bases = colnames(Frequencies)
  if ( any(IDX)) { iprint("Not all gene symbols are found: ", IDZ[IDX]) }
  LS.Frequencies = colsplit(Frequencies)
  Freqz.IDZ = Frequencies[IDZ,]
  LS.Freqz.IDZ = colsplit(Freqz.IDZ)

  MN =  p0(l(IDZ), " selected transcripts within all bases' distribution")
  MWW.results = MWW.pval = list.fromNames(bases)
  i=1
  for (i in 1:5 ) {
    MWW.results[[i]] = suppressWarnings(wilcox.test( LS.Freqz.IDZ[[i]], LS.Frequencies[[i]]))
    MWW.pval[[i]] = iround(as.numeric(MWW.results[[i]][3]),2)
    SBx=p0('p-value for ',bases[i],' @ MWW: ', MWW.pval[[i]])
    iprint (SBx)
  } #for
  MWW.pval= unlist(MWW.pval)
  EXTREME = which_names(MWW.pval[1:4] < .1)
  if (l(EXTREME)) { SB = iprint("Significant p-values: ", paste(EXTREME, MWW.pval[EXTREME], sep = ":"), " (at least marg.sign.) ") } else { SB= ""}


  wvioplot_list(LS.Frequencies, plotname =  MN, sub = SB, ylb = "Base Frequency [%]")
  FQ = cbind(
    jitter(sort(rep(1:5, l(IDZ))), amount = .15),
    unlist(Freqz.IDZ)
  )
  points(FQ, lwd=2, bg="gold", pch=21, cex=1)
  LL = p0(bases,": ", MWW.pval)
  wlegend(fill_ = 2:5, legend = LL, poz = 2, cex=.75, title = "p-values @ Wilcoxon")

}
# compare.base.distr.MWW()


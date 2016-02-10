autosomez_exp = sort(unlist(ChrSums[-23,]))
autosomez_exp = autosomez_exp [autosomez_exp > LowerLimit_ReadCount]

borderz = seq (0,l(autosomez_exp), by=100); borderz						# Get bin sizes with 100 entries, by sorting the data!
slices = na.omit(unique(autosomez_exp[borderz]))						# Read count values of the 100th, 200th etc entry
slices[l(slices)] = max(autosomez_exp)									# Extend the upper limit to the maximal data point
slices = c(LowerLimit_ReadCount, slices); slices							# Add the lower limit

hitzzz =  non_hitzzz =  list(0)
percentiles_of_chrX_ALL = hitzzz =  non_hitzzz =  numeric(0)
hitzzz =  non_hitzzz = NullDistr = TestChr = percentiles_of_chrX_ = h = nh = LevelOfSiginif = list(NA)
for (s in 1:(l(slices)-1)) {
# 	 	s=2
	bin_diag = c(slices[s],slices[s+1]);	print (bin_diag)
	InDaBin = (ChrSums >= bin_diag[1]) & (ChrSums < bin_diag[2])
	any_print(sum(InDaBin), "data points in da bin.")
	xIDB 	= InDaBin[23,]
	autoIDB = InDaBin[1:22,]

	NullDistr[[s]] = PP_chr_level[1:22,][autoIDB]
	LevelOfSiginif[[s]] = quantile(NullDistr[[s]], c(0.05, 0.95), na.rm = T)      			# Highly Siginif; # 	LevelOfSiginif[[s]] = quantile(NullDistr[[s]], c(.05, .95), na.rm = T)      				# Marginal; 	# 	LevelOfSiginif[[s]] = quantile(NullDistr[[s]], c(.025, .975), na.rm = T);LevelOfSiginif[[s]]
	TestChr[[s]]   = PP_chr_level[23,][xIDB]
	if (l(TestChr[[s]])){
		percentiles_of_chrX_[[s]] = 100*ecdf(NullDistr[[s]])(TestChr[[s]]);	names(percentiles_of_chrX_[[s]]) = names (TestChr[[s]]);percentiles_of_chrX_[[s]] 		# calculate the concrete percentiles belonging to chrX-data point
		SignifSkewedAllelicExp = TestChr[[s]] < LevelOfSiginif[[s]][1] | TestChr[[s]] > LevelOfSiginif[[s]][2]
		h[[s]] =	names(SignifSkewedAllelicExp[SignifSkewedAllelicExp==T])
		nh[[s]] = names(SignifSkewedAllelicExp[SignifSkewedAllelicExp==F])
	} # if positive length
	else { h[[s]] = nh[[s]] = percentiles_of_chrX_[[s]] = NA }
	any_print("Nr of chrX in the plots: ", l(SignifSkewedAllelicExp))
	any_print("SignifSkewedAllelicExp: ", h[[s]]); 	# any_print("NOT SignifSkewedAllelicExp: ", nh[[s]])
} # for

percentiles_of_chrX_ALL =  sort(unlist (percentiles_of_chrX_))
hitzzz = as.vector(na.omit(unlist (h)))
non_hitzzz =  as.vector(na.omit(unlist (nh)))

a = sort(hitzzz[!Control[hitzzz]])
bb = percentiles_of_chrX_ALL[a]

ccc = sort(hitzzz[Control[hitzzz]]);ccc

NrOfControlCellsAbove100Reads = sum(ChrSums[23, names(which(Control_HQ))] > 100); NrOfControlCellsAbove100Reads

# -------------------------------------------------------------------------------------------------------------------------------------------------
# Plot the outliers the ------------------------------------------------------------------------------------------------------------------------------------

thr  =1500
"not_on_plot= sum(c(MatImpMat, MatImpPat, PatImpMat, PatImpPat)>thr);not_on_plot"

neraly = a


fname = "X_bias"
plot(1, type="n", main=fname,
	 sub = paste ("not_on_plot: ",not_on_plot),
	 xlab="Total Maternal Expression", ylab="Total Paternal Expression", xlim=c(0, thr), ylim=c(0, thr))
points (ChrSums_Mat[23,], ChrSums_Pat[23,], pch =c(3,4), col="darkgreen")


text (ChrSums_Mat[23,][neraly], ChrSums_Pat[23,][neraly], label =neraly, col="darkblue",  adj = c(0,0), cex=1)

for (sl in slices){ 							# Add Slices
	abline (a=sl, b=-1, col ="darkgrey")
}
points (ChrSums_Mat[-23,], ChrSums_Pat[-23,], pch =".", cex =2 )				# Add null distribution
x = c(0, LowerLimit_ReadCount, 0); y = c(0, 0, LowerLimit_ReadCount) 			# Draw triangle in unevaluated region
polygon(x, y, density = NULL, angle = 45, border = NA, col = rgb(0, 0, 0, 0.3))
bins = T
if (bins) {
	slices = slices[-l(slices)] # remove last element
	RC_Floors = slices[-l(slices)]
	RC_Ceilings =  slices[-1]
	NrOfSlices = l(slices)-1
	Limitz = sort(c(RC_Floors, RC_Ceilings) )
	xx = as.data.frame(LevelOfSiginif) #
	colnames(xx) = slices; xx

	Xorder = c(1,2,2,1) # to draw the rectangle, you need the corners anti-clockwise
	for (r in 1:NrOfSlices) { print (r);
								q =r+1
							  Xpz = xx[ Xorder, r]
							  Ypz = (1-xx)[ Xorder, r]
							  mplyr = slices[c(r,r,q,q)]; mplyr;
							  x= mplyr * Xpz
							  y= mplyr * Ypz; x+y
							  polygon (x,y, density = NULL, angle = 45, border = NA, col = rgb(1, 1, 0, 0.2) )
	} # for
} # if (bins)

wplot_save_this(plotname = fname, mdlink = T)


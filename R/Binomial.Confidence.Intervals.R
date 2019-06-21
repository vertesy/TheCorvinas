
# Modified // From : https://www.r-bloggers.com/binomial-confidence-intervals/
PlotIt = F

try (source ('~/GitHub/TheCorvinas/R/CodeAndRoll.R'),silent= F)

irequire(binom)
set.seed(0)
nsims <- 10000
maxn <- 1000
n <- seq(2,maxn, by=2)

my.method <- c("exact", "ac", "asymptotic", "wilson", "prop.test", "bayes", "logit", "cloglog", "probit")
my.method <- my.method[sort.list(my.method)]
# coverage <- matrix(NA, nrow=length(n), ncol=length(my.method))
coverage <- matrix.fromNames(rowname_vec = n, colname_vec = my.method)
# ci.lower <- ci.upper <- matrix(NA, ncol=length(my.method), nrow=length(n))
ci.lower <- ci.upper <- matrix.fromNames(colname_vec = my.method, rowname_vec = n)
  
try.dev.off()

if (PlotIt) { pdfA4plot_on("Binomial.Confidence.Intervals", rows = 2, cols=1) } #if

i=1
for(i in 1:length(n)){
  m <- n[i]/2
  y <- rbinom(n = nsims, size = n[i], prob = .5)
  ll <- binom.confint(m,n[i], conf.level=.95, method=my.method)$lower
  ul <- binom.confint(m,n[i], conf.level=.95, method=my.method)$upper
  ci.lower[i,] <- ll
  ci.upper[i,] <- ul
  for(j in 1:length(my.method)){
    sig <- length(y[y/n[i]<=ul[j] & y/n[i]>=ll[j]])
    coverage[i,j] <- sig/nsims
  }
}

if (PlotIt) { 
  plot(n,NULL, xlim=c(1,nrow(coverage)+1), ylim=c(.83,1),
       col=1, pch=16, ylab="Percent", xlab="N",
       main="95% Confidence Intervals for p=.5")
  
  points(replicate(ncol(coverage),n),coverage, col=c(1:9),
         pch=16, cex=.5)
  
  abline(h=seq(.93,.97, by=.01), col="grey")
  abline(h=.95, col="#000000", lwd=2)
  abline(v=seq(2,maxn, by=20), col="grey")
  legend("bottomright", my.method, col=c(1:9), pch=16, title="Interval Types", bg="#FFFFFF")
  plot(n,NULL, xlim=c(1,100), ylim=c(0,1),
       col=1, pch=16, ylab="Percent", xlab="N",
       main="95% Confidence Interval for p=.5")

  for(k in 1:ncol(coverage)){
    lines(n, ci.lower[,k], col=k, lwd=1)
    lines(n, ci.upper[,k], col=k, lwd=1)
  }
  legend("bottomright", my.method, col=c(1:9), ncol=2, lwd=1, title="Interval Types", bg="#FFFFFF")
  pdfA4plot_off()
  
  } #if PlotIt


Metod = "exact"
Binomial.confidence.intervals.exact = cbind(n, ci.lower[,Metod], ci.upper[,Metod])
write.simple.tsv(Binomial.confidence.intervals.exact)

plot(Binomial.confidence.intervals.exact[1:50,1:2], type = "l")


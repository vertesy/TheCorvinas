######################################################################
# Fit.Exponential.Decay
######################################################################
# source ('~/GitHub/TheCorvinas/R/Fit.Exponential.Decay.R')
# rm(list=ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
try (source ('~/GitHub/CodeAndRoll/CodeAndRoll.R'),silent= F)

# Setup ------------------------
# OutDir = ""
# setup_MarkdownReports(OutDir = OutDir, scriptname = "")
# OutDirOrig = OutDir

# Metadata ------------------------

# Parameters ------------------------
set.seed(1223)
len <- 24

# Data ------------------------
x <- runif(len)
y <- 1/x^2 + rnorm(len, 0, 0.16)
ds <- data.frame(x = x, y = y)
str(ds)
plot(y ~ x, main = "Known square root, with noise")
s <- seq(0, 1, length = 100)

# Fit an exponential model ------------------------
m.e <- nls(y ~ I(exp(1)^(a + b * x)), data = ds, start = list(a = 0, b = 1), trace = T)
summary(m.e)$coefficients
summary(m.e)

# Plot ------------------------
a <- round(summary(m.e)$coefficients[1, 1], 4)
b <- round(summary(m.e)$coefficients[2, 1], 4)
plot(y ~ x, main = "Fitted exponential function", sub = "Blue: nls() fit | Red: lm() fit", xlim=c(0, max(x)) )
  abline (h=seq(0,120, by = 10), lty=3, col="grey1")
  abline (v=2:4*t0.5, lty=3)
s <- seq(0, 1, length = 100)

# Draw fitted model and equation ------------------------
lines(s, predict(m.e, list(x = s)), lty = 1, col = "blue")

# Equations
text(0.2, 100, paste("y =e^ (", a, " + ", b, " * x)", sep = ""), pos = 4)
text(0.2, 80, paste("y =", iround(exp(1)^a), "* e^ ", iround(b), "*x", sep = ""), pos = 4)

t0.5 = -iround(0.693/b)
text(0.2, 60, paste("t[1/2]: ",  t0.5, sep = ""), pos = 4)
# wlegend.old("topleft", legend = eq, cex=.7, title="Exponential decay model")

# No R-squared ------------------------
"Getting an R-squared is tricky: you can only calculate it with a linear model, however that looks very different"


# Fit an exponential model ------------------------
exponential.model <- lm( log(y) ~ x, data = ds)
summary(exponential.model)
# plot(m.lm) #stats
# plot(y ~ x, main = "Fitted exponential function", sub = "Blue: fit", xlim=c(0, max(x)) )
Counts.exponential2 <- exp(predict(exponential.model,list(Time=sort(x))))
# add fitted curve
lines(sort(x), Counts.exponential2[order(x)],lwd=2, col = "red", xlab = "Time (s)", ylab = "Counts")

wplot_save_this(plotname = "Fit.Exponential.Decay")


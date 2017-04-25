# source("/Users/abelvertesy/Github_repos/TheCorvinas/R/diffexpnb_function.2017.03.28.Dominic.New.R")

diffexpnb <- function(x, A, B, DESeq=FALSE, method="pooled", norm=FALSE, vfit=NULL, locreg=FALSE, ...){
  if ( ! method %in% c("per-condition", "pooled") ) stop("invalid method: choose pooled or per-condition")
  x <- x[, c(A, B)]
  if ( DESeq ){
    require(DESeq2)
    # run on sc@expdata
    des <- data.frame( row.names = colnames(x), condition = factor(c( rep(1, length(A)), rep(2, length(B)) )), libType = rep("single-end", dim(x)[2]))
    cds <- DESeqDataSetFromMatrix(countData=round(x, 0), colData=des, design =~ condition, fitType='local', ...)
    res <- results(cds)
    list(des=des, cds=cds, res=res)
  }else{
    if (norm) x <- as.data.frame( t(t(x)/apply(x, 2, sum))*min(apply(x, 2, sum, na.rm=TRUE)) )
    fit <- list()
    m   <- list()
    v   <- list()
    for ( i in 1:2 ){
      g <- if ( i == 1 ) A else B
      m[[i]] <- if ( length(g) > 1 ) apply(x[, g], 1, mean) else x[, g]
      v[[i]] <- if ( length(g) > 1 ) apply(x[, g], 1, var)  else apply(x, 1, var)
      if ( method == "pooled"){
        mg <- apply(x, 1, mean)
        vg <- apply(x, 1, var)
        vl <- log2(vg)
        ml <- log2(mg)
      }else{
        vl <- log2(v[[i]])
        ml <- log2(m[[i]])
      }

      if ( locreg ){
        f <- order(ml, decreasing=FALSE)
        u <- 2**ml[f]
        y <- 2**vl[f]
        lf <- locfit(y~lp(u, nn=.7), family="gamma", maxk=500)
        fit[[i]] <- approxfun(u, fitted(lf), method = "const")
      }else{
        if ( is.null(vfit) ){
          f <- ml > -Inf & vl > -Inf
          ml <- ml[f]
          vl <- vl[f]
          mm <- -8
          repeat{
            fit[[i]] <- lm(vl ~ ml + I(ml^2))
            if( coef(fit[[i]])[3] >= 0 | mm >= -1){
              break
            }
            mm <- mm + .5
            f <- ml > mm
            ml <- ml[f]
            vl <- vl[f]
          }
        }else{
          fit[[i]] <- vfit
        }
      }
    }

    if ( locreg ){
      vf  <- function(x, i) fit[[i]](x)
    }else{
      vf  <- function(x, i) 2**(coef(fit[[i]])[1] + log2(x)*coef(fit[[i]])[2] + coef(fit[[i]])[3] * log2(x)**2)
    }
    sf  <- function(x, i) x**2/(max(x + 1e-6, vf(x, i)) - x)

    psp <- 1e-99
    pv <- apply(data.frame(m[[1]], m[[2]]), 1, function(x){ p12 <- (dnbinom(0:round(x[1]*length(A) + x[2]*length(B), 0), mu=mean(x)*length(A), size=length(A)*sf(mean(x), 1)) + psp)*(dnbinom(round(x[1]*length(A) + x[2]*length(B), 0):0, mu=mean(x)*length(B), size=length(B)*sf(mean(x), 2)) + psp); sum(p12[p12 <= p12[round(x[1]*length(A), 0) + 1]])/sum(p12)} )

    res <- data.frame(baseMean=(m[[1]] + m[[2]])/2, baseMeanA=m[[1]], baseMeanB=m[[2]], foldChange=m[[2]]/m[[1]], log2FoldChange=log2(m[[2]]/m[[1]]), pval=pv, padj=p.adjust(pv, method="BH"))
    vf1 <- data.frame(m=m[[1]], v=v[[1]], vm=vf(m[[1]], 1))
    vf2 <- data.frame(m=m[[2]], v=v[[2]], vm=vf(m[[2]], 2))
    rownames(res) <- rownames(x)
    rownames(vf1) <- rownames(x)
    rownames(vf2) <- rownames(x)
    list(vf1=data.frame(m=m[[1]], v=v[[1]], vm=vf(m[[1]], 1)), vf2=data.frame(m=m[[2]], v=v[[2]], vm=vf(m[[2]], 2)), res=res)
  }
}

plotdiffgenesnb <- function(x, pthr=.05, padj=TRUE, lthr=1, mthr=-Inf, Aname=NULL, Bname=NULL, show_names=TRUE, lthr_line=TRUE, ...){
  y <- as.data.frame(x$res)
  if ( is.null(Aname) ) Aname <- "baseMeanA"
  if ( is.null(Bname) ) Bname <- "baseMeanB"

  plot(log2(y$baseMean), y$log2FoldChange, pch=20, ..., col="grey",
       xlab=paste("log2 ( ( #mRNA[", Aname, "] + #mRNA[", Bname, "] )/2 )", sep=""),
       ylab=paste("log2 #mRNA[", Bname, "] - log2 #mRNA[", Aname, "]", sep=""))
  abline(0, 0)
  if(lthr_line){abline(h=c(-lthr,lthr), lty=3, col="grey")}
  if ( ! is.null(pthr) ){
    if ( padj ) f <- y$padj < pthr else f <- y$pval < pthr
    points(log2(y$baseMean)[f], y$log2FoldChange[f], col="red", pch=20)
  }
  if ( !is.null(lthr) ) f <- f & abs( y$log2FoldChange ) > lthr
  if ( !is.null(mthr) ) f <- f & log2(y$baseMean) > mthr
  if ( show_names ){
    if ( sum(f) > 0 ) text(log2(y$baseMean)[f], y$log2FoldChange[f], rownames(y)[f], cex=.5)
  }
}

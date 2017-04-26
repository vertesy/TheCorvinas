# source("~/Github_repos/TheCorvinas/R/diffexpnb_AV.r")

diffexpnb_AV <- function(x, cells_of_interest, cells_background , norm=TRUE, DESeq=FALSE, method="per-condition", vfit=NULL, locreg=FALSE){
  if ( ! method %in% c("per-condition","pooled","pooled-CR") ) stop("invalid method")
  x <- x[,c(cells_background,cells_of_interest)]
  if ( DESeq ){
    des <- data.frame( "row.names" = colnames(x),
                       "condition" = c( rep(1,length(cells_background)), rep(2,length(cells_of_interest)) ),
                       "libType" = rep("single-end", dim(x)[2]))
    cds <- newCountDataSet( round(x,0), des$condition )
    cds <- estimateSizeFactors( cds )
    cds <- estimateDispersions( cds, method=method, fitType="local" )
    res <- nbinomTest( cds, 1, 2 )
    rownames(res) <- res$id
    res <- res[,-1]
    list("des"=des, "cds"=cds, "res"=res)

  } else {
    if (norm) x <- as.data.frame( t(t(x)/apply(x,2,sum))*median(apply(x,2,sum,na.rm=TRUE)) )
    fit = Meanz = Medianz = v = list()

    for ( i in 1:2 ){
      g <- if ( i == 1 ) cells_background else cells_of_interest

      Meanz[[i]] <- if ( length(g) > 1 ) apply(x[,g],1,mean) else x[,g]
      Medianz[[i]] <- if ( length(g) > 1 ) apply(x[,g],1,median) else x[,g]

      v[[i]] <- if ( length(g) > 1 ) apply(x[,g],1,var)  else apply(x,1,var)

      if ( method == "pooled"){
        mg <- apply(x,1,mean)
        vg <- apply(x,1,var)
        f <- vg > 0 & mg > .5
        logv <- log2(vg[f])
        logm <- log2(mg[f])
      } else {
        f <- v[[i]] > 0 & Meanz[[i]] > .5
        logv <- log2(v[[i]][f])
        logm <- log2(Meanz[[i]][f])
      }

      if ( locreg ){
        f <- order(logm,decreasing=FALSE)
        u <- 2**logm[f]
        y <- 2**logv[f]
        lf <- locfit(y~lp(u,nn=.7),family="gamma",maxk=500)
        fit[[i]] <- approxfun(u, fitted(lf), method = "const")
      }else{
        fit[[i]] <- if ( is.null(vfit) ) lm(logv ~ logm + I(logm^2)) else vfit
      }
    }

    if ( locreg ){
      vf  <- function(x,i) fit[[i]](x)
    }else{
      vf  <- function(x,i) 2**(coef(fit[[i]])[1] + log2(x)*coef(fit[[i]])[2] + coef(fit[[i]])[3] * log2(x)**2)
    }
    sf  <- function(x,i) x**2/(max(x + 1e-6,vf(x,i)) - x)

    pv <- apply(data.frame(Meanz[[1]],Meanz[[2]]),1,
                function(x){
                  p12 <-  dnbinom("x" = 0:round(x[1]*length(cells_background) + x[2]*length(cells_of_interest),0), "mu" = mean(x)*length(cells_background), "size" = length(cells_background)*sf(mean(x),1))*
                    dnbinom("x" = round(x[1]*length(cells_background) + x[2]*length(cells_of_interest),0):0, "mu" = mean(x)*length(cells_of_interest), "size" = length(cells_of_interest)*sf(mean(x),2));
                  sum(p12[p12 <= p12[round(x[1]*length(cells_background),0) + 1]])/sum(p12)
                }
    ) # end of apply

    cat(paste(length(which(pv=="NaN")),"NaNs were replaced \n")) # Mauro's fix in these 3 lines
    cat(rownames(x)[which(pv=="NaN \n")])
    pv[which(pv=="NaN")] <-0.000000e+00

    res <- data.frame("baseMean" =mean(c(cells_of_interest, cells_background)),
                      "baseMean_BG"=Meanz[[1]],
                      "baseMean_target"=Meanz[[2]],
                      "foldChange"=Meanz[[2]]/Meanz[[1]],
                      "log2FoldChange"=log2(Meanz[[2]]/Meanz[[1]]),
                      "baseMedian_BG"=Medianz[[1]],
                      "baseMean_target"=Medianz[[2]],
                      "baseMean_target"=Medianz[[2]],
                      "pval"=pv,
                      "padj"=p.adjust(pv,method="BH"))
    vf1 <- data.frame("Meanz"=Meanz[[1]], "v"=v[[1]], "vm"=vf(Meanz[[1]],1))
    vf2 <- data.frame("Meanz"=Meanz[[2]], "v"=v[[2]], "vm"=vf(Meanz[[2]],2))
    rownames(res) <- rownames(vf1) <- rownames(vf2) <- rownames(x)
    list("vf1"= vf1, "vf2" = vf2, "res" = res)
  }
}


plotdiffgenesnb <- function(diffexpnb_object, pthr=.05, lthr=1, mthr=0, xname="A", yname="B", bgcol="grey", ppch = 20,
                            lcol = "grey33", show_names=TRUE, padj=TRUE, draw_mthr =T, draw_lthr =T, dontplotbelow=1, ...){
  y <- diffexpnb_object$res

  xlname =  paste0("log2 ((mRNA[",xname,"] + mRNA[",yname,"])/2)")
  ylname =  paste0("log2 (mRNA[",yname,"] / mRNA[",xname,"])")
  meanExpr = log2( (y$baseMeanA + y$baseMeanB)/2 )
  if (!is.na(dontplotbelow)) { meanExpr[meanExpr<dontplotbelow] = NA} # exlcude lowly expressed

  plot(meanExpr, y$log2FoldChange, pch=ppch, xlab=xlname , ylab=ylname, col=bgcol, ...)
  if (draw_mthr) { abline(v=mthr, lty=3, col=lcol) }
  if (draw_lthr) { abline(h=c(lthr, -lthr), lty=3, col=lcol)  }

  if ( ! is.null(pthr) ){
    if ( padj ) f <- y$padj < pthr else f <- y$pval < pthr
    points(meanExpr[f], y$log2FoldChange[f], col="red",pch=20)
  }
  if ( !is.null(lthr) ) f <- f & abs( y$log2FoldChange ) > lthr
  if ( !is.null(mthr) ) f <- f &  (log2( (y$baseMeanA + y$baseMeanB)/2 ) ) > mthr
  if ( show_names )  text(meanExpr[f],y$log2FoldChange[f],id2name(rownames(y))[f],cex=.5)

}
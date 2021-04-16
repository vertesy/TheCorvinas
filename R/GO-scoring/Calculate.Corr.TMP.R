
pc_TRUE(combined.obj@misc$expr.q90 >0)

# Calc.Cor.Seurat <- function(assay.use = "RNA", slot.use = "data"
#                             , quantileX = 0.95, max.cells =  40000, seed = p$"seed"
#                             , digits = 2, obj = combined.obj) {
#   expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
#   if (ncol(expr.mat) > max.cells) {
#     set.seed(seed = seed)
#     cells.use <- sample(x = colnames(expr.mat), size = max.cells)
#   }
#
#   qname = p0("q", quantileX * 100)
#   quantile_name = kpp("expr", qname)
#
#   if (is.null(obj@misc[[quantile_name]])) iprint("Call: combined.obj <- Calcq90Expression(combined.obj, quantileX =",quantileX," first )")
#   genes.HE = which_names(obj@misc[[quantile_name]] > 0)
#   iprint("Pearson correlation is calculated for", l(genes.HE), "HE genes with expr.",qname,": > 0.")
#   tic(); ls.cor <- sparse.cor(smat = t(expr.mat[genes.HE, cells.use])); toc()
#   ls.cor <- lapply(ls.cor, round, digits = 2)
#
#   slot__name <- kpp(slot.use, assay.use, quantile_name)
#   obj@misc[[kpp('cor', slot__name)]] <- ls.cor$'cor'
#   obj@misc[[kpp('cov', slot__name)]] <- ls.cor$'cov'
#   iprint("Stored under obj@misc$", kpp('cor', slot.use, assay.use), "or cov... ." )
#   return(obj)
# }

# combined.obj <- Calc.Cor.Seurat(assay.use = "RNA", slot.use = "data", digits = 2, obj = combined.obj, ...)

combined.obj <- Calcq90Expression(combined.obj, quantileX = 0.95, max.cells =  400000, set.all.genes = F)
memory.biggest.objects()
combined.obj <- Calc.Cor.Seurat(assay.use = "RNA", slot.use = "data"
                                , quantile = 0.95, max.cells = 10000, digits = 3)
print(object.size(combined.obj), units = "auto")
# idim(combined.obj@assays$RNA@scale.data)
# idim(combined.obj@misc$cor.data.RNA.expr.q95)

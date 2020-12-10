##################################################################################
# A custom R functions for one-liner ggplot similar to those in Markdownreports
##################################################################################
# try(source("~/GitHub/TheCorvinas/R/ggMarkdownreports.R"), silent = T)



stop("Do not use this script")
wA4 = 8.27 # A4 inches
hA4 = 11.69

## setup -------------------------------------------------------------------------------------------------
require(ggplot2)
# require(tibble)

## Themes -------------------------------------------------------------------------------------------------
#
# pie_theme <- theme_minimal() +
#  theme(
#   axis.title.x = element_blank(),
#   axis.title.y = element_blank(),
#   panel.border = element_blank(),
#   panel.grid = element_blank(),
#   axis.text.y = element_blank(),
#   axis.text.x = element_text(colour = 'black'),
#   axis.ticks = element_blank(),
#   plot.title = element_text(size = 14, face = "bold"),
#   axis.title = element_blank()
#  )


## qgpie -------------------------------------------------------------------------------------------------
"Now use ggpubr / ggExpress"
# ggpie <- function(dat, by, totals) {
#  # https://mathematicalcoffee.blogspot.com/2014/06/ggpie-pie-graphs-in-ggplot2.html
#  ggplot(dat, aes_string(x = factor(1), y = totals, fill = by)) +
#   geom_bar(stat = 'identity', color = 'black') +
#   guides(fill = guide_legend(override.aes = list(colour = NA))) + # removes black borders from legend
#   coord_polar(theta = 'y') +
#   pie_theme +
#   scale_y_continuous(breaks = cumsum(dat[[totals]]) - dat[[totals]] / 2, labels = dat[[by]])
# }
# # ggpie(dat = x, by = 'name', totals = 'value')
#
# qgpie <- function(NamedVector, percentage = TRUE, both_pc_and_value = FALSE,
#          plotname = substitute(NamedVector)) {
#  df = tibble::enframe(as.named.vector(NamedVector) )
#  ggpie(dat = df, by = 'name', totals = 'value') + ggtitle(plotname)
# }

# perGenotype = table(rep(LETTERS[1:3], 3)); qgpie(perGenotype)


## qgbar -------------------------------------------------------------------------------------------------

qgbar <- function(NamedVector, percentage = TRUE, both_pc_and_value = FALSE
         , plotname = substitute(NamedVector)
         , ylim = range(NamedVector)
         , ylab = substitute(NamedVector)
         , xlab = "names"
         , col = unless.specified("b.def.colors", "gold1")
) {
 df = tibble::enframe(NamedVector, name = xlab)

 ggplot(df) + geom_bar(stat = "identity", aes_string (x = (colnames(df)[1]), y = 'value'), fill = col) +
  ggtitle(plotname) +
  coord_cartesian(ylim = ylim) +
  labs(y = ylab) +
  theme_cowplot() +
  background_grid("y")
}

# (xtable <- table(rpois(100, 5)))
# qgbar(NamedVector = xtable)




## -------------------------------------------------------------------------------------------------

ggbar <- function(dat, by, totals) {
 # https://mathematicalcoffee.blogspot.com/2014/06/ggpie-pie-graphs-in-ggplot2.html
 ggplot(data = dat, aes_string(x = factor(1), y = totals, fill = by)) +
  geom_bar(stat = 'identity') +
  guides(fill = guide_legend(override.aes = list(colour = NA))) + # removes black borders from legend
  coord_polar(theta = 'y') +
  pie_theme +
  scale_y_continuous(breaks = cumsum(dat[[totals]]) - dat[[totals]] / 2, labels = dat[[by]])
}
# ggpie(dat = x, by = 'name', totals = 'value')


# qsave_plot.A4v2(filename = "UMAPs.batches.and.samples.png", plot = p1, base_height = 12, ncol = 1, nrow = 1) #Figure 2
# qsave_plot.A4h4
# qsave_plot.A4v2
## -------------------------------------------------------------------------------------------------
# Currently in Seurat.utils
# qqsave <- function(ggplot.obj, h = 7, PNG = F, title = NULL, plotit = F) { # Quickly save a ggplot object, and optionally display it
#  pname = substitute(ggplot.obj)
#  fname = ww.FnP_parser(pname, if (PNG) "png" else "pdf")
#  save_plot(filename = fname, plot = ggplot.obj, base_height = h) #, ncol = 1, nrow = 1
#  if (plotit) ggplot.obj
# }


## -------------------------------------------------------------------------------------------------
# Currently in Seurat.utils
# qqSaveGridA4 <- function(plotlist = pl # Save 2 or 4 ggplot objects using plot_grid() on an A4 page
#             , plots = 1:2, NrPlots = length(plots), height = hA4, width = wA4
#             , fname = "Fractions.Organoid-to-organoid variation.png") {
#  stopifnot(NrPlots %in% c(2,4))
#  iprint(NrPlots,"plots found,", plots,"are saved.")
#  pg.cf = plot_grid(plotlist = plotlist[plots], nrow = 2, ncol = NrPlots/2, labels = LETTERS[1:NrPlots] )
#  if (NrPlots  == 4) list2env(list(height = width, width = height), envir = as.environment(environment()))
#  save_plot(filename = fname,
#       plot = pg.cf, base_height = height, base_width = width)
#  ww.FnP_parser(fname)
# }
# # qqSaveGridA4(plotlist = pl, plots = 1:2, fname = "Fractions.per.Cl.png")
# # qqSaveGridA4(plotlist = pl, plots = 1:4, fname = "Fractions.per.Cl.4.png")


## -------------------------------------------------------------------------------------------------
# p <- ggplot(df, aes(x = weight)) +

# p

## qghist -------------------------------------------------------------------------------------------------
qghist <- function(NamedVector, percentage = TRUE, both_pc_and_value = FALSE
          , nrbins = 40
          , plotname = substitute(NamedVector)
          # , ylim = range(NamedVector)
          , ylab = substitute(NamedVector)
          , xlab = "names"
          , logYaxis = TRUE
          , vline = 10
          , filtercol = 1
          , col = unless.specified("b.def.colors", "gold1")
          , savefile = unless.specified("b.save.wplots")
) {
 df = tibble::enframe(NamedVector)


 p <- ggplot(df) +
  # geom_histogram(stat = "identity", color = "black", fill = "white") +
  geom_histogram(aes_string(x = 'value'), bins = nrbins, color = "black", fill = col) +
  ggtitle(plotname) +
  # coord_cartesian(ylim = ylim) +
  labs(y = ylab) +
  theme_cowplot() +
  background_grid("y")
 if (logYaxis) p <- p + scale_y_log10()
 if (vline) p <- p + geom_vline(xintercept = vline, linetype = "dashed", color = "red", size = 1)
 # if (filtercol) p <- p + 1
 save_plot(filename = ppp(plotname, "pdf"), plot = p)
 p
}
# qghist(NamedVector)
# ## -------------------------------------------------------------------------------------------------
# qqSaveA4p <- function(ggobj, ext =c("png", "pdf")[1]) {
#   title = substitute(ggobj)
#   save_plot(plot = ggobj, filename = kpp(title, ext), base_height = hA4, base_width = wA4)
# }
#
#
# qqSaveA4l <- function(ggobj, ext =c("png", "pdf")[1]) {
#   title = substitute(ggobj)
#   save_plot(plot = ggobj, filename = kpp(title, ext), base_height = wA4, base_width =hA4 )
# }
#
#
# qqSaveA5l <- function(ggobj, ext =c("png", "pdf")[1]) {
#   title = substitute(ggobj)
#   save_plot(plot = ggobj, filename = kpp(title, ext), base_height = hA4/2, base_width = 8.27 )
# }

# qqSave <- function(ggobj, ext =c("png", "pdf")[1], w =5, h = w) {
#   title = substitute(ggobj)
#   save_plot(plot = ggobj, filename = kpp(title, ext), base_height = w, base_width = h)
# }


## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------

## Functions for tidyverse interaction -------------------------------------------------------------------------------------------------

# col2named.vec.tbl <- function(tbl.2col) {
#  nvec = tbl.2col[[2]]
#  names(nvec) = tbl.2col[[1]]
#  nvec
# }
#
#
#
#
#

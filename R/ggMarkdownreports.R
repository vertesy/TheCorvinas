##################################################################################
# A custom R functions for one-liner ggplot similar to those in Markdownreports
##################################################################################
# try(source("~/GitHub/TheCorvinas/R/ggMarkdownreports.R"), silent = T)



wA4 = 8.27 # A4 inches
hA4 =11.69


## setup -------------------------------------------------------------------------------------------------
require(ggplot2)
# require(tibble)

## Themes -------------------------------------------------------------------------------------------------

pie_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_text(colour='black'),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    axis.title=element_blank()
  )


## qgpie -------------------------------------------------------------------------------------------------

ggpie <- function (dat, by, totals) {
  # https://mathematicalcoffee.blogspot.com/2014/06/ggpie-pie-graphs-in-ggplot2.html
  ggplot(dat, aes_string(x=factor(1), y=totals, fill=by)) +
    geom_bar(stat='identity', color='black') +
    guides(fill=guide_legend(override.aes=list(colour=NA))) + # removes black borders from legend
    coord_polar(theta='y')  +
    pie_theme +
    scale_y_continuous(breaks=cumsum(dat[[totals]]) - dat[[totals]] / 2, labels=dat[[by]])
}
# ggpie(dat = x, by = 'name', totals = 'value')

qgpie <- function(NamedVector, percentage = TRUE, both_pc_and_value = FALSE,
                  plotname = substitute(NamedVector)) {
  df = tibble::enframe(NamedVector )
  df = tibble::enframe(NamedVector)
  ggpie(dat = x, by = 'name', totals = 'value') + ggtitle(plotname)
}
# perGenotype = table(rep(LETTERS[1:3], 3)); qgpie(perGenotype)


## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------------

## Functions for tidyverse interaction -------------------------------------------------------------------------------------------------

col2named.vec.tbl <- function(tbl.2col) {
  nvec = tbl.2col[[2]]
  names(nvec) = tbl.2col[[1]]
  nvec
}






#
#
# lapply(df.MapData.per.Worm, l)
# idim(df.MapData.per.Worm)
# xsapply(words[1:12], str_split, pattern=",")
require(stringr)
datMAT=matrix(data = words[1:12], nrow = 3)

fred <- sapply(datMAT, str_split, pattern="")
# fred <- as.list(1:12)


dim(fred) <- c(3, 4); fred
dimnames(fred) = list(letters[1:3], LETTERS[1:4])
fred[['a','A']]

lapply(fred,unlist)
sapply(fred,unlist)[2,2]
arr.of.lists.2.df(mdlapply(fred,l))


copy.dimension.and.dimnames <- function(list.1D, obj.2D) {  # copy dimension and dimnames
  dim(list.1D) <- dim(obj.2D)
  dimnames(list.1D) <- dimnames(obj.2D)
  list.1D
  }


mdlapply <- function(list_2D, ...) {  # multi dimensional lapply
  x = lapply(list_2D, ...)
  copy.dimension.and.dimnames(x,list_2D)
 }


mdlapply(fred, l)

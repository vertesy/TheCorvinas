# Object Oriented Programming in R
# http://adv-r.had.co.nz/S3.html



x <- 1
attr(x, "class") <- "foo"
# Or in one line
x <- structure(1, class = "foo")



letters
x <- structure(1, class = letters)

bar <- function(x) UseMethod("bar", x)
bar.z <- function(x) "it is z"
bar(x)
# [1] "z"


bar.x <- function(x) "it is XXX"

# You can call methods directly, but you shouldn't!
bar.x(x)
# [1] "x"
bar.z(x)
# [1] "z"

# ------------------------------------------------------------------------------------------------------------
# NextMethod

baz <- function(x) UseMethod("baz", x)
baz.A <- function(x) "A"
baz.B <- function(x) "B"

ab <- structure(1, class = c("A", "B"))
ba <- structure(1, class = c("B", "A"))
baz(ab)
baz(ba)

baz.C <- function(x) c("C", NextMethod())

ca <- structure(1, class = c("C", "A"))
cb <- structure(1, class = c("C", "B"))

baz(ca)
baz(cb)


# ------------------------------------------------------------------------------------------------------------

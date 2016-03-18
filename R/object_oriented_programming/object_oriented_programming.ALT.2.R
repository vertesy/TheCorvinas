# Object Oriented Programming in R


xx <- rnorm(1000)
class(xx)
plot(xx)
yy <- ecdf(xx)
class(yy)
plot(yy)
plot
plot.ecdf
plot.default
methods("plot")
getS3method("plot", "histogram")

"What plot does, depends on the class of the x argument. It is a method. plot.ecdf is the ecdf method for plot."



## Constructing a new S3 Class

jim <- list(height = 2.54 * 12 * 6/100, weight = 180/2.2, name = "James")
jim
class(jim) <- "person2"
class(jim)
jim

"We have now made an object of class person. We now define a print method."


print(jim)
print.person2 <- function(x) {
	cat("name:", x$name, "\n")
	cat("height:", x$height, "meter s", "\n")
	cat("weight:", x$weight, "kilograms", "\n")
}
print(jim)


'Note the method/class has the "dot" naming convention of method.class.'

## S3 classes are not robust

fit <- lm(rnorm(100) ~ 1)
class(fit)
print(fit)
class(fit) <- "something"
print(fit)
class(fit) <- "person"
print(fit)


"In case print does not have a method for the class, it dispatches to the default method, print.default."
"S3 does not have the concept of type checking { there is no way to formally define a class and ensure that the object conform to the definition."

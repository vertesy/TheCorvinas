# object oriented programming in R

bubba <- list(first="one", second="two", third="third")
class(bubba) <- append(class(bubba),"Flamboyancy")

bubba


	GetFirst <- function(x)
		+ {
			+     UseMethod("GetFirst",x)
			+ }
>
	> GetFirst.Flamboyancy <- function(x)
		+ {
			+    return(x$first)
			+ }
>
	> GetFirst(bubba)


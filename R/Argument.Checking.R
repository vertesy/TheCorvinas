install.packages('assertthat')

identifier <- function(name, version, seed) {
  
  # throw error if any of the following assertions fail
  stopifnot(
    length(name) == 1,    # name must be a scalar
    length(version) == 1, # version must be a scalar
    length(seed) == 1,    # seed must be a scalar
    rlang::is_integerish(version),  # version must be a whole number
    rlang::is_integerish(seed),     # seed must be a whole number
    !stringr::str_detect(name, "[[:space:]._]"), # name can't have spaces, periods, or underscores 
    seed > 0,      # seed must be positive
    seed < 10000,  # seed must be less than 10000
    version > 0,   # version must be positive
    version < 100  # version must be less than 100
  )
  
  # the actual work of the function
  version <- stringr::str_pad(version, width = 2, pad = "0")
  seed <- stringr::str_pad(seed, width = 4, pad = "0")
  paste(name, version, seed, sep = "_") 
}


rnd4l <- function(set = c(LETTERS, 0:9), n = 4) {
  print(paste0(paste0( sample(x = set, size = n), collapse = ''), '__'))
}
rnd4l()


identi <- function(name) {
  stopifnot( 'Not sure why' = length(name) == 2, length(name) == 3)
  print(name)
}

identi(1:2)
identi(1)


# ----------------------------------------------------------------------------------------------------
#  checkmate
# ----------------------------------------------------------------------------------------------------
require(checkmate)


# Define a function that takes a list as an argument
my_function <- function(x) {
  # Check if x is a list with at least 2 elements
  assert_list(x, max.len = 2, names = checkSubset(x, c(LETTERS)) )
              
  # Check if the first element of x is a numeric vector
  assert_numeric(x[[1]])
  # assert_character(x[[2]])
            
              
}

# Call the function with valid arguments
nana <-  list( 'a' = 1:10, 'b' = c("a", "b"))
my_function(x = nana)

# Call the function with invalid arguments
my_function(list("hello", 2))



my_function <- function(x, y) {
  assertNumeric(x,lower = 2, )
  assertNumeric(y,lower = 1)
  
  
  # Do something with x and y
  x*y
}

# Call the function with valid arguments
my_function(x = 3:2, y = 2)




# ----------------------------------------------------------------------------------------------------
#  assertthat
# ----------------------------------------------------------------------------------------------------
require(assertthat)

# Define a function that takes two arguments
my_function <- function(x, y) {
  # Check if x is a numeric vector
  assert_that(length(x) == 1, is.numeric(x), is.numeric(y))
  
  # Check if y is a character vector
  assert_that(is.character(y))
  
  # Do something with x and y
  x*y
}

# Call the function with valid arguments
nanana <- c("a", "b", "c")
my_function(1:10, y = nana)

# Call the function with invalid arguments
my_function("hello", 2)



# ----------------------------------------------------------------------------------------------------
#  assertthat
# ----------------------------------------------------------------------------------------------------



library(checkmate)

# Define a function that takes two arguments
my_function <- function(x, y) {
  # Check if x and y are overlapping
  assert( max(x[1], y[1]) <= min(x[2], y[2]), "The inputs are not overlapping" )
  
  # Do something with x and y
}

# Call the function with overlapping inputs
my_function(c(1, 5), c(3, 7))

# Call the function with non-overlapping inputs
my_function(c(1, 5), c(6, 10))

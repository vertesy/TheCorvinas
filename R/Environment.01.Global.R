# Environment.01.Global.R

# env_frame <- new.env()
my_env <- new.env(parent = baseenv())


rm('y'); y
x <- 1
my_env$x <- x
source('~/GitHub/TheCorvinas/R/Environment.02.Local.R', local = my_env) 
y
my_env$y


# Now when you attach my_env, it will not see .GlobalEnv
attach(my_env)
y
x
z=2
detach(my_env)
y
z


# ----
baseenv() <- my_env
my_env$parent_env <- .GlobalEnv
rm(list = ls(), envir = .GlobalEnv )
x
y
z
attach(my_env, pos = 3L)
x
y
z

print(get("y", envir = my_env))

# Environment.01.Global.R

rm(list = ls(all.names = TRUE))

ls()
# env_frame <- new.env()
my_env <- new.env(parent = baseenv())
?search()
my_env <- .GlobalEnv

rm('y'); y
x <- 1
my_env$x <- x

my.script <- '~/GitHub/Packages/isoENV/Examples/Environment.02.Local.R'
source(file = my.script, local = my_env) 
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



rm(list = setdiff(ls(envir = .GlobalEnv), "my_env"), envir = .GlobalEnv)
attach(my_env)
x
rm(list = setdiff(ls(), lsf.str()))

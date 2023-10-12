print("Loading Abel's .Rprofile")

### Custom Options -------------------------------------------------------
options(bitmapType = 'cairo') # otherwise png does not work on CBE.
options(max.print = 150)

### Custom Global Variables -------------------------------------------------------
onCBE <- TRUE # A global env. variable telling R that you are on CBE

### Custom Functions -------------------------------------------------------
gggcon <- function() print('~/GitHub/Projects/CON/');
oo <- function() {print(list.files(getwd())); print(getwd())}

### Custom Packages -------------------------------------------------------
try(require("colorout"), silent = TRUE)
# install_github(repo = "jalvesaq/colorout")

try(source('https://raw.githubusercontent.com/vertesy/Rocinante/main/R/Rocinante.R'), silent = T)




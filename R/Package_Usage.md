# [How many people use a package?](http://www.r-statistics.com/2013/06/answering-how-many-people-use-my-r-package/)


```
# if (!require('devtools')) install.packages('devtools'); require('devtools')
# make sure you have Rtools installed first! if not, then run:
#install_Rtools()
#install_github('installr', 'talgalili') # get the latest installr R package
# or run the code from here:
# https://github.com/talgalili/installr/blob/master/R/RStudio_CRAN_data.r
 
if(packageVersion("installr") %in% c("0.8","0.9","0.9.2")) install.packages('installr') #If you have one of the older installr versions, install the latest one....
 
require(installr)
 
# The first two functions might take a good deal of time to run (depending on the date range)
RStudio_CRAN_data_folder <- download_RStudio_CRAN_data(START = '2013-04-02', END = '2013-04-05') # around the time R 3.0.0 was released
my_RStudio_CRAN_data <- read_RStudio_CRAN_data(RStudio_CRAN_data_folder)
 
 # barplots: (more functions can easily be added in the future)
barplot_package_users_per_day("plyr", my_RStudio_CRAN_data)
barplot_package_users_per_day("installr", my_RStudio_CRAN_data)
```
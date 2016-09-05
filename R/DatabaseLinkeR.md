
# Parse links to databases from your list of  gene symbols

[Get DatabaseLinke.R](https://github.com/vertesy/TheCorvinas/blob/master/R/DatabaseLinkeR)


### Install

1.) Download and source
https://github.com/vertesy/TheCorvinas/blob/master/R/DatabaseLinke.R
2.) Source from the web
source(https://raw.githubusercontent.com/vertesy/TheCorvinas/master/R/DatabaseLinke.R)
Read: https://github.com/vertesy/TheCorvinas/blob/master/R/DatabaseLinkeR.md




## Components

1. R script containing functions each of which operates on vectors of gene symbols.
2. An executable bash script that you need to make. You can make it anywhere, but you need to specify it in the R-script's `BashScriptLocation` variable.


## Usage

You can use the functions in 3 ways:


`link_String("Mecom")`

Open the link in your web browser.

`link_String("Mecom", writeOut = T)`

Writes the link in (an executable) bash script run.sh (if you have too many links).

It will actually write the links to `BashScriptLocation` in an executable format, so once you run the script, it opens all the links in your default browser (on OS X).

Using `link_String("Mecom", writeOut = T)`, in your bash script, you will find:

```
open 'http://string-db.org/newstring_cgi/show_network_section.pl?identifier=Tuba3a&species=10090'
open 'http://string-db.org/newstring_cgi/show_network_section.pl?identifier=Piwil1&species=10090'
```


`link_String("Mecom", writeOut = F, Open=F)`

Writes the link to the screen. You get:

```
http://string-db.org/newstring_cgi/show_network_section.pl?identifier=Tuba3a&species=10090
http://string-db.org/newstring_cgi/show_network_section.pl?identifier=Piwil1&species=10090
```

as a character vector, so you can write out in a column of your gene-table.


#### For detailed usage, browse the code of the functions.

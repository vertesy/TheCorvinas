
# Parse links to databases from your list of  gene symbols

[Get DatabaseLinke.R](https://github.com/vertesy/TheCorvinas/blob/master/R/DatabaseLinkeR)



## Components

1. R script containing functions each of which operates on vectors of gene symbols.
2. An executable bash script that you need to make. You can make it anywhere, but you need to specify it in the R-script's `BashScriptLocation` variable.


## Usage

Each function work like this:

`link_String(YourGeneSymbols)`

And it will write the links `BashScriptLocation` in an executable format, so once you run the script, it opens all the links in your default browser (on OS X).

In your bash script, you find:

```
open 'http://string-db.org/newstring_cgi/show_network_section.pl?identifier=Tuba3a&species=10090'
open 'http://string-db.org/newstring_cgi/show_network_section.pl?identifier=Piwil1&species=10090'
```


If you specify `writeOut=FALSE`, it will return the list of links in the terminal

`link_String(YourGeneSymbols, writeOut=FALSE)`

You get:

```
http://string-db.org/newstring_cgi/show_network_section.pl?identifier=Tuba3a&species=10090
http://string-db.org/newstring_cgi/show_network_section.pl?identifier=Piwil1&species=10090
```

as a character vector, so you can write out in a column of your gene-table.


#### For detailed usage, browse the code of the functions.

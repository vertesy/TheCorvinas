
# Parse links to databases from your list of  gene symbols

[Get the DatabaseLinke.R](https://github.com/vertesy/TheCorvinas/blob/master/R/DatabaseLinke.R)

### Summary

A set of functions to parse links to genomics related websites.

It can do 3 things:

1. Parse (and store) the links for a vector of gene symbols
2. Parse the links, and **open** them **in the default browser** (beware for 20+ genes)
3. Parse the link, and **write** `open` commands **into an** (executable) **bash script**, which you can later open in batches.

*Note: It typically parses query (search) links, instead direct links*

### List of Databases

- Ensemble
- Uniprot
- STRING: a protein interaction database
- Pubmed
- Wikipedia

#### Species specific
- Wormbase
- HGNC (Human Gene Nomenclature)

### Install

1.) Download and source: https://github.com/vertesy/TheCorvinas/blob/master/R/DatabaseLinke.R


2.) Source from the web:

`source("https://raw.githubusercontent.com/vertesy/TheCorvinas/master/R/DatabaseLinke.R")`

Read: https://github.com/vertesy/TheCorvinas/blob/master/R/DatabaseLinkeR.md

## Components

1. R script containing functions each of which operates on vectors of gene symbols.
2. An executable bash script that you need to make. You can make it anywhere, but you need to specify it in the R-script's `BashScriptLocation` variable.


## Usage
*Default species is typically mice*

You can use the functions in 3 ways:

1. Open the link in your web browser:

   `link_String("Mecom")`

2. Writes the link in a bash script, called run.sh:   `link_String("Mecom", writeOut = T)`
  - use, if  you have too many links to open at once
  - When you run the script, it opens all the links as tabs in your default browser (on OS X /  *nix).
  - **Comment out some lines if its too much**
  - More precisely, it will write the links to `BashScriptLocation` in an executable format.

### Writing to file

Using `link_String("Mecom", writeOut = T)`, in your bash script, you will find:

```
open 'http://string-db.org/newstring_cgi/show_network_section.pl?identifier=Mecom&species=10090'
```


### Print links to the screen

`link_String("Mecom", writeOut = F, Open=F)`

Writes the link to the screen. You get:

```
http://string-db.org/newstring_cgi/show_network_section.pl?identifier=Mecom&species=10090
```

as a character vector, so you can write out in a column of your gene-table.



#### [Google search URL / search query](http://www.our-picks.com/archives/2007/01/30/google-search-urls-revealed-or-how-to-create-your-own-search-url/)

| String        | Meaning                                  |
| ------------- | ---------------------------------------- |
| as_oq         | This tells Google to find pages in which at least instance of nintendo **OR**  wii is found |
| as_q          | This means that you look for **both** nintendo and wii in the same page |
| as_epq        | Google translates this as a Google search of “nintendo wii”, searches the  exact phrase ‘nintendo wii’ |
| num           | The number of results you want displayed, it ranges from 0 to 100. If you  set num to 0 you will get the ‘No match found” message |
| safe          | If you set this to active the Google Safe Search is on and the adult  material will be filtered |
| as_eq         | Use this to exclude a term from your search |
| as_qdr        | Shows only results that have been updated in the given time interval.  Possible values: y (year), m6 (6 months), m3 (3 months). |
| as_sitesearch | Limits the search to a specific domain or TLD (.us; .gov; .co.uk; .ro;  etc) |
| as_occt       | This is set by default to ‘any’ but if you change it you can search in:  title, url, links |





---------




### If you encounter a bug, or something doesn't work,

1. please let me know by raising an issue on [TheCorvinas](https://github.com/vertesy/TheCorvinas/issues/new?milestone=DatabaseLinke.R)
2. Fix by providing the missing function.
      - Try to `source("https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R")`
        [See details for CodeAndRoll](https://github.com/vertesy/TheCorvinas/blob/master/R/CodeAndRoll.md)
      - If still some functons are missing, try to install MarkdownReports:
        `install.packages("devtools")`
        `devtools::install_github(repo = "vertesy/MarkdownReports/MarkdownReports")`
        `require("MarkdownReports")`
        [See details for MarkdownReports](https://vertesy.github.io/MarkdownReports/)





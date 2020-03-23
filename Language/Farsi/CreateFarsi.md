## Conversion pipeline of the most frequent farsi words to ANKI flash cards


1. #### Get the N most frequent farsi word you want to use from [Behnam Esfahbod's Persian Words Frequency Database](https://github.com/behnam/persian-words-frequency)

	I used the 1000 most frequent words from [wikipedia (6mb file!)](https://github.com/behnam/persian-words-frequency/blob/master/persian-wikipedia.txt)

- #### Paste into [Google Translate](https://translate.google.com/#fa/en/) and add translation column to your tab separated file.


- #### [Transliterate/Romanise](http://mylanguages.org/farsi_romanization.php) the original farsi words on mylanguages.org

- #### After romanisation *reverse* the letters in each word!
	Because Farsi is read right to left. *See R script at the end.*

- #### Now go to [cleanup](http://mylanguages.org/romanization_cleanup.php) on mylanguages.org
	Manual checking of the results would be great. Without this, its very approximate.

- #### Paste back to the .tsv file and edit it

	- #### Use `md.csv_to_md_parser.py` shell function to  create a markdown reference table`most.freq.words.in.persian-wikipedia.UTF8.tsv.md`

		See: [md.csv_to_md_parser.py](/Users/abelvertesy/Github_repos/TheCorvinas/Python/MarkdownParsers/md.tsv_to_md_parser.py)

		Optionally, in bash:
		`alias tsv2markdown='/Users/XXXX/Github_repos/TheCorvinas/Python/MarkdownParsers/md.tsv_to_md_parser.py' `

		- Add [wiktionary](https://www.wiktionary.org/) links. These are parsed for each word, like:
		`https://en.wiktionary.org/wiki/`ﻮ`#Persian`


- #### [View](most.freq.words.in.persian-wikipedia.UTF8.tsv.md)

- #### Formatting before importing into Anki
	See code at the end.

- #### Import into 4 decks of ~240 words to [Anki](https://apps.ankiweb.net/) and export.




------------------------------------------------------------------------------------------------------------------------------------------------

### R Code for reversing letters in each word

```
source("https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R") # get fromClipboard.as_vec() for OS X
library(IRanges) # for reverse()

fs =fromClipboard.as_vec() # for non-OS X use read() from a file with your words
sf = reverse(fs)
toClipboard(sf)
```

### R Code: formatting for Anki

Use R's make.unique() of the latin words

```
en = make.unique(fromClipboard.as_vec())
toClipboard(en))
```
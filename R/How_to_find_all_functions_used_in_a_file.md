# How_to_find_all_functions_used_in_a_file
> /Users/abelvertesy/TheCorvinas/R/How_to_find_all_functions_used_in_a_file.md


1. Open Script
2. in **Sublime2** search for `[a-zA-Z]+\(` with regex, select all
3. copy into a new file / to clipboard.
4. Process in R & search for funs 
			
	```
	inline_vec.char.from_Clipboard() # from CodeAndRoll
	funz= c( 'kollapse', ... blablabla ...'kollapse', 'l', 'write')
	funz2 = sort(table(funz), decreasing = T)
	funz2

	```
- Find your custom functions, and paste them back in the file

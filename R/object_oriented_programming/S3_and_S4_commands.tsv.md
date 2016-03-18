| S3 	| S4 	| Function |
|---|---|---|
| class (x) <- yourClass 	| setClass (x, "class") 	| Define the class of an object |
| yourFunction.yourClass <- function(x, ...) 	| setMethod (x, "method") 	| Define the method belonging to a certain class of objects. |
| function (yourFunction) { UseMethod("yourFunction") } 	| setGeneric( yourFunction, standardGeneric("yourFunction") ) 	| Define a generic function |
| UseMethod("yourFunction") 	| standardGeneric("yourFunction") 	| Dispatching a method (inside a generic function). |

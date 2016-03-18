# Object oriented programming in R

1. You define a `class` for a `variable`. A variable with a class / classes is an `object`.
	> The class will be sort of a *tag* that tells the functions, what kind of variable, with what kind of dimensions should it expect.
2. You define a `method` for your function
	> A method is sort of a *specialized function* that makes assumptions based on the `class` of the object.
3. You call your function on your `object`,
	> You either directly call your **method** (specialized function), or your let a **generic function** decided which method to call.


Or, as [Hadley says](http://adv-r.had.co.nz/S3.html)

*"Central to any object-oriented system are the concepts of class and method. 
A **class** defines a type of object, describing what properties it possesses, how it behaves, and how it relates to other types of objects. 
Every **object** must be an instance of some class. 
A **method** is a function associated with a particular type of object."*

...*"[In S3,] a special type of function called a **generic function decides which method to call.**"*

...*"**R looks for methods in the order** in which they appear in the class vector."*



```
x <- structure(1, class = "foo")

mean <- function (x, ...) {
   UseMethod("mean", x)
 }

mean(x) # x of class foo 
# Equals to
mean.foo(x)

```

## Comparison S3 /S4


| S3 	| S4 	| Function |
|---|---|---|
| class (x) <- yourClass 	| setClass (x, "class") 	| Define the class of an object |
| yourFunction.yourClass <- function(x, ...) 	| setMethod (x, "method") 	| Define the method belonging to a certain class of objects. |
| function (yourFunction) { UseMethod("yourFunction") } 	| setGeneric( yourFunction, standardGeneric("yourFunction") ) 	| Define a generic function |
| UseMethod("yourFunction") 	| standardGeneric("yourFunction") 	| Dispatching a method (inside a generic function). |
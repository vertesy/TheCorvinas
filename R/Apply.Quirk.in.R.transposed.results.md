In R, the `apply()` function is used to apply a function across **rows** or **columns** of a matrix. The behavior depends on the `MARGIN` argument:

1. `MARGIN = 1`: apply the function **row-wise**
2. MARGIN = 2`: apply the function **column-wise**

However there is a pretty unexpected quirk in what they return

1.  `apply(..., MARGIN = 1)` returns a **transposed result**, while 
2. `apply(..., MARGIN = 2)` does **not**. 

This means that **when applying a function row-wise, the result matrix has rows and columns flipped unless corrected with `t()`.** 
Column-wise application behaves as expected and preserves shape.


Here’s a simple example to demonstrate:

```r
m <- matrix(1:9, nrow = 3, byrow = TRUE)
#      [,1] [,2] [,3]
# [1,]    1    2    3
# [2,]    4    5    6
# [3,]    7    8    9
```

Row-wise (`MARGIN = 1`):

```r
apply(m, 1, function(x) x + 1)
#      [,1] [,2] [,3]
# [1,]    2    5    8
# [2,]    3    6    9
# [3,]    4    7   10
```

→ Output is **transposed** — each row result becomes a column. To fix:

```r
t(apply(m, 1, function(x) x + 1))
```

Column-wise (`MARGIN = 2`):

```r
apply(m, 2, function(x) x + 1)
#      [,1] [,2] [,3]
# [1,]    2    3    4
# [2,]    5    6    7
# [3,]    8    9   10
```

→ Output is in expected shape — no transposition needed.

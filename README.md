# `disize`: A tool for size factor estimation

Leverage experimental design for size factor estimation by providing a `design_formula` to `disize`.

Currently, `disize` accepts either a count matrix `counts` and a `metadata` dataframe containing observation-level metadata(the predictors in your design!), or an already formatted `model_data` dataframe.

Note: the rows of `counts` should have the same names as the `obs_name` column in `metadata`. If `counts` has no row names, `disize` assumes it is ordered such that the row indices correspond to the row indices of `metadata`.

# Examples
```R
> dim(counts)
[1] 900 500

> metadata
# A tibble: 900 × 4
   obs_id batch_id donor cell_type
   <fct>  <fct>    <fct> <fct>
 1 1      1        1     1
 2 2      1        1     2
 3 3      1        1     3
 4 4      1        1     1
 5 5      1        1     2
 6 6      1        1     3
 7 7      1        1     1
 8 8      1        1     2
 9 9      1        1     3
10 10     1        1     1
# ℹ 890 more rows
# ℹ Use `print(n = ...)` to see more rows

> (size_factors <- disize(
+     design_formula = ~ cell_type + (1 | cell_type:donor),
+     counts = counts,
+     metadata = metadata
+ ))
            1             2             3             4             5            6
 0.0020845603  0.0013157564 -0.0025107518 -0.0025500306  0.0006061994  0.0010440993
```

# TODO
- Figure out multithreading

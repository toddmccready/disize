# **disize**: A tool for size factor estimation

Leverage experimental design for size factor estimation by providing a `design_formula` to `disize`.

Currently, `disize` accepts either a count matrix `counts` and a `metadata` dataframe containing observation-level metadata(the predictors in your design!), or an already formatted `model_data` dataframe.

Note: the rows of `counts` should have the same names as the `obs_name` column in `metadata`. If `counts` has no column names, `disize` assumes it is ordered such that the column indices correspond to the row indices of `metadata`.

# Examples
```r
> data
# A tibble: 30,000 × 6
    gene donor batch cell_type cell_barcode counts
   <int> <int> <int>     <int> <chr>         <dbl>
 1     1     1     1         1 1:1_1             4
 2     1     1     1         1 1:1_2             8
 3     1     1     1         1 1:1_3             8
 4     1     1     1         1 1:1_4             7
 5     1     1     1         1 1:1_5             6
 6     1     1     1         1 1:1_6            10
 7     1     1     1         1 1:1_7            11
 8     1     1     1         1 1:1_8             4
 9     1     1     1         1 1:1_9            17
10     1     1     1         1 1:1_10            7
# ℹ 29,990 more rows
# ℹ Use `print(n = ...)` to see more rows
> (size_factors <- disize(
+     design_formula = ~ cell_type + (1 | donor:cell_type),
+     model_data = data,
+     n_threads = 7,
+     verbose = FALSE
+ ))
         1          2
-0.9821126  0.4858037
```

# TODO
- Do custom convergence checks for size factor parameters, essentially compile the model once and use the "init" argument for rstan::optimization to optimize in shorter passes, then add custom stopping conditions for the size factors.
- Offer downsampling (i.i.d) observations as we don't care about precision of expression estimates.

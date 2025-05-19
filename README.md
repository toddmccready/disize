# **disize**: A tool for size factor estimation

Thought of a way to use sparsity-inducing priors to leverage experimental design when estimating size factors.

The only function exported is `disize::disize`, provide your count matrix and a dataframe providing metadata for the samples (or pass a explicit `model_data` argument) and you're good to go.

# Examples
```r
> data
# A tibble: 200 × 5
    gene donor batch cell_barcode counts
   <int> <int> <int> <chr>         <dbl>
 1     1     1     1 1_1               5
 2     1     2     2 2_1               7
 3     2     1     1 1_1               0
 4     2     2     2 2_1               4
 5     3     1     1 1_1               2
 6     3     2     2 2_1              16
 7     4     1     1 1_1               0
 8     4     2     2 2_1               5
 9     5     1     1 1_1               1
10     5     2     2 2_1               3
# ℹ 190 more rows
# ℹ Use `print(n = ...)` to see more rows
> (size_factors <- disize(design_formula = ~ (1 | donor), model_data = data))
sf_batch1  sf_batch2
-0.2957847  0.2279738
```

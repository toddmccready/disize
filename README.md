# disize: A tool for size factor estimation

Leverage information from the experimental design during size factor estimation.

# Installation

As `disize` is not yet on CRAN, installation is not a one-liner with `install.packages`:

## With `remotes`
```R
# Install disize
remotes::install_github("https://github.com/toddmccready/disize")

# Set up CmdStan toolchain
cmdstanr::install_cmdstan()
```

## With `rv`

Add the following entries to your `rproject.toml` file(if not already present):
```
repositories = [
    # ...
    { alias = "STAN", url = "https://stan-dev.r-universe.dev" },
    # ...
]

# ...

dependencies = [
    # ...
    { name = "disize", git = "https://github.com/toddmccready/disize", tag = "v0.4.14" },
    # ...
]
```

Then sync your lockfile:
```sh
rv sync
```

And finally install the CmdStan toolchain in `R`:
```R
cmdstanr::install_cmdstan()
```

# Usage

The only export is `disize::disize`! The required arguments are:

- `design_formula`: An R formula that specifies the experimental design. This is the same R formula you would pass to something like `DESeq2` or `edgeR` including predictors like `condition`, `sex`, etc to estimate your expression quantities of interest(except we allow for random-effects). All terms used in this formula should be present in `metadata`.

- `counts`: A (observation x features) matrix containing the transcript counts. This can be dense or sparse; an internal coercion to a dense matrix will be done after subsetting relevant features and observations.

- `metadata`: A dataframe containing the observation-level information(with observations as rows).

The batch- and observation identifiers are specified by `batch_name` and `obs_name`, respectively. Ensure that your batch identifier is specified as a column in `metadata`, however `obs_name` is not required if the row indices of `counts` and `metadata` correspond to the same observation(i.e., if either `rownames(counts)` or `metadata[[obs_name]]` is `NULL` then the row indices of `counts` and `metadata` are assumed to correspond to the same observation).


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

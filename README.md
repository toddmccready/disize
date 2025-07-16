# disize

The existing methods for RNAseq normalization are
`DESeq2`'s [median of ratios](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106) (MoR)
and `edgeR`'s [trimmed mean of M values](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25) (TMM) .
These methods however do not include information about the experimental design when trying to estimate size factors, and can fail for more complex study designs.
**D**esign **i**nformed **size** factor estimation (or `disize`) is an alternative normalization method that [jointly models gene expression and batch-effects](https://toddmccready.github.io/disize/articles/implementation.html) following a specified design to gain precision on size factor estimates.

# Usage

Take a look at the [Get started](https://toddmccready.github.io/disize/articles/disize.html) page to familiarize yourself with `disize`.

# Installation

As `disize` is not yet on CRAN, installation is not a one-liner with `install.packages`:

## With `remotes`
```R
# Install disize
remotes::install_github("https://github.com/toddmccready/disize")

# Set up CmdStan toolchain
cmdstanr::install_cmdstan()
```

## With [`rv`](https://a2-ai.github.io/rv-docs/)

Add the following entry to your `rproject.toml` file (if not already present):
```
dependencies = [
    # ...
    { name = "disize", git = "https://github.com/toddmccready/disize", tag = "v0.4.38" },
    # ...
]
```

Then install the CmdStan toolchain in `R`:
```R
cmdstanr::install_cmdstan()
```

# disize

**D**esign **i**nformed **size** factor estimation (or `disize`) is a normalization method meant to be an alternative to `DESeq2`'s [median of rations](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106) and `edgeR`'s [trimmed mean of M values](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)

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
    { name = "disize", git = "https://github.com/toddmccready/disize", tag = "v0.4.18" },
    # ...
]
```

Then sync your project:
```sh
rv sync
```

And finally install the CmdStan toolchain in `R`:
```R
cmdstanr::install_cmdstan()
```

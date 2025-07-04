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
    { name = "disize", git = "https://github.com/toddmccready/disize", tag = "v0.4.20" },
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

# Implementation

Internally, `disize` uses Stan to fit a Bayesian model that jointly estimates the effect of covariates(structured according to `design_formula`) on expression *and* any confounding batch-effects:

$$\begin{aligned}
    \mathbf{y}_g &\sim \text{NegBinom}(\mathbf{\mu}_g, \phi) \\
    \log \mathbf{\mu}_g &= \mathbf{\alpha} + \mathbf{X} \mathbf{\beta}_g + \mathbf{Z}\mathbf{b}_g + \mathbf{o} \\
        &= \mathbf{\alpha} + \mathbf{X} \mathbf{\beta}_g + \mathbf{Z}\mathbf{b}_g + \mathbf{B} \mathbf{s} \\
    \mathbf{\beta}_g &\sim \text{Normal}(\mathbf{0}, \mathbf{G})
\end{aligned}$$

Where $\mathbf{y}_g$ denotes the vector of counts of a gene $g$ for all observations, which is realized from the distribution parameterized by the effect of the covariates ($\mathbf{\alpha} + \mathbf{X} \mathbf{\beta}_g + \mathbf{Z}\mathbf{b}_g$) and any batch-effects ($\mathbf{B} \vec{s}$).

The experimental design is specified by an R formula (`design_formula`) that constructs a "fixed-effects" model matrix $\mathbf{X}$ (without an intercept) and a "random-effects" model matrix $\mathbf{Z}$; this by itself is a regular GLMM.

The confounding batch-effect is essentially an unknown offset $\mathbf{o} = \mathbf{B} \mathbf{s}$, where $\mathbf{B}$ specifies the batch membership for each observation and $\mathbf{s}$ contains the "size factors" which scale the true magnitude of expression.

At its face, however, this model should not be identifiable for most experimental designs (often the batch ID is perfectly collinear with a predictor or interaction between predictors). This identifiability issue is overcome by assuming only a fraction of features are significantly affected by the covariates measured in the experiment; in other words, the estimated coefficients for the model matrices specifying the experimental design are *sparse*.

This assumption is encoded in the final model by placing distinct horseshoe priors on each of the model coefficients (both the fixed- and random-effects). This allows the priors to be learned independently of each other using the large number of features measured in RNAseq experiments.

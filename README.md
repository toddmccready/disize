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
    { name = "disize", git = "https://github.com/toddmccready/disize", tag = "v0.4.21" },
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

Internally, `disize` uses Stan to fit a Bayesian model that jointly estimates the effect of covariates(structured according to `design_formula`) on gene expression *and* any confounding batch-effects:

$$\begin{aligned}
    \mathbf{y}_g &\sim \text{NegBinom}(\mathbf{\mu}_g, \phi) \\
    \log \mathbf{\mu}_g &= \mathbf{\alpha}_g + \mathbf{X} \mathbf{\beta}_g + \mathbf{Z}\mathbf{b}_g + \mathbf{o} \\
        &= \mathbf{\alpha}_g + \mathbf{X} \mathbf{\beta}_g + \mathbf{Z}\mathbf{b}_g + (\mathbf{B} \mathbf{s}) \\
    \mathbf{b}_g &\sim \text{Normal}(\mathbf{0}, \mathbf{G})
\end{aligned}$$

Where $\mathbf{y}_g$ denotes the vector of counts of a gene $g$ for all observations, which is realized from the distribution parameterized by the effect of the covariates ($\mathbf{\alpha}_g + \mathbf{X} \mathbf{\beta}_g + \mathbf{Z}\mathbf{b}_g$) and any batch-effects ($\mathbf{B} \mathbf{s}$).

The experimental design is specified by an R formula (`design_formula`) that constructs a "fixed-effects" model matrix $\mathbf{X}$ (without an intercept) and a "random-effects" model matrix $\mathbf{Z}$; this by itself is a regular GLMM.

The confounding batch-effect is essentially an unknown offset $\mathbf{o} = \mathbf{B} \mathbf{s}$, where $\mathbf{B}$ specifies the batch membership for each observation and $\mathbf{s}$ contains the "size factors" which scale the true magnitude of expression.

At its face, however, this model should not be identifiable for most experimental designs (often the batch ID is perfectly collinear with a predictor or interaction between predictors). This identifiability issue is overcome by constraining $\mathbf{s}$ and assuming only a fraction of features are significantly affected by the covariates measured in the experiment; in other words, the estimated coefficients $\mathbf{\beta}_g, \mathbf{b}_g$ (excluding the intercept) are *sparse* across genes.

This assumption is encoded in the model by placing distinct horseshoe priors on each of the model coefficients (both the fixed- and random-effects). This allows the priors to be learned independently of each other using the large number of features measured in RNAseq experiments.

## Estimation

Since the posterior distribution for $\mathbf{s}$ is heavily concentrated for typical datasets with thousands of features, it is sufficient to quickly provide a point estimate of $\mathbf{s}$ and delegate estimation of the model coefficients to existing tools like `DESeq2` or `edgeR` (only replacing the small normalization step of their workflows).

`disize` uses Stan's [L-BFGS optimization algorithm](https://mc-stan.org/docs/reference-manual/optimization.html) to find the model's *maximum a posteriori* (MAP) for $\mathbf{s}$. We end up doing this in fewer iterations than needed for all parameters to converge by using a heuristic to guess how long the procedure should run for; this is followed up with a diagnostic to ensure the size factors have converged.

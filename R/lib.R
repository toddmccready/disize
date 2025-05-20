#' Modify the design by adding an interaction by genes
#'
#' @param formula this do stuff
#' @param gene_name this also do stuff
modify_design <- function(formula, gene_name) {
    # Extract terms
    terms <- attr(terms(formula), "term.labels")

    # Wrap random effects in brackets
    terms <- gsub("\\| ", "\\| \\(", terms)
    terms[grepl("\\|", terms)] <- paste0(terms[grepl("\\|", terms)], ")")

    # Add interaction
    terms <- paste0("(", terms, ":", gene_name, ")")

    # Reconstruct formula
    formula <- paste(terms, collapse = " + ")

    list(
        formula = formula,
        terms = terms,
        random = grepl("\\|", terms)
    )
}

#' @title Design-infored size factor estimation.
#'
#'
#' @param design_formula The formula describing the experimental design.
#' @param counts A (obsservation x feature) count matrix.
#' @param metadata A dataframe containing sample metadata.
#' @param model_data The model data if already constructed. Must contain the
#'  'batch_name', 'obs_name', 'gene_name' identifiers as columns and any
#'  predictors used in 'design_formula'.
#' @param batch_name The identifier for the batch column in 'metadata',
#'  defaults to "batch".
#' @param obs_name The identifier for the observation column in 'metadata',
#'  defaults to "obs".
#' @param gene_name The identifier for the gene column in 'model_data',
#'  defaults to "gene".
#' @param n_genes The number of genes used during estimation, defaults to 500.
#'  Increasing this value will result in this function taking longer but more
#'  confidence in the size factors.
#' @param n_threads The number of threads to be used during estimation,
#'  defaults to 1. Increasing this value will generally decrease runtime.
#' @param backend How to call Stan, defaults to "rstan".
#'
#' @returns A named numeric vector containing the size factor point estimates.
#'
#' @export
disize <- function(
    design_formula,
    counts = NULL,
    metadata = NULL,
    model_data = NULL,
    batch_name = "batch",
    obs_name = "obs",
    gene_name = "gene",
    n_genes = 500,
    n_threads = NULL,
    verbose = TRUE
) {
    # Check design formula is correct
    if (!is(design_formula, "formula")) {
        stop("'design_formula' should be an R formula")
    } else if (2 < length(design_formula)) {
        stop("'design_formula' should be of the form '~ x + ...'")
    }

    # Check data is inputted correctly
    if (is.null(model_data) && (!is.null(counts) && !is.null(metadata))) {
        # Check for the same number of samples
        if (nrow(metadata) != nrow(counts)) {
            stop(
                "'counts' and 'metadata' should have the same # of ",
                "samples(rows)."
            )
        }

        # Subset genes for model
        n_genes <- min(n_genes, nrow(counts))
        counts <- counts[, order(colMeans(counts != 0), decreasing = TRUE)[
            1:n_genes
        ]]

        # Format counts to include sample-level and gene-level in long format
        counts <- reshape2::melt(
            counts,
            c(sample_name, gene_name),
            value.name = "counts"
        )
        counts[[sample_name]] <- factor(counts[[sample_name]])
        counts[[gene_name]] <- factor(counts[[gene_name]])

        # Merge counts and metadata
        model_data <- merge(counts, metadata, by = sample_name)
    } else if (!is.null(model_data) && (is.null(counts) && is.null(metadata))) {
        # Ensure relevant terms are factors
        if (!is(model_data[[batch_name]], "factor")) {
            model_data[[batch_name]] <- factor(model_data[[batch_name]])
        }
        if (!is(model_data[[gene_name]], "factor")) {
            model_data[[gene_name]] <- factor(model_data[[gene_name]])
        }
    } else {
        stop(
            "either 'counts', 'metadata' can be specified(and 'model_data' ",
            "left NULL) or 'model_data' can be specified(and 'counts', ",
            "'metadata' left NULL)"
        )
    }

    # Save batch names
    batches <- levels(model_data[[batch_name]])

    # Modify the design formula
    design <- modify_design(design_formula, gene_name)
    design_formula <- paste0("design ~ 0 + ", design$formula)

    # Extract number of batches
    n_batches <- length(batches)

    # Construct formula for brms
    formula <- paste0(
        "counts ~ intercept + design + log(sf) + log(",
        n_batches,
        ")"
    )
    formula <- brms::bf(formula, nl = TRUE) +
        brms::lf(design_formula) +
        brms::lf(paste0("sf ~ 0 + ", batch_name)) +
        brms::lf(paste0("intercept ~ 0 + ", gene_name))

    # Construct priors
    priors <- brms::prior_string(
        paste0("dirichlet(rep_vector(1, ", n_batches, "))"),
        nlpar = "sf"
    )

    # Include optional priors
    if (any(!design$random)) {
        # Induce sparsity in fixed effects
        priors <- priors +
            brms::prior(
                "horseshoe(main = TRUE)",
                class = "b",
                nlpar = "design"
            )
    }
    if (any(design$random)) {
        # Induce sparsity in random effects
        priors <- priors +
            brms::prior("horseshoe(1)", class = "sd", nlpar = "design")
    }

    # Construct Stan code
    if (verbose) message("Compiling Stan model...")
    my_code <- brms::stancode(
        formula,
        data = model_data,
        family = "negbinomial",
        prior = priors,
        threads = brms::threading(n_threads)
    )
    my_data <- brms::standata(
        formula,
        data = model_data,
        family = "negbinomial",
        prior = priors,
        threads = brms::threading(n_threads)
    )

    # Construct model
    model <- rstan::stan_model(model_code = my_code)

    # Estimate model parameters
    if (verbose) message("Estimating size factors...")
    fit <- rstan::optimizing(model, data = my_data)

    # Extract and name size factors
    size_factors <- fit$theta_tilde[grepl("sf", colnames(fit$theta_tilde))]
    names(size_factors) <- batches

    # Transform size factors to more useful scale
    size_factors <- log(size_factors) + log(n_batches)

    size_factors
}

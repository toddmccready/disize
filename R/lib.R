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

#' @title Estimate size factors
#'
#' @param design_formula The formula describing the experimental design.
#' @param counts A samples x genes count matrix.
#' @param metadata A dataframe containing sample metadata.
#' @param model_data The model data if already formatted.
#' @param batch_name The identifier for the batch column in 'metadata'.
#' @param sample_name The identifier for the sample column in 'metadata'.
#' @param gene_name The identifier for the gene column in 'model_data'.
#' @param n_genes The number of genes used to estimate size factors.
#'
#' @export
disize <- function(
    design_formula,
    counts = NULL,
    metadata = NULL,
    model_data = NULL,
    batch_name = "batch",
    sample_name = "sample",
    gene_name = "gene",
    n_genes = 500
) {
    # Check design formula is correct
    if (!is(design_formula, "formula")) {
        stop("'design_formula' should be an R formula.")
    } else if (2 < length(design_formula)) {
        stop("'design_formula' should be of the form '~ x + ...'.")
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
        model_data <- model_data
    } else {
        stop(
            "either 'counts', 'metadata' can be specified(and 'model_data' ",
            "left NULL) or 'model_data' can be specified."
        )
    }

    # Modify the design formula
    design <- modify_design(design_formula)
    design_formula <- paste0("design ~ 0 + ", design$formula)

    # Extract number of batches
    n_batches <- length(unique(model_data[[batch_name]]))

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
    if (any(design$random)) {
        # Induce sparsity in random effects
        priors <- priors +
            brms::prior("horseshoe(1)", class = "sd", nlpar = "design")
    }
    if (any(!design$random)) {
        # Induce sparsity in fixed effects
        priors <- priors +
            brms::prior("horseshoe(1)", class = "b", nlpar = "design")
    }

    # Estimate model parameters
    model <- brms::brm(
        formula,
        model_data,
        family = stats::poisson(),
        priors,
        algorithm = "meanfield",
        iter = 1e5
    )

    # Extract size factors
    size_factors <- log(brms::fixef(model)[1:n_batches, 1]) + log(n_batches)

    size_factors
}

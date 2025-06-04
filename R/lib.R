#' Modify the design by adding an interaction by genes
#'
#' @param formula The design formula
modify_design <- function(formula) {
    # Extract terms
    terms <- attr(terms(formula), "term.labels")

    # Wrap random effects in brackets
    terms <- gsub("\\| ", "\\| \\(", terms)
    terms[grepl("\\|", terms)] <- paste0(terms[grepl("\\|", terms)], ")")

    # Add interaction
    terms <- paste0("(", terms, ":feat_idx)")

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
#'  'batch_name' and "feat_idx" identifiers as columns, and any
#'  predictors used in 'design_formula'.
#' @param batch_name The identifier for the batch column in 'metadata',
#'  defaults to "batch".
#' @param n_feats The number of genes used during estimation, defaults to 500.
#'  Increasing this value will result in this function taking longer but more
#'  confidence in the size factors.
#' @param n_subset The number of observations per experimental unit used during
#'  estimation, defaults to 50.
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
    n_feats = 500,
    n_subset = 50,
    n_threads = NULL,
    verbose = 3,
    n_passes = 20,
    n_iters = 100,
    tolerance = 1e-3
) {
    # Argument checks ----
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
                "observations(rows)."
            )
        }

        # Remove dimnames in 'counts' if present
        dimnames(data) <- list(NULL, NULL)

        # Include explicit observation index variable
        metadata[["obs_idx"]] <- 1:nrow(metadata)

        # Ensure valid number of features selected
        n_feats <- min(n_feats, ncol(counts))

        # Subset features
        ordering <- order(Matrix::colMeans(counts), decreasing = TRUE)
        counts <- counts[, ordering[1:n_feats]]

        # Include original feature indices
        colnames(counts) <- ordering[1:n_feats]

        # Subset observations
        predictors <- all.vars(design_formula)
        metadata <- metadata |>
            dplyr::group_by(dplyr::across(dplyr::all_of(predictors))) |>
            dplyr::slice_sample(n = n_subset, replace = FALSE) |>
            dplyr::ungroup() |>
            dplyr::select(dplyr::all_of(c(predictors, obs_name, batch_name)))

        counts <- counts[metadata[[obs_name]], ]

        # Convert to dense matrix if needed
        if (is(counts, "sparseMatrix")) {
            counts <- as.matrix(counts)
        }

        # Format counts to long format
        counts <- reshape2::melt(
            counts,
            c("obs_idx", "feat_idx"),
            value.name = "counts"
        )

        # Merge counts and metadata
        model_data <- merge(counts, metadata, by = obs_name)
    } else if (!is.null(model_data) && (is.null(counts) && is.null(metadata))) {
    } else {
        stop(
            "either 'counts', 'metadata' can be specified(and 'model_data' ",
            "left NULL) or 'model_data' can be specified(and 'counts', ",
            "'metadata' left NULL)"
        )
    }

    # Formatting ----
    # Ensure relevant columns are factors
    if (!is(model_data[[batch_name]], "factor")) {
        model_data[[batch_name]] <- factor(model_data[[batch_name]])
    }
    if (!is(model_data[["feat_idx"]], "factor")) {
        model_data[["feat_idx"]] <- factor(model_data[["feat_idx"]])
    }

    # Format characters to factors
    model_data <- model_data |>
        dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

    # Save batch names
    batches <- levels(model_data[[batch_name]])

    # Modify the design formula
    design <- modify_design(design_formula)
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
        brms::lf(paste0("intercept ~ 0 + feat_idx"))

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
                "horseshoe(main = TRUE, scale_slab = 10000)",
                class = "b",
                nlpar = "design"
            )
    }
    if (any(design$random)) {
        # Induce sparsity in random effects
        priors <- priors +
            brms::prior(
                "horseshoe(1, scale_slab = 10000)",
                class = "sd",
                nlpar = "design"
            )
    }

    # Construct Stan code
    if (verbose) message("Compiling Stan model...")
    model_code <- brms::stancode(
        formula,
        data = model_data,
        family = "negbinomial",
        prior = priors,
        threads = brms::threading(n_threads)
    )
    model_data <- brms::standata(
        formula,
        data = model_data,
        family = "negbinomial",
        prior = priors,
        threads = brms::threading(n_threads)
    )

    # Compile model
    model <- rstan::stan_model(
        model_code = model_code
    )

    # Estimate model parameters ----
    if (verbose) message("Estimating size factors...")

    # Free up space
    gc()

    # Construct progress bar
    if (2 < verbose) {
        pb <- progress::progress_bar$new(total = n_passes)
        pb$tick(0)
    }

    # Compute initial fit
    cur_fit <- rstan::optimizing(
        model,
        data = model_data,
        iter = n_iters,
        as_vector = FALSE
    )

    # Extract size factors
    sf_hist <- list()
    sf_hist[[1]] <- log(cur_fit$par$b_sf) + log(n_batches)

    if (2 < verbose) pb$tick()
    for (i in 2:n_passes) {
        # Compute next fit
        cur_fit <- rstan::optimizing(
            model,
            data = model_data,
            iter = n_iters,
            init = cur_fit$par,
            as_vector = FALSE
        )

        # Extract size factors
        sf_hist[[i]] <- log(cur_fit$par$b_sf) + log(n_batches)

        # Free up space
        gc()

        # Evaluate convergence
        # TODO: do something smart with the history
        if (all(abs(sf_hist[[i]] - sf_hist[[i - 1]]) < tolerance)) {
            # Name and return size factors
            sf <- sf_hist[[i]]
            names(sf) <- batches

            if (2 < verbose) pb$terminate()
            return(sf)
        }

        if (2 < verbose) pb$tick()
    }

    if (1 < verbose) {
        warning("Model did not converge, size factors may be imprecise.")
    }

    # Terminate progress bar
    if (2 < verbose) pb$terminate()

    # Name and return size factors
    sf <- sf_hist[[i]]
    names(sf) <- batches

    sf
}

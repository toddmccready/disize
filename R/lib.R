#' Modify the design by adding an interaction by genes
#'
#' @param design_formula The design formula
modify_design <- function(design_formula, feat_name) {
    # Extract terms
    terms <- attr(terms(design_formula), "term.labels")

    # Identify random effects
    re <- grepl("\\| ", terms)

    # Wrap random effects in brackets
    terms <- gsub("\\| ", "\\| \\(", terms)
    terms[re] <- paste0(terms[re], ")")

    # Add interaction
    terms <- paste0("(", terms, ":", feat_name, ")")

    fixed <- NULL
    if (length(terms[!re]) != 0) {
        fixed <- formula(paste0(" ~ 0 + ", paste(terms[!re], collapse = " + ")))
    }

    random <- NULL
    if (length(terms[re]) != 0) {
        random <- formula(paste0(" ~ 0 + ", paste(terms[re], collapse = " + ")))
    }

    list(
        formula = as.formula(paste0(" ~ 0 + ", paste(terms, collapse = " + "))),
        fixed = fixed,
        random = random
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
    obs_name = "obs_id",
    feat_name = "feat_id",
    n_feats = 500,
    n_subset = 50,
    n_threads = NULL,
    verbose = 3,
    n_passes = 20,
    n_iters = 100,
    tolerance = 1e-3
) {
    # Check design formula is correct
    if (!is(design_formula, "formula")) {
        stop("'design_formula' should be an R formula")
    } else if (2 < length(design_formula)) {
        stop("'design_formula' should be of the form '~ x + ...'")
    }

    if (is.null(model_data) && (!is.null(counts) && !is.null(metadata))) {
        # Check for the same number of samples
        if (nrow(metadata) != nrow(counts)) {
            stop(
                "'counts' and 'metadata' should have the same # of ",
                "observations(rows)."
            )
        }

        # Include explicit observation names if not present
        if (is.null(rownames(counts))) {
            rownames(counts) <- 1:nrow(counts)
        }
        if (is.null(metadata[[obs_name]])) {
            metadata[[obs_name]] <- 1:nrow(counts)
        }

        # Include explicit feature names if not present
        if (is.null(colnames(counts))) {
            colnames(counts) <- 1:ncol(counts)
        }

        # Ensure valid number of features selected
        n_feats <- min(n_feats, ncol(counts))

        # Subset features
        ordering <- order(Matrix::colMeans(counts), decreasing = TRUE)
        counts <- counts[, ordering[1:n_feats]]

        # Extract predictor terms
        predictors <- all.vars(design_formula)

        # Subset observations
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
            c(obs_name, feat_name),
            value.name = "count"
        )

        # Merge counts and metadata
        model_data <- merge(counts, metadata, by = obs_name)
    } else if (!is.null(model_data) && (is.null(counts) && is.null(metadata))) {
        # TODO make n_subset work here
    } else {
        stop(
            "either 'counts', 'metadata' can be specified(and 'model_data' ",
            "left NULL) or 'model_data' can be specified(and 'counts', ",
            "'metadata' left NULL)"
        )
    }

    # Ensure relevant columns are factors
    model_data[[batch_name]] <- as.factor(model_data[[batch_name]])
    model_data[[feat_name]] <- as.factor(model_data[[feat_name]])

    # Format characters to factors
    model_data <- model_data |>
        dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

    # Save batch names and number
    batches <- levels(model_data[[batch_name]])
    n_batches <- length(batches)

    # Modify the design formula
    design <- modify_design(design_formula, feat_name)

    # Construct model matrix
    if (!is.null(design$fixed)) {
        fixed_matrix <- Matrix::sparse.model.matrix(design$fixed, model_data)
    }

    # Construct random effects matrix
    if (!is.null(design$random)) {
        random_matrices <- lapply(
            X = reformulas::mkReTrms(
                bars = reformulas::findbars(design$random),
                fr = model_data,
                calc.lambdat = FALSE
            )$Ztlist,
            FUN = function(Z) {
                Matrix::t(Z)
            }
        )
    }

    # TODO: Extract matrix data in CSR format

    # TODO: See whether we need multiple Stan models (f, r, f + r) or just one

    # Construct Stan model
    if (verbose) message("Compiling Stan model...")

    # Compile model

    # Estimate model parameters ----
    if (verbose) message("Estimating size factors...")

    # Construct progress bar
    if (2 < verbose) {
        pb <- progress::progress_bar$new(total = n_passes)
        pb$tick(0)
    }

    # Extract size factors
    sf_hist <- list()

    if (2 < verbose) pb$tick()
    for (i in 2:n_passes) {
        # Compute next fit
        cur_fit <- rstan::optimizing(
            model,
            data = stan_data,
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

#' Split design formula into fixed and random effects.
#'
#' @param design_formula The design formula
split_formula <- function(design_formula) {
    # Extract terms
    terms <- attr(terms(design_formula), "term.labels")

    # Identify random effects
    re <- grepl("\\| ", terms)

    fixed <- NULL
    if (length(terms[!re]) != 0) {
        fixed <- formula(paste0(" ~ 0 + ", paste(terms[!re], collapse = " + ")))
    }

    random <- NULL
    if (length(terms[re]) != 0) {
        random <- formula(paste0(" ~ 0 + ", paste(terms[re], collapse = " + ")))
    }

    list(
        formula = design_formula,
        fixed = fixed,
        random = random
    )
}

#' Extract the parameters out of an optimized model.
#'
#' @param cur_fit A CmdStanMLE object
#' @param stan_data The formatted Stan data used to fit the model.
extract_sf <- function(cur_fit) {
    params <- cur_fit$mle()

    unname(params[grepl("^sf", names(params))])
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
#'  defaults to "batch_id".
#' @param obs_name The identifier for the observation column in 'metadata',
#'  defaults to "obs_id".
#' @param feat_name The identifier for the feature column in 'metadata',
#'  defaults to "feat_id".
#' @param n_feats The number of genes used during estimation, defaults to 500.
#'  Increasing this value will result in this function taking longer but more
#'  confidence in the size factors.
#' @param n_subset The number of observations per experimental unit used during
#'  estimation, defaults to 50.
#' @param n_threads The number of threads to use for parallel processing.
#' @param verbose The verbosity level.
#' @param n_passes The number of optimization passes to go through.
#' @param n_iters The number of iterations used for a single optimization pass.
#' @param tolerance The tolerance used to evaluate convergence of the size factors.
#'
#' @returns A named numeric vector containing the size factor point estimates.
#'
#' @export
disize <- function(
    design_formula,
    counts = NULL,
    metadata = NULL,
    batch_name = "batch_id",
    obs_name = "obs_id",
    feat_name = "feat_id",
    n_feats = 1000,
    n_subset = 50,
    n_threads = 1,
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

    # Check for the same number of samples
    if (nrow(metadata) != nrow(counts)) {
        stop(
            "'counts' and 'metadata' should have the same # of ",
            "observations(rows)."
        )
    }

    # Include explicit observation names if not present
    if (is.null(rownames(counts)) & is.null(metadata[[obs_name]])) {
        rownames(counts) <- 1:nrow(counts)
        metadata[[obs_name]] <- 1:nrow(counts)
    } else if (!is.null(rownames(counts)) & is.null(metadata[[obs_name]])) {
        metadata[[obs_name]] <- rownames(counts)
    } else if (is.null(rownames(counts)) & !is.null(metadata[[obs_name]])) {
        rownames(counts) <- metadata[[obs_name]]
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

    # Ensure relevant columns are factors
    metadata[[batch_name]] <- as.factor(metadata[[batch_name]])

    # Formating Data For Stan ----
    if (2 < verbose) {
        message("Formatting data...")
    }

    # Allocate named list for Stan
    stan_data <- list(n_obs = nrow(counts), n_feats = ncol(counts))

    # Include batch-level for Stan
    stan_data[["n_batches"]] <- length(levels(metadata[[batch_name]]))
    stan_data[["batch_id"]] <- as.integer(metadata[[batch_name]])

    # Modify the design formula
    design <- split_formula(design_formula)

    # Construct fixed-effects model matrix if present
    if (!is.null(design$fixed)) {
        fe_design <- model.matrix(design$fixed, metadata)

        stan_data[["n_fe"]] <- ncol(fe_design)
        stan_data[["fe_design"]] <- fe_design
    } else {
        stan_data[["n_fe"]] <- 0L
        stan_data[["fe_design"]] <- array(0, dim = c(stan_data$n_obs, 0))
    }

    # Construct random-effects matrix if present
    if (!is.null(design$random)) {
        remm <- reformulas::mkReTrms(
            bars = reformulas::findbars(design$random),
            fr = metadata,
            calc.lambdat = FALSE,
            sparse = TRUE
        )
        re_design <- Matrix::t(remm$Zt) |> as("RsparseMatrix")

        # Check if all random-effects terms are scalar normals
        all_scalar <- lapply(remm$cnms, function(b) {
            length(b) == 1
        }) |>
            unlist() |>
            all()
        if (!all_scalar) {
            stop(
                "Only include one predictor on the LHS of a random-effects bar."
            )
        }

        # Include data for Stan
        stan_data[["n_re"]] <- ncol(re_design)
        stan_data[["n_nz_re"]] <- length(re_design@x)
        stan_data[["re_design_x"]] <- re_design@x
        stan_data[["re_design_j"]] <- re_design@j + 1L
        stan_data[["re_design_p"]] <- re_design@p + 1L
        stan_data[["n_re_terms"]] <- length(remm$cnms)
        stan_data[["re_id"]] <- rep(1:length(remm$cnms), times = diff(remm$Gp))
    } else {
        stan_data[["n_re"]] <- 0L
        stan_data[["n_nz_re"]] <- 0L
        stan_data[["re_design_x"]] <- numeric(0)
        stan_data[["re_design_j"]] <- integer(0)
        stan_data[["re_design_p"]] <- integer(stan_data[["n_obs"]] + 1)
        stan_data[["n_re_terms"]] <- 0L
        stan_data[["re_id"]] <- integer(0)
    }

    # Include counts for Stan
    stan_data[["counts"]] <- counts |> t()

    # Construct Stan model
    model <- instantiate::stan_package_model(
        name = "disize",
        package = "disize",
        compile = TRUE,
        force = TRUE,
        cpp_options = list(stan_threads = TRUE)
    )

    # Estimate model parameters ----
    if (verbose) {
        message("Estimating size factors...")
    }

    # Construct progress bar
    if (2 < verbose) {
        pb <- progress::progress_bar$new(total = n_passes)
        pb$tick(0)
    }

    # Estimate initial fit
    options(cmdstanr_warn_inits = FALSE)
    cur_fit <- model$optimize(
        data = stan_data,
        iter = n_iters,
        show_messages = FALSE,
        sig_figs = 18,
        threads = n_threads
    )

    # Extract size factors
    sf_hist <- list()
    sf_hist[[1]] <- extract_sf(cur_fit)

    if (2 < verbose) {
        pb$tick()
    }
    for (i in 2:n_passes) {
        # Compute next fit
        cur_fit <- model$optimize(
            data = stan_data,
            init = cur_fit,
            iter = n_iters,
            show_messages = FALSE,
            sig_figs = 18,
            threads = n_threads
        )

        # Extract size factors
        sf_hist[[i]] <- extract_sf(cur_fit)

        # Evaluate convergence
        # TODO: do something smart with the history
        if (all(abs(sf_hist[[i]] - sf_hist[[i - 1]]) < tolerance)) {
            # Name and return size factors
            sf <- as.vector(sf_hist[[i]])
            names(sf) <- levels(metadata[[batch_name]])

            if (2 < verbose) {
                pb$terminate()
            }
            return(sf)
        }

        if (2 < verbose) pb$tick()
    }

    if (1 < verbose) {
        warning("Model did not converge, size factors may be imprecise.")
    }

    # Terminate progress bar
    if (2 < verbose) {
        pb$terminate()
    }

    # Name and return size factors
    sf <- as.vector(sf_hist[[i]])
    names(sf) <- levels(metadata[[batch_name]])

    sf
}

#' Split design formula into fixed and random effects.
#'
#' @param design_formula The design formula
#'
#' @noRd
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

#' @title Design-informed size factor estimation
#'
#' @param design_formula The formula describing the experimental design.
#' @param counts A (obsservation x feature) count matrix.
#' @param metadata A dataframe containing sample metadata.
#' @param batch_name The identifier for the batch column in 'metadata',
#'  defaults to "batch_id".
#' @param obs_name The identifier for the observation column in 'metadata',
#'  defaults to "obs_id".
#' @param n_feats The number of features used during estimation, defaults to
#'  automatically choosing based on a heuristic.
#' @param n_subset The number of observations per experimental unit used during
#'  estimation, defaults to 50.
#' @param n_iters The number of iterations used for a single optimization pass.
#' @param n_threads The number of threads to use for parallel processing.
#' @param verbose The verbosity level.
#'
#' @returns A named numeric vector containing the size factor point estimates.
#'
#' @export
disize <- function(
    design_formula,
    counts,
    metadata,
    batch_name = "batch_id",
    obs_name = "obs_id",
    n_feats = min(10000L, ncol(counts)),
    n_subset = 50L,
    n_iters = "auto",
    n_threads = ceiling(parallel::detectCores() / 2),
    grainsize = ceiling(n_feats / n_threads),
    init_alpha = 1e-6,
    verbose = 3L
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
    if (is.null(rownames(counts)) && is.null(metadata[[obs_name]])) {
        rownames(counts) <- 1:nrow(counts)
        metadata[[obs_name]] <- 1:nrow(counts)
    } else if (!is.null(rownames(counts)) && is.null(metadata[[obs_name]])) {
        metadata[[obs_name]] <- rownames(counts)
    } else if (is.null(rownames(counts)) && !is.null(metadata[[obs_name]])) {
        rownames(counts) <- metadata[[obs_name]]
    }

    # Include explicit feature names if not present
    if (is.null(colnames(counts))) {
        colnames(counts) <- 1:ncol(counts)
    }

    # Extract predictor terms
    predictors <- all.vars(design_formula)

    # Subset observations
    metadata <- metadata |>
        dplyr::group_by(dplyr::across(dplyr::all_of(c(
            predictors,
            batch_name
        )))) |>
        dplyr::slice_sample(n = n_subset, replace = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(c(predictors, obs_name, batch_name)))

    counts <- counts[metadata[[obs_name]], ]

    # Subset features out features with no counts
    subset <- Matrix::colSums(counts) != 0
    counts <- counts[, subset]

    # Ensure valid number of features selected
    if (1 < verbose && ncol(counts) < n_feats) {
        warning(
            "Insufficient number of features after subsetting observations (",
            ncol(counts),
            ") to satisfy n_feats = ",
            n_feats,
            "."
        )
    }
    n_feats <- min(n_feats, ncol(counts))

    # Preferentially subset features with sufficient counts
    ordering <- order(Matrix::colMeans(counts), decreasing = TRUE)
    counts <- counts[, ordering[1:n_feats]]

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

    # Include grainsize
    stan_data[["grainsize"]] <- grainsize

    # Compute heuristic for maximum # of iterations
    if (n_iters == "auto") {
        n_iters <- as.integer(
            100 *
                sqrt(
                    stan_data[["n_batches"]] +
                        stan_data[["n_fe"]] +
                        stan_data[["n_re"]]
                ) +
                1000 * log10(stan_data[["n_feats"]])
        )
    }

    # Construct Stan model
    model <- instantiate::stan_package_model(
        name = "disize",
        package = "disize",
        compile = TRUE,
        force = TRUE,
        cpp_options = list(stan_threads = TRUE),
        compile_model_methods = TRUE
    )

    # Estimate maximum time-to-fit
    fit <- model$optimize(
        data = stan_data,
        iter = 50L,
        threads = n_threads,
        algorithm = "lbfgs",
        init_alpha = init_alpha,
        history_size = 8L,
        sig_figs = 16L,
        show_messages = FALSE
    )

    # Estimate model parameters ----
    if (3L <= verbose) {
        message(
            "Estimating size factors... (Max ETA: ~",
            round(n_iters * fit$time()$total / 50, 1),
            "s)"
        )
    } else if (4L <= verbose) {}

    # Estimate fit
    fit <- model$optimize(
        data = stan_data,
        iter = n_iters,
        threads = n_threads,
        algorithm = "lbfgs",
        init_alpha = init_alpha,
        history_size = 8L,
        sig_figs = 16L,
        show_messages = (4L <= verbose),
        refresh = ceiling(n_iters / 10)
    )

    # Extract parameter estimates
    params <- fit$mle()

    # Remove transformed size factors
    params <- params[!grepl("^sf", names(params))]

    # Extract variable skeleton
    skeleton <- fit$variable_skeleton(FALSE)

    # Extract unconstrained variable
    uc_params <- fit$unconstrain_variables(relist(params, skeleton))

    # Compute gradient
    gradient <- fit$grad_log_prob(uc_params)
    names(gradient) <- names(params)[
        names(params) != paste0("raw_sf[", stan_data[["n_batches"]], "]")
    ]

    # Extract size factor terms
    end <- length(uc_params) -
        (length(skeleton$fe_tau) +
            length(skeleton$re_tau) +
            length(skeleton$lambda) +
            1)
    start <- end - (length(skeleton$raw_sf) - 1) + 1
    sf_gradient <- gradient[start:end]

    # Calculate norm
    norm <- sum(sf_gradient^2) /
        max(abs(attr(gradient, "log_prob")), 1.0)

    # Check gradient for size factors
    if (1 < norm && 2L <= verbose) {
        warning(
            "Relative size factor gradient (",
            norm,
            ") not sufficiently small, estimates may not have converged. Try ",
            "increasing 'n_iters' from ",
            "current value:.",
            n_iters,
            "."
        )
    } else if (4L <= verbose) {
        message("Relative size factor gradient: ", norm)
    }

    # Extract size factors
    sf <- fit$mle("sf")
    names(sf) <- levels(metadata[[batch_name]])

    sf
}

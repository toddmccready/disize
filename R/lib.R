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
        intercept = as.formula(paste0("~ 0 + ", feat_name)),
        fixed = fixed,
        random = random
    )
}

extract_params <- function(cur_fit) {
    params <- cur_fit$mle()

    params_list <- list(
        tau = unname(params["tau"]),
        iodisp = unname(params["iodisp"]),
        raw_sf = unname(params[grepl("^raw_sf", names(params))]),
        int_coefs = unname(params[grepl("^int_coefs", names(params))]),
        sf = unname(params[grepl("^sf", names(params))])
    )

    if (any(grepl("fe_", names(params)))) {
        params_list[["lambda"]] <- unname(params[grepl("^lambda", names(params))])
        params_list[["fe_coefs"]] <- unname(params[grepl("^fe_coefs", names(params))])
    } else {
        params_list[["lambda"]] <- numeric(0)
        params_list[["fe_coefs"]] <- numeric(0)
    }

    if (any(grepl("re_", names(params)))) {
        params_list[["z_re"]] <- unname(params[grepl("^z_re", names(params))])
        params_list[["re_coefs"]] <- unname(params[grepl("^re_coefs", names(params))])
        params_list[["re_sigma"]] <- unname(params[grepl("^re_sigma", names(params))])
    } else {
        params_list[["z_re"]] <- numeric(0)
        params_list[["re_coefs"]] <- numeric(0)
        params_list[["re_sigma"]] <- numeric(0)
    }

    params_list
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
    batch_name = "batch_id",
    obs_name = "obs_id",
    feat_name = "feat_id",
    n_feats = 500,
    n_subset = 50,
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
            value.name = "counts"
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
    if (2 < verbose) message("Formatting data...")

    # Allocate list for formatted Stan data
    stan_data <- list(n_obs = nrow(model_data))

    # Ensure relevant columns are factors
    model_data[[batch_name]] <- as.factor(model_data[[batch_name]])
    model_data[[feat_name]] <- as.factor(model_data[[feat_name]])

    # Include number of features
    stan_data[["n_feats"]] <- length(levels(model_data[[feat_name]]))
    stan_data[["feat_id"]] <- as.integer(model_data[[feat_name]])

    # Format characters to factors
    model_data <- model_data |>
        dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

    # Include batch-level for Stan
    stan_data[["n_batches"]] <- levels(model_data[[batch_name]]) |>
        length()
    stan_data[["batch_id"]] <- as.integer(model_data[[batch_name]])

    # Modify the design formula
    design <- modify_design(design_formula, feat_name)

    # Construct intercept model matrix
    int_design <- Matrix::sparse.model.matrix(design$intercept, model_data) |>
        as("RsparseMatrix")

    stan_data[["n_int"]] <- ncol(int_design)
    stan_data[["n_nz_int"]] <- length(int_design@x)
    stan_data[["int_design_x"]] <- int_design@x
    stan_data[["int_design_j"]] <- int_design@j + 1L
    stan_data[["int_design_p"]] <- int_design@p + 1L

    # Construct fixed-effects model matrix if present
    if (!is.null(design$fixed)) {
        fe_design <- Matrix::sparse.model.matrix(design$fixed, model_data) |>
            as("RsparseMatrix")

        stan_data[["n_fe"]] <- ncol(fe_design)
        stan_data[["n_nz_fe"]] <- length(fe_design@x)
        stan_data[["fe_design_x"]] <- fe_design@x
        stan_data[["fe_design_j"]] <- fe_design@j + 1L
        stan_data[["fe_design_p"]] <- fe_design@p + 1L
    } else {
        stan_data[["n_fe"]] <- 0L
        stan_data[["n_nz_fe"]] <- 0L
        stan_data[["fe_design_x"]] <- numeric(0)
        stan_data[["fe_design_j"]] <- integer(0)
        stan_data[["fe_design_p"]] <- integer(stan_data[["n_obs"]] + 1)
    }

    # Construct random effects matrix
    if (!is.null(design$random)) {
        remm <- reformulas::mkReTrms(
            bars = reformulas::findbars(design$random),
            fr = model_data,
            calc.lambdat = FALSE,
            sparse = TRUE
        )
        re_design <- Matrix::t(remm$Zt) |> as("RsparseMatrix")

        # Check if all random-effects terms are scalar normals
        all_scalar <- lapply(remm$cnms, function(b) {length(b) == 1}) |>
            unlist() |>
            all()
        if (!all_scalar) {
            stop("Only include one predictor on the LHS of a random-effects pipe.")
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
    stan_data[["counts"]] <- as.integer(model_data[["counts"]])

    # Construct Stan model
    model <- instantiate::stan_package_model(name = "disize", package = "disize")

    # Estimate model parameters ----
    if (verbose) message("Estimating size factors...")

    # Construct progress bar
    if (2 < verbose) {
        pb <- progress::progress_bar$new(total = n_passes)
        pb$tick(0)
    }

    # Estimate initial fit
    cur_fit <- model$optimize(stan_data, iter = n_iters, show_messages = F, sig_figs = 18)

    # Extract parameters
    cur_params <- extract_params(cur_fit)

    # Extract size factors
    sf_hist <- list()
    sf_hist[[1]] <- cur_params[["sf"]]

    if (2 < verbose) pb$tick()
    for (i in 2:n_passes) {
        # Compute next fit
        cur_fit <- model$optimize(stan_data, init = list(cur_params), iter = n_iters, show_messages = F, sig_figs = 18)

        # Extract parameters
        cur_params <- extract_params(cur_fit)

        # Extract size factors
        sf_hist[[i]] <- cur_params[["sf"]]

        # Evaluate convergence
        # TODO: do something smart with the history
        if (all(abs(sf_hist[[i]] - sf_hist[[i - 1]]) < tolerance)) {
            # Name and return size factors
            sf <- as.vector(sf_hist[[i]])
            names(sf) <- levels(model_data[[batch_name]])

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
    sf <- as.vector(sf_hist[[i]])
    names(sf) <- levels(model_data[[batch_name]])

    sf
}

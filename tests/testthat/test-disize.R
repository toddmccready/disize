library(tibble)
library(dplyr)
library(tidyr)
library(purrr)

test_that("small-simple-bulk", {
    # Simulate data
    n_g <- 1000
    n_d <- 2
    n_o <- 1

    data <- tibble(feat_id = 1:n_g) |>
        # Simulate true expression quantities
        mutate(
            d = pmap(list(feat_id), function(g) {
                nonzero <- rbinom(1, 1, 0.1) == 1

                if (nonzero) {
                    return(tibble(
                        donor = 1:n_d,
                        q = rlnorm(n_d, meanlog = rnorm(1), sdlog = 0.5)
                    ))
                } else {
                    return(tibble(
                        donor = 1:n_d,
                        q = rlnorm(n_d, meanlog = rnorm(1), sdlog = 0.0)
                    ))
                }
            })
        ) |>
        unnest(d)

    # Simulate batch-effect
    true_sf <- tibble(
        donor = 1:n_d,
        batch_id = 1:n_d,
        sf = runif(n_d, 1.0, 10.0)
    )

    data <- data |>
        left_join(
            true_sf,
            by = "donor"
        ) |>
        # Generate counts
        mutate(
            d = pmap(list(feat_id, donor, q, sf), function(g, d, q, sf) {
                tibble(
                    obs_id = paste0(d, "_", 1:n_o),
                    counts = as.integer(rnbinom(n_o, mu = q * sf, size = 100))
                )
            })
        ) |>
        unnest(d) |>
        select(feat_id, donor, batch_id, obs_id, counts)

    # Coerce relevant columns to a factor
    data[["donor"]] <- factor(data[["donor"]])

    # Format into 'counts' and 'metadata'
    counts <- data |>
        select(feat_id, obs_id, counts) |>
        pivot_wider(names_from = obs_id, values_from = counts) |>
        column_to_rownames("feat_id") |>
        as.matrix() |>
        Matrix::t()

    metadata <- data %>%
        select(donor, batch_id, obs_id) %>%
        distinct()

    # Compute size factors
    size_factors <- exp(disize::disize(
        design_formula = ~ (1 | donor),
        counts = counts,
        metadata = metadata
    ))

    expect_equal(
        object = as.vector(unname(size_factors)),
        expected = (true_sf$sf / sum(true_sf$sf)) * n_d,
        tolerance = 0.1
    )
})

test_that("large-simple-bulk", {
    # Simulate data
    n_g <- 1000
    n_d <- 12
    n_o <- 1

    data <- tibble(feat_id = 1:n_g) |>
        # Simulate true expression quantities
        mutate(
            d = pmap(list(feat_id), function(g) {
                nonzero <- rbinom(1, 1, 0.1) == 1

                if (nonzero) {
                    return(tibble(
                        donor = 1:n_d,
                        q = rlnorm(n_d, meanlog = rnorm(1), sdlog = 0.5)
                    ))
                } else {
                    return(tibble(
                        donor = 1:n_d,
                        q = rlnorm(n_d, meanlog = rnorm(1), sdlog = 0.0)
                    ))
                }
            })
        ) |>
        unnest(d)

    # Simulate batch-effect
    true_sf <- tibble(
        donor = 1:n_d,
        batch_id = 1:n_d,
        sf = runif(n_d, 1.0, 10.0)
    )

    data <- data |>
        left_join(
            true_sf,
            by = "donor"
        ) |>
        # Generate counts
        mutate(
            d = pmap(list(feat_id, donor, q, sf), function(g, d, q, sf) {
                tibble(
                    obs_id = paste0(d, "_", 1:n_o),
                    counts = as.integer(rnbinom(n_o, mu = q * sf, size = 100))
                )
            })
        ) |>
        unnest(d) |>
        select(feat_id, donor, batch_id, obs_id, counts)

    # Coerce relevant columns to a factor
    data[["donor"]] <- factor(data[["donor"]])

    # Format into 'counts' and 'metadata'
    counts <- data |>
        select(feat_id, obs_id, counts) |>
        pivot_wider(names_from = obs_id, values_from = counts) |>
        column_to_rownames("feat_id") |>
        as.matrix() |>
        Matrix::t()

    metadata <- data %>%
        select(donor, batch_id, obs_id) %>%
        distinct()

    # Compute size factors
    size_factors <- exp(disize::disize(
        design_formula = ~ (1 | donor),
        counts = counts,
        metadata = metadata
    ))

    expect_equal(
        object = as.vector(unname(size_factors)),
        expected = (true_sf$sf / sum(true_sf$sf)) * n_d,
        tolerance = 0.1
    )
})

test_that("small-simple-sc", {
    # Simulate data
    n_g <- 1000
    n_d <- 4
    n_p <- 3
    n_o <- 50

    data <- tibble(feat_id = 1:n_g) |>
        crossing(tibble(cell_type = 1:n_p)) |>
        # Simulate true expression quantities
        mutate(
            d = pmap(list(feat_id, cell_type), function(g, p) {
                nonzero <- rbinom(1, 1, 0.1) == 1

                if (nonzero) {
                    return(tibble(
                        donor = 1:n_d,
                        q = rlnorm(n_d, meanlog = rnorm(1), sdlog = 0.5)
                    ))
                } else {
                    return(tibble(
                        donor = 1:n_d,
                        q = rlnorm(n_d, meanlog = rnorm(1), sdlog = 0.0)
                    ))
                }
            })
        ) |>
        unnest(d)

    # Simulate batch-effect
    true_sf <- tibble(
        donor = 1:n_d,
        batch_id = 1:n_d,
        sf = runif(n_d, 1.0, 10.0)
    )

    data <- data |>
        left_join(
            true_sf,
            by = "donor"
        ) |>
        # Generate counts
        mutate(
            d = pmap(
                list(feat_id, donor, cell_type, q, sf),
                function(g, d, p, q, sf) {
                    tibble(
                        obs_id = paste0(d, ":", p, "_", 1:n_o),
                        counts = as.integer(rnbinom(
                            n_o,
                            mu = q * sf,
                            size = 100
                        ))
                    )
                }
            )
        ) |>
        unnest(d) |>
        select(feat_id, donor, batch_id, cell_type, obs_id, counts)

    # Coerce relevant columns to a factor
    data[["donor"]] <- factor(data[["donor"]])
    data[["cell_type"]] <- factor(data[["cell_type"]])

    # Format into 'counts' and 'metadata'
    counts <- data |>
        select(feat_id, obs_id, counts) |>
        pivot_wider(names_from = obs_id, values_from = counts) |>
        column_to_rownames("feat_id") |>
        as.matrix() |>
        Matrix::t()

    metadata <- data %>%
        select(donor, batch_id, cell_type, obs_id) %>%
        distinct()

    # Compute size factors
    size_factors <- exp(disize::disize(
        design_formula = ~ cell_type + (1 | donor:cell_type),
        counts = counts,
        metadata = metadata
    ))

    expect_equal(
        object = as.vector(unname(size_factors)),
        expected = (true_sf$sf / sum(true_sf$sf)) * n_d,
        tolerance = 0.1
    )
})


test_that("large-simple-sc", {
    # Simulate data
    n_g <- 1000
    n_d <- 12
    n_p <- 3
    n_o <- 50

    data <- tibble(feat_id = 1:n_g) |>
        crossing(tibble(cell_type = 1:n_p)) |>
        # Simulate true expression quantities
        mutate(
            d = pmap(list(feat_id, cell_type), function(g, p) {
                nonzero <- rbinom(1, 1, 0.1) == 1

                if (nonzero) {
                    return(tibble(
                        donor = 1:n_d,
                        q = rlnorm(n_d, meanlog = rnorm(1), sdlog = 0.5)
                    ))
                } else {
                    return(tibble(
                        donor = 1:n_d,
                        q = rlnorm(n_d, meanlog = rnorm(1), sdlog = 0.0)
                    ))
                }
            })
        ) |>
        unnest(d)

    # Simulate batch-effect
    true_sf <- tibble(
        donor = 1:n_d,
        batch_id = 1:n_d,
        sf = runif(n_d, 1.0, 10.0)
    )

    data <- data |>
        left_join(
            true_sf,
            by = "donor"
        ) |>
        # Generate counts
        mutate(
            d = pmap(
                list(feat_id, donor, cell_type, q, sf),
                function(g, d, p, q, sf) {
                    tibble(
                        obs_id = paste0(d, ":", p, "_", 1:n_o),
                        counts = as.integer(rnbinom(
                            n_o,
                            mu = q * sf,
                            size = 100
                        ))
                    )
                }
            )
        ) |>
        unnest(d) |>
        select(feat_id, donor, batch_id, cell_type, obs_id, counts)

    # Coerce relevant columns to a factor
    data[["donor"]] <- factor(data[["donor"]])
    data[["cell_type"]] <- factor(data[["cell_type"]])

    # Format into 'counts' and 'metadata'
    counts <- data |>
        select(feat_id, obs_id, counts) |>
        pivot_wider(names_from = obs_id, values_from = counts) |>
        column_to_rownames("feat_id") |>
        as.matrix() |>
        Matrix::t()

    metadata <- data %>%
        select(donor, batch_id, cell_type, obs_id) %>%
        distinct()

    # Compute size factors
    size_factors <- exp(disize::disize(
        design_formula = ~ cell_type + (1 | cell_type:donor),
        counts = counts,
        metadata = metadata
    ))

    expect_equal(
        object = as.vector(unname(size_factors)),
        expected = (true_sf$sf / sum(true_sf$sf)) * n_d,
        tolerance = 0.1
    )
})

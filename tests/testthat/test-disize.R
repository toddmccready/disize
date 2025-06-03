library(dplyr)
library(tidyr)
library(purrr)

test_that("small-simple-bulk", {
    # Simulate data
    n_g <- 500
    n_d <- 2
    n_o <- 1

    data <- tibble(gene = 1:n_g) |>
        # Simulate true expression quantities
        mutate(
            d = pmap(list(gene), function(g) {
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
        batch = 1:n_d,
        sf = runif(n_d, 1.0, 10.0)
    )

    data <- data |>
        left_join(
            true_sf,
            by = "donor"
        ) |>
        # Generate counts
        mutate(
            d = pmap(list(gene, donor, q, sf), function(g, d, q, sf) {
                tibble(
                    cell_barcode = paste0(d, "_", 1:n_o),
                    counts = as.integer(rnbinom(n_o, mu = q * sf, size = 100))
                )
            })
        ) |>
        unnest(d) |>
        select(gene, donor, batch, counts)

    # Compute size factors
    size_factors <- exp(disize(
        design_formula = ~ (1 | donor),
        model_data = data,
        n_threads = ceiling(parallel::detectCores() / 2)
    ))

    expect_equal(
        object = as.vector(unname(size_factors)),
        expected = (true_sf$sf / sum(true_sf$sf)) * n_d,
        tolerance = 0.1
    )
})

test_that("large-simple-bulk", {
    # Simulate data
    n_g <- 500
    n_d <- 12
    n_o <- 1

    data <- tibble(gene = 1:n_g) |>
        # Simulate true expression quantities
        mutate(
            d = pmap(list(gene), function(g) {
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
        batch = 1:n_d,
        sf = runif(n_d, 1.0, 10.0)
    )

    data <- data |>
        left_join(
            true_sf,
            by = "donor"
        ) |>
        # Generate counts
        mutate(
            d = pmap(list(gene, donor, q, sf), function(g, d, q, sf) {
                tibble(
                    cell_barcode = paste0(d, "_", 1:n_o),
                    counts = as.integer(rnbinom(n_o, mu = q * sf, size = 100))
                )
            })
        ) |>
        unnest(d) |>
        select(gene, donor, batch, counts)

    # Compute size factors
    size_factors <- exp(disize(
        design_formula = ~ (1 | donor),
        model_data = data,
        n_threads = ceiling(parallel::detectCores() / 2)
    ))

    expect_equal(
        object = as.vector(unname(size_factors)),
        expected = (true_sf$sf / sum(true_sf$sf)) * n_d,
        tolerance = 0.1
    )
})

test_that("small-simple-sc", {
    # Simulate data
    n_g <- 500
    n_d <- 2
    n_p <- 3
    n_o <- 50

    data <- tibble(gene = 1:n_g) |>
        crossing(tibble(cell_type = 1:n_p)) |>
        # Simulate true expression quantities
        mutate(
            d = pmap(list(gene, cell_type), function(g, p) {
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
        batch = 1:n_d,
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
                list(gene, donor, cell_type, q, sf),
                function(g, d, p, q, sf) {
                    tibble(
                        cell_barcode = paste0(d, ":", p, "_", 1:n_o),
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
        select(gene, donor, batch, cell_type, counts)

    # Compute size factors
    size_factors <- exp(disize(
        design_formula = ~ cell_type + (1 | donor:cell_type),
        model_data = data,
        n_threads = ceiling(parallel::detectCores() / 2)
    ))

    expect_equal(
        object = as.vector(unname(size_factors)),
        expected = (true_sf$sf / sum(true_sf$sf)) * n_d,
        tolerance = 0.1
    )
})


test_that("large-simple-sc", {
    # Simulate data
    n_g <- 500
    n_d <- 6
    n_p <- 3
    n_o <- 50

    data <- tibble(gene = 1:n_g) |>
        crossing(tibble(cell_type = 1:n_p)) |>
        # Simulate true expression quantities
        mutate(
            d = pmap(list(gene, cell_type), function(g, p) {
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
        batch = 1:n_d,
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
                list(gene, donor, cell_type, q, sf),
                function(g, d, p, q, sf) {
                    tibble(
                        cell_barcode = paste0(d, ":", p, "_", 1:n_o),
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
        select(gene, donor, batch, cell_type, counts)

    # Compute size factors
    size_factors <- exp(disize(
        design_formula = ~ cell_type + (1 | donor:cell_type),
        model_data = data,
        n_threads = 8
    ))

    expect_equal(
        object = as.vector(unname(size_factors)),
        expected = (true_sf$sf / sum(true_sf$sf)) * n_d,
        tolerance = 0.1
    )
})

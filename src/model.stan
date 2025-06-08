data {
    // Dimensions ----
    int<lower=1> n_obs; // # of observations !

    int<lower=0> n_int; // # of distinct intercept parameters !
    int<lower=0> n_fe; // # of distinct fixed-effects parameters !
    int<lower=0> n_re; // # of distinct random-effects parameters !

    int<lower=0> n_nz_int; // # of nonzero elements in the intercept design matrix !
    int<lower=0> n_nz_fe; // # of nonzero elements in the fixed-effects design matrix !
    int<lower=0> n_nz_re; // # of nonzero elements in the random-effects design matrix !

    int<lower=0> n_re_terms; // # of random-effects terms in design formula !

    int<lower=1> n_batches; // # of batches !
    int<lower=1> n_feats; // # of original features !


    // Index Variables ----
    array[n_re] int<lower=1, upper=n_re_terms> re_id; // Random-effects term membership !
    array[n_obs] int<lower=1, upper=n_batches> batch_id; // Batch membership !
    array[n_obs] int<lower=1, upper=n_feats> feat_id; // Feature membership !


    // Response ----
    array[n_obs] int<lower=0> counts;


    // Design Matrices ----
    tuple(vector[n_nz_int], array[n_nz_int] int<lower=1>, array[n_obs + 1] int<lower=1>) int_design; // (data, col_idx, row_start_idx) !
    tuple(vector[n_nz_fe], array[n_nz_fe] int<lower=1>, array[n_obs + 1] int<lower=1>) fe_design; // (data, col_idx, row_start_idx) !
    tuple(vector[n_nz_re], array[n_nz_re] int<lower=1>, array[n_obs + 1] int<lower=1>) re_design; // (data, col_idx, row_start_idx) !
}
parameters {
    // Shrinkage ----
    real<lower=0> tau; // Global shrinkage parameter
    vector<lower=0>[n_fe] lambda; // Local shrinkage parameters


    // Feature Expression ----
    vector[n_int] int_coefs; // Intercept coefficients
    vector[n_fe] fe_coefs; // Fixed-effect coefficients
    vector[n_re] re_coefs; // Random-effect coefficients

    vector<lower=0>[n_re_terms] re_sigma; // Random-effect std-devs


    // Size Factors ----
    simplex[n_batches] raw_sf;


    // Feature-level Dispersion
    real<lower=0> iodisp; // TODO
}
transformed parameters {
    // Size Factors ----
    vector[n_batches] sf = log(raw_sf) + log(n_batches);
}
model {
    vector[n_obs] log_mu;

    // Priors ----
    // Higher Horseshoe
    lambda ~ cauchy(0, tau);
    re_sigma ~ cauchy(0, tau);

    // Random-effects
    for (i in 1:n_re) {
        re_coefs[i] ~ normal(0, re_sigma[re_id[i]]);
    }
    // Lower Horseshoe
    fe_coefs ~ normal(0, lambda);


    // Computing gene expression quantity
    log_mu = csr_matrix_times_vector(n_obs, n_int, int_design.1, int_design.2, int_design.3, int_coefs);
    log_mu += csr_matrix_times_vector(n_obs, n_fe, fe_design.1, fe_design.2, fe_design.3, fe_coefs);
    log_mu += csr_matrix_times_vector(n_obs, n_re, re_design.1, re_design.2, re_design.3, re_coefs);


    for (i in 1:n_obs) {
        // Adjusting for batch-effect
        log_mu += sf[batch_id[i]];

        // Likelihood ----
        counts[i] ~ neg_binomial_2_log(log_mu[i], iodisp);
    }
}

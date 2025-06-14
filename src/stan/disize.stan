data {
    // Dimensions ----
    int<lower=1> n_obs; // # of observations

    int<lower=0> n_int; // # of distinct intercept parameters
    int<lower=0> n_fe; // # of distinct fixed-effects parameters
    int<lower=0> n_re; // # of distinct random-effects parameters

    int<lower=0> n_nz_int; // # of nonzero elements in the intercept design matrix
    int<lower=0> n_nz_fe; // # of nonzero elements in the fixed-effects design matrix
    int<lower=0> n_nz_re; // # of nonzero elements in the random-effects design matrix

    int<lower=0> n_re_terms; // # of random-effects terms in design formula

    int<lower=1> n_batches; // # of batches
    int<lower=1> n_feats; // # of original features


    // Index Variables ----
    array[n_re] int<lower=1, upper=n_re_terms> re_id; // Random-effects term membership
    array[n_obs] int<lower=1, upper=n_batches> batch_id; // Batch membership
    array[n_obs] int<lower=1, upper=n_feats> feat_id; // Feature membership


    // Response ----
    array[n_obs] int<lower=0> counts;


    // Design Matrices ----
    vector[n_nz_int] int_design_x;
    array[n_nz_int] int int_design_j;
    array[n_obs + 1] int int_design_p;

    vector[n_nz_fe] fe_design_x;
    array[n_nz_fe] int fe_design_j;
    array[n_obs + 1] int fe_design_p;

    vector[n_nz_re] re_design_x;
    array[n_nz_re] int re_design_j;
    array[n_obs + 1] int re_design_p;
}
parameters {
    // Shrinkage for Fixed Effects (Horseshoe) ----
    real<lower=0> tau;          // Global shrinkage parameter
    vector<lower=0>[n_fe] lambda; // Local shrinkage parameters

    // Feature Expression ----
    vector[n_int] int_coefs;  // Intercept coefficients
    vector[n_fe] fe_coefs;    // Fixed-effect coefficients

    // Non-Centered Random Effects ----
    vector[n_re] z_re;        // Standardized ("raw") random-effect coefficients
    vector<lower=0>[n_re_terms] re_sigma; // Random-effect std-devs

    // Size Factors ----
    simplex[n_batches] raw_sf;

    // Feature-level Dispersion ----
    real<lower=0> iodisp;
}
transformed parameters {
    // Size Factors ----
    vector[n_batches] sf = log(raw_sf) + log(n_batches);

    // Actual Random Effects (Scaled) ----
    vector[n_re] re_coefs;
    for (i in 1:n_re) {
        re_coefs[i] = z_re[i] * (re_sigma[re_id[i]] * tau);
    }
}
model {
    vector[n_obs] log_mu;

    // --- Priors ---
    z_re ~ std_normal();

    // Horseshoe prior for fixed effects
    lambda ~ cauchy(0, 1);
    fe_coefs ~ normal(0, lambda * tau);

    // --- Linear Predictor ---
    log_mu = csr_matrix_times_vector(n_obs, n_int, int_design_x, int_design_j, int_design_p, int_coefs);
    if (n_fe != 0) {
        log_mu += csr_matrix_times_vector(n_obs, n_fe, fe_design_x, fe_design_j, fe_design_p, fe_coefs);
    }
    if (n_re != 0) {
        log_mu += csr_matrix_times_vector(n_obs, n_re, re_design_x, re_design_j, re_design_p, re_coefs);
    }

    // --- Likelihood ----
    for (i in 1:n_obs) {
        // Adjusting for batch-effect
        log_mu[i] += sf[batch_id[i]];
        counts[i] ~ neg_binomial_2_log(log_mu[i], iodisp);
    }
}

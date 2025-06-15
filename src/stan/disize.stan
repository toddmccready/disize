data {
    // Dimensions ----
    int<lower=1> n_obs; // # of observations
    int<lower=1> n_feats; // # of features

    int<lower=0> n_fe; // # of fixed-effects per-gene
    int<lower=0> n_re; // # of random-effects per-gene

    int<lower=0> n_nz_re; // # of nonzero elements in the random-effects design matrix
    int<lower=0> n_re_terms; // # of random-effects terms

    int<lower=1> n_batches; // # of batches

    // Design Matrices ----
    matrix[n_obs, n_fe] fe_design; // fixed-effects model matrix

    // Random-effects model matrix
    vector[n_nz_re] re_design_x;
    array[n_nz_re] int re_design_j;
    array[n_obs + 1] int re_design_p;

    // Index Variables ----
    array[n_re] int<lower=1, upper=n_re_terms> re_id; // random-effects term membership
    array[n_obs] int<lower=1, upper=n_batches> batch_id; // batch membership

    // Response ----
    array[n_feats, n_obs] int<lower=0> counts; // counts matrix
}
parameters {
    // Feature Expression ----
    vector[n_feats] intercept;

    matrix[n_fe, n_feats] fe_coefs;

    matrix[n_re, n_feats] raw_re_coefs;
    matrix<lower=0>[n_re_terms, n_feats] re_sigma;

    // Batch Effect ----
    simplex[n_batches] raw_sf;

    // Shrinkage ----
    real<lower=0> tau;

    // Inverse Overdispersion ----
    real<lower=0> iodisp;
}
transformed parameters {
    // Size Factors ----
    vector[n_batches] sf = log(raw_sf) + log(n_batches);

    // Random Effects ----
    matrix[n_re, n_feats] re_coefs;
    for (f_i in 1:n_feats) {
        for (r_i in 1:n_re) {
            re_coefs[r_i, f_i] = raw_re_coefs[r_i, f_i] * (re_sigma[re_id[r_i], f_i] * tau);
        }
    }
}
model {
    vector[n_obs] log_mu;

    // Priors ----
    raw_re_coefs ~ std_normal();

    // Horseshoe prior for fixed effects
    lambda ~ cauchy(0, 1);
    fe_coefs ~ normal(0, lambda * tau);

    for (f_i in 1:n_feats) {
        // Estimated Feature Expression ----
        log_mu = intercept[f_i];

        if (n_fe != 0) {
            log_mu += fe_design * col(fe_coefs, f_i);
        }

        if (n_re != 0) {
            log_mu += csr_matrix_times_vector(n_obs, n_re, re_design_x, re_design_j, re_design_p, col(re_coefs, f_i));
        }

        // Batch-effect Adjustment
        for (o_i in 1:n_obs) {
            log_mu[o_i] += sf[batch_id[o_i]];
        }

        // Likelihood ----
        counts[f_i] ~ neg_binomial_2_log(log_mu, iodisp);
    }
}

data {
    // Dimensions ----
    int<lower=0> n_int; // # of distinct intercept parameters
    int<lower=0> n_fe; // # of distinct fixed-effects parameters
    int<lower=0> n_re; // # of distinct random-effects parameters

    int<lower=0> n_re_terms; // # of random-effects terms in design formula

    int<lower=1> n_batches; // # of batches
    int<lower=1> n_feats; // # of original features

    int<lower=1> n_threads; // # of threads to use


    // Index Variables ----
    array[n_re] int<lower=1, upper=n_re_terms> re_id; // Random-effects term membership
    array[n_obs] int<lower=1, upper=n_batches> batch_id; // Batch membership
    array[n_obs] int<lower=1, upper=n_feats> feat_id; // Feature membership


    // Response ----
    array[n_obs] int<lower=0> counts;


    // Thread-specific Data ----
    int<lower=1> reals_length;
    int<lower=1> ints_length;

    array[n_threads, reals_length] real reals;
    array[n_threads, ints_length] real ints;


    int<lower=1> phi_size = 1 // tau
                   + n_fe // lambda
                   + n_int // int_coefs
                   + n_fe // fe_coefs
                   + n_re // re_coefs
                   + n_re_terms // re_sigma
                   + n_batches // sf
                   + n_feats // iodisp
                   + 1 // iodisp_mu
                   + 1; // iodisp_sigma
}
parameters {
    // Shrinkage ----
    real<lower=0> tau; // Global shrinkage parameter
    vector<lower=0>[n_fe] lambda; // Local shrinkage parameters for fixed-effects

    // Feature Expression ----
    vector[n_int] int_coefs;  // Intercept coefficients
    vector[n_fe] fe_coefs;    // Fixed-effects coefficients
    vector[n_re] z_re; // Raw random-effects coefficients
    vector<lower=0>[n_re_terms] re_sigma; // Random-effect std-devs

    // Size Factors ----
    simplex[n_batches] raw_sf;

    // Feature-level Dispersion ----
    vector<lower=0>[n_feats] iodisp;

    real iodisp_mu;
    real<lower=0> iodisp_sigma;
}
transformed parameters {
    // Size Factors ----
    vector[n_batches] sf = log(raw_sf) + log(n_batches);

    // Feature Expression ----
    vector[n_re] re_coefs; // Random-effects coefficients
    for (i in 1:n_re) {
        re_coefs[i] = z_re[i] * (re_sigma[re_id[i]] * tau);
    }


    // Pack parameters ----
    vector[phi_size] = phi;

    int phi_idx = 1;

    phi[phi_idx] = tau; // Pack tau
    phi_idx += 1;

    phi[phi_idx:(phi_idx + n_fe - 1)] = lambda; // Pack lambda
    phi_idx += n_fe;

    phi[phi_idx:(phi_idx + n_int - 1)] = int_coefs; // Pack int_coefs
    phi_idx += n_int;

    phi[phi_idx:(phi_idx + n_fe - 1)] = fe_coefs; // Pack fe_coefs
    phi_idx += n_fe;

    phi[phi_idx:(phi_idx + n_re - 1)] = re_coefs; // Pack re_coefs
    phi_idx += n_re;

    phi[phi_idx:(phi_idx + n_re_terms - 1)] = re_sigma; // Pack re_sigma
    phi_idx += n_re_terms;

    phi[phi_idx:(phi_idx + n_batches - 1)] = sf; // Pack sf
    phi_idx += n_batches;

    phi[phi_idx:(phi_idx + n_feats - 1)] = iodisp; // Pack iodisp
    phi_idx += n_feats;

    phi[phi_idx] = iodisp_mu; // Pack iodisp_mu
    phi_idx += 1;

    phi[phi_idx] = iodisp_sigma; // Pack iodisp_sigma
    phi_idx += 1;
}
model {
    vector[n_obs] log_mu;

    // Priors ----
    z_re ~ std_normal();

    // Horseshoe prior for fixed effects
    lambda ~ cauchy(0, 1);
    fe_coefs ~ normal(0, lambda * tau);

    // Inverse overdispersion regularization
    iodisp ~ lognormal(iodisp_mu, iodisp_sigma);

    // Computing log feature expression quantities
    log_mu = csr_matrix_times_vector(n_obs, n_int, int_design_x, int_design_j, int_design_p, int_coefs);
    if (n_fe != 0) {
        log_mu += csr_matrix_times_vector(n_obs, n_fe, fe_design_x, fe_design_j, fe_design_p, fe_coefs);
    }
    if (n_re != 0) {
        log_mu += csr_matrix_times_vector(n_obs, n_re, re_design_x, re_design_j, re_design_p, re_coefs);
    }

    // Likelihood ----
    for (i in 1:n_obs) {
        // Adjusting for batch-effect
        log_mu[i] += sf[batch_id[i]];

        counts[i] ~ neg_binomial_2_log(log_mu[i], iodisp[feat_id[i]]);
    }
}

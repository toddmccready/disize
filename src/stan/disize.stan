functions {
    vector shard_ll(vector phi,
                        vector theta,
                        array[] real x_r,
                        array[] int x_i) {
        // Unpack Metadata ----
        // Shard Dimensions ----
        int n_obs_shard = x_i[1];
        int n_int = x_i[2];
        int n_fe = x_i[3];
        int n_re = x_i[4];
        int n_nz_int_shard = x_i[5];
        int n_nz_fe_shard = x_i[6];
        int n_nz_re_shard = x_i[7];
        int n_batches = x_i[8];
        int n_feats = x_i[9];

        // Unpack Parameters ----
        int pos = 1;
        vector[n_int] int_coefs = phi[pos:(pos + n_int - 1)];
        pos += n_int;
        vector[n_fe] fe_coefs = phi[pos:(pos + n_fe - 1)];
        pos += n_fe;
        vector[n_re] re_coefs = phi[pos:(pos + n_re - 1)];

        pos += n_re;
        vector[n_batches] sf = phi[pos:(pos + n_batches - 1)];
        pos += n_batches;
        real iodisp = phi[pos];

        // Unpack Shard Data ----
        // Integers ----
        int i_pos = 10;
        array[n_obs_shard] int batch_id_shard = x_i[i_pos:(i_pos + n_obs_shard - 1)];
        i_pos += n_obs_shard;
        array[n_obs_shard] int feat_id_shard = x_i[i_pos:(i_pos + n_obs_shard - 1)];
        i_pos += n_obs_shard;
        array[n_obs_shard] int counts_shard = x_i[i_pos:(i_pos + n_obs_shard - 1)];
        i_pos += n_obs_shard;

        array[n_nz_int_shard] int int_design_j = x_i[i_pos:(i_pos + n_nz_int_shard - 1)];
        i_pos += n_nz_int_shard;
        array[n_obs_shard + 1] int int_design_p = x_i[i_pos:(i_pos + n_obs_shard)];
        i_pos += n_obs_shard + 1;

        array[n_nz_fe_shard] int fe_design_j;
        array[n_obs_shard + 1] int fe_design_p;
        if (n_fe > 0) {
            fe_design_j = x_i[i_pos:(i_pos + n_nz_fe_shard - 1)];
            i_pos += n_nz_fe_shard;
            fe_design_p = x_i[i_pos:(i_pos + n_obs_shard)];
            i_pos += n_obs_shard + 1;
        }

        array[n_nz_re_shard] int re_design_j;
        array[n_obs_shard + 1] int re_design_p;
        if (n_re > 0) {
            re_design_j = x_i[i_pos:(i_pos + n_nz_re_shard - 1)];
            i_pos += n_nz_re_shard;
            re_design_p = x_i[i_pos:(i_pos + n_obs_shard)];
        }

        // Reals ----
        int r_pos = 1;
        vector[n_nz_int_shard] int_design_x = to_vector(x_r[r_pos:(r_pos + n_nz_int_shard - 1)]);
        r_pos += n_nz_int_shard;

        if (n_fe > 0) {
            vector[n_nz_fe_shard] fe_design_x = to_vector(x_r[r_pos:(r_pos + n_nz_fe_shard - 1)]);
            r_pos += n_nz_fe_shard;
        }

        if (n_re > 0) {
            vector[n_nz_re_shard] re_design_x = to_vector(x_r[r_pos:(r_pos + n_nz_re_shard - 1)]);
        }

        // Compute negative binomial likelihood
        vector[n_obs_shard] log_mu;
        real log_lik = 0;

        // Compute log_mu from design matrices
        log_mu = csr_matrix_times_vector(n_obs_shard, n_int, int_design_x, int_design_j, int_design_p, int_coefs);
        if (n_fe > 0) {
            log_mu += csr_matrix_times_vector(n_obs_shard, n_fe, fe_design_x, fe_design_j, fe_design_p, fe_coefs);
        }
        if (n_re > 0) {
            log_mu += csr_matrix_times_vector(n_obs_shard, n_re, re_design_x, re_design_j, re_design_p, re_coefs);
        }

        // Likelihood loop
        for (i in 1:n_obs_shard) {
            // Adjust for batch-effect size factor
            log_mu[i] += sf[batch_id_shard[i]];

            log_lik += neg_binomial_2_log_lpmf(counts_shard[i] | log_mu[i], iodisp);
        }

        return [log_lik]';
    }
}
data {
    // Dimensions ----
    int<lower=0> n_int;       // # of distinct intercept parameters
    int<lower=0> n_fe;        // # of distinct fixed-effects parameters
    int<lower=0> n_re;        // # of distinct random-effects parameters
    int<lower=0> n_re_terms;  // # of random-effects terms in design formula
    int<lower=1> n_batches;   // # of batches
    int<lower=1> n_feats;     // # of original features
    int<lower=1> n_threads;   // # of threads to use

    // Index Variables ----
    array[n_re] int<lower=1, upper=n_re_terms> re_id;

    // Thread-specific Data ----
    int<lower=1> reals_per_thread;
    int<lower=1> ints_per_thread;
    array[n_threads, reals_per_thread] real x_r;
    array[n_threads, ints_per_thread] int x_i;
}
transformed data {
    array[n_threads] vector[0] theta;
}
parameters {
    // Shrinkage ----
    real<lower=0> tau;
    vector<lower=0>[n_fe] lambda;

    // Feature Expression ----
    vector[n_int] int_coefs;
    vector[n_fe] fe_coefs;
    vector[n_re] z_re; // Raw random-effects coefficients
    vector<lower=0>[n_re_terms] re_sigma;

    // Size Factors ----
    simplex[n_batches] raw_sf;

    // Feature-level Dispersion ----
    real<lower=0> iodisp;
}
transformed parameters {
    // Size Factors ----
    vector[n_batches] sf = log(raw_sf) + log(n_batches);

    // Feature Expression ----
    vector[n_re] re_coefs;
    for (i in 1:n_re) {
        re_coefs[i] = z_re[i] * (re_sigma[re_id[i]] * tau);
    }

    vector[n_int + n_fe + n_re + n_batches + 1] phi =
        append_row(int_coefs,
            append_row(fe_coefs,
                append_row(re_coefs,
                    append_row(sf, iodisp))));
}
model {
    // Priors ----
    z_re ~ std_normal();

    // Horseshoe prior for fixed effects
    lambda ~ cauchy(0, 1);
    fe_coefs ~ normal(0, lambda * tau);

    // Likelihood ----
    target += sum(map_rect(shard_ll, phi, theta, x_r, x_i));
}

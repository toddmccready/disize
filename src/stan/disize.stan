functions {
    vector shard_ll(vector phi,
                        vector theta,
                        array[] real x_r,
                        array[] int x_i) {
        // Unpack Metadata ----
        // Shard Dimensions ----
        int n_obs = x_i[1];
        int n_int = x_i[2];
        int n_fe = x_i[3];
        int n_re = x_i[4];
        int n_nz_int = x_i[5];
        int n_nz_fe = x_i[6];
        int n_nz_re = x_i[7];
        int n_batches = x_i[8];
        int n_feats = x_i[9];

        // Index Parameters ----
        int pos = 1;

        tuple(int, int) int_coefs = (pos, pos + n_int - 1);
        pos += n_int;

        tuple(int, int) fe_coefs = (pos, pos + n_fe - 1);
        pos += n_fe;

        tuple(int, int) re_coefs = (pos, pos + n_re - 1);
        pos += n_re;

        tuple(int, int) sf = (pos, pos + n_batches - 1);
        pos += n_batches;

        int iodisp = pos;

        // Index Shard Data ----
        // Integers ----
        int i_pos = 10;
        tuple(int, int) batch_id = (i_pos, i_pos + n_obs - 1);
        i_pos += n_obs;

        tuple(int, int) feat_id = (i_pos, i_pos + n_obs - 1);
        i_pos += n_obs;

        tuple(int, int) counts = (i_pos, i_pos + n_obs - 1);
        i_pos += n_obs;

        tuple(int, int) int_design_j = (i_pos, i_pos + n_nz_int - 1);
        i_pos += n_nz_int;

        tuple(int, int) int_design_p = (i_pos, i_pos + n_obs);
        i_pos += n_obs + 1;

        tuple (int, int) fe_design_j;
        tuple (int, int) fe_design_p;
        if (n_fe > 0) {
            fe_design_j.1 = i_pos; fe_design_j.2 = i_pos + n_nz_fe - 1;
            i_pos += n_nz_fe;

            fe_design_p.1 = i_pos; fe_design_p.2 = i_pos + n_obs;
            i_pos += n_obs + 1;
        }

        tuple (int, int) re_design_j;
        tuple (int, int) re_design_p;
        if (n_re > 0) {
            re_design_j.1 = i_pos; re_design_j.2 = i_pos + n_nz_re - 1;
            i_pos += n_nz_re;

            re_design_p.1 = i_pos; re_design_p.2 = i_pos + n_obs;
        }

        // Reals ----
        int r_pos = 1;
        tuple(int, int) int_design_x = (r_pos, r_pos + n_nz_int - 1);
        r_pos += n_nz_int;

        tuple(int, int) fe_design_x;
        if (n_fe > 0) {
            fe_design_x.1 = r_pos; fe_design_x.2 = r_pos + n_nz_fe - 1;
            r_pos += n_nz_fe;
        }

        tuple(int, int) re_design_x;
        if (n_re > 0) {
            re_design_x.1 = r_pos; re_design_x.2 = r_pos + n_nz_re - 1;
        }

        // Compute negative binomial likelihood
        vector[n_obs] log_mu;
        real log_lik = 0;

        // Compute log_mu from design matrices
        log_mu = csr_matrix_times_vector(
            n_obs,
            n_int,
            to_vector(x_r[(int_design_x.1):(int_design_x.2)]),
            x_i[(int_design_j.1):(int_design_j.2)],
            x_i[(int_design_p.1):(int_design_p.2)],
            phi[(int_coefs.1):(int_coefs.2)]
        );
        if (n_fe > 0) {
            log_mu += csr_matrix_times_vector(
                n_obs,
                n_fe,
                to_vector(x_r[(fe_design_x.1):(fe_design_x.2)]),
                x_i[(fe_design_j.1):(fe_design_j.2)],
                x_i[(fe_design_p.1):(fe_design_p.2)],
                phi[(fe_coefs.1):(fe_coefs.2)]
            );
        }
        if (n_re > 0) {
            log_mu += csr_matrix_times_vector(
                n_obs,
                n_re,
                to_vector(x_r[(re_design_x.1):(re_design_x.2)]),
                x_i[(re_design_j.1):(re_design_j.2)],
                x_i[(re_design_p.1):(re_design_p.2)],
                phi[(re_coefs.1):(re_coefs.2)]
            );
        }

        // Likelihood loop
        for (i in 1:n_obs) {
            // Adjust for batch-effect size factor
            log_mu[i] += phi[(sf.1):(sf.2)][x_i[(batch_id.1):(batch_id.2)][i]];

            log_lik += neg_binomial_2_log_lpmf(x_i[(counts.1):(counts.2)][i] | log_mu[i], phi[iodisp]);
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

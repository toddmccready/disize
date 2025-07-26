functions {
  real partial_posterior(// Thread-specific ----
                         array[,] int counts_slice, int start, int end,
                         // Shared ----
                         int n_obs, int n_fe, matrix fe_design, int n_re,
                         vector re_design_x, array[] int re_design_j,
                         array[] int re_design_p, array[] int re_id,
                         array[] int batch_id, vector intercept,
                         matrix raw_fe_coefs, matrix raw_re_coefs,
                         matrix re_lambda, matrix fe_lambda, vector fe_tau,
                         vector re_tau, vector sf, real iodisp) {
    real log_prob = 0;
    vector[n_obs] log_mu;
    vector[n_re] re_coefs_col;
    
    int n_slice_feats = end - start + 1; // # of features in slice
    for (i in 1 : n_slice_feats) {
      int feat_i = start + i - 1;
      
      // Priors ----
      // Half-cauchy prior over random-effects variance (to fulfill horseshoe)
      log_prob += cauchy_lpdf(re_lambda[ : , feat_i] | 0, 1);
      
      // Horseshoe prior over fixed-effects
      log_prob += cauchy_lpdf(fe_lambda[ : , feat_i] | 0, 1);
      log_prob += normal_lpdf(raw_fe_coefs[ : , feat_i] | 0, 1);
      
      // Normal prior over (raw) random-effects
      log_prob += std_normal_lpdf(raw_re_coefs[ : , feat_i]);
      
      // Covariate Effects ----
      log_mu = rep_vector(intercept[feat_i], n_obs);
      
      if (n_fe != 0) {
        log_mu += fe_design
                  * (raw_fe_coefs[ : , feat_i]
                     .* (fe_lambda[ : , feat_i] .* fe_tau));
      }
      
      if (n_re != 0) {
        // Scaling standardized random-effects coefficients
        for (re_i in 1 : n_re) {
          re_coefs_col[re_i] = raw_re_coefs[re_i, feat_i]
                               * (re_lambda[re_id[re_i], feat_i]
                                  * re_tau[re_id[re_i]]);
        }
        
        log_mu += csr_matrix_times_vector(n_obs, n_re, re_design_x,
                                          re_design_j, re_design_p,
                                          re_coefs_col);
      }
      
      // Batch-effect Adjustment ----
      for (obs_i in 1 : n_obs) {
        log_mu[obs_i] += sf[batch_id[obs_i]];
      }
      
      // Likelihood ----
      log_prob += neg_binomial_2_log_lpmf(counts_slice[i] | log_mu, iodisp);
    }
    return log_prob;
  }
}
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
  // fixed-effects
  matrix[n_obs, n_fe] fe_design;
  
  // random-effects
  vector[n_nz_re] re_design_x;
  array[n_nz_re] int re_design_j;
  array[n_obs + 1] int re_design_p;
  
  // Index Variables ----
  array[n_re] int<lower=1, upper=n_re_terms> re_id; // random-effects term membership
  array[n_obs] int<lower=1, upper=n_batches> batch_id; // batch membership
  
  // Response ----
  array[n_feats, n_obs] int<lower=0> counts; // counts matrix
  
  // Configuration for Multi-threading ----
  int<lower=1> grainsize;
}
parameters {
  // Feature Expression ----
  vector[n_feats] intercept;
  
  matrix[n_fe, n_feats] raw_fe_coefs;
  
  matrix[n_re, n_feats] raw_re_coefs;
  matrix<lower=0>[n_re_terms, n_feats] re_lambda;
  
  // Batch Effect ----
  simplex[n_batches] raw_sf;
  
  // Shrinkage ----
  vector<lower=0>[n_fe] fe_tau;
  vector<lower=0>[n_re_terms] re_tau;
  matrix<lower=0>[n_fe, n_feats] fe_lambda;
  
  // Inverse Overdispersion ----
  real<lower=0> iodisp;
}
transformed parameters {
  // Size Factors ----
  vector[n_batches] sf = log(raw_sf) + log(n_batches);
}
model {
  // Parallel posterior eval ----
  target += reduce_sum(partial_posterior, counts, grainsize,
                       // Shared ----
                       n_obs, n_fe, fe_design, n_re, re_design_x,
                       re_design_j, re_design_p, re_id, batch_id, intercept,
                       raw_fe_coefs, raw_re_coefs, re_lambda, fe_lambda,
                       fe_tau, re_tau, sf, iodisp);
}

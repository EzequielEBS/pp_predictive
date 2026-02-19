functions {
  matrix post_par_conj_lm(vector mu, matrix S, real a, real b, matrix V, matrix X, vector y) {
    // V: nxn data covariance matrix
    // mu: px1 prior mean
    // S: pxp prior covariance matrix
    // a: scalar, shape parameter of the inverse-gamma distribution
    // b: scalar, scale parameter of the inverse-gamma distribution
    // X: nxp design matrix
    // y: nx1 response vector
    
    int n = rows(V);
    int p = rows(mu);
    matrix[p,p] S_inv = inverse(S);
    matrix[n,n] V_inv = inverse(V);
    matrix[p,p] S_star = inverse(S_inv + X' * V_inv * X);
    vector[p] mu_star = S_star * (S_inv * mu + X' * V_inv * y);
    real a_star = a + 0.5*n;
    real b_star = b + 0.5 * (y' * V_inv * y + mu' * S_inv * mu - 
                           mu_star' * inverse(S_star) * mu_star);
    matrix[p,p+3] out;
    out[1:p, 1] = mu_star;
    out[1:p, 2:(p+1)] = S_star;
    out[1, p+2] = a_star;
    out[2, p+3] = b_star;
    
    return out;
  }
  
  matrix pp_hyper_conj_lm(real eta, vector mu, matrix S, real a, real b, matrix X0, vector y0) {
    // eta: scalar, power parameter
    // mu: px1 prior mean
    // S: pxp prior covariance matrix
    // a: scalar, shape parameter of the inverse-gamma distribution
    // b: scalar, scale parameter of the inverse-gamma distribution
    // X0: nxp design matrix for the historical data
    // y0: nx1 response vector for the historical data
    
    int n0 = rows(X0);
    int p = rows(mu);
    matrix[p,p+3] out;
    if (eta == 0) {
      out[1:p, 1] = mu;
      out[1:p, 2:(p+1)] = S;
      out[1, p+2] = a;
      out[2, p+3] = b;
    } else {
      matrix[n0, n0] V0 = 1/eta * diag_matrix(rep_vector(1.0, n0));
      real a0 = a + 0.5*n0*(eta-1);
      out = post_par_conj_lm(mu, S, a0, b, V0, X0, y0);
    }
    return out;
  }
}

data {
  // int<lower=1> K; // number of historical data
  int<lower=0> n0; // number of observations of each hist data
  int<lower=0> n; // number of current observations 
  int<lower=0> p; // number of current covariates
  matrix[n0,p] X0;
  vector[n0] y0;
  matrix[n,p] X;
  vector[n] y;
  real<lower=0> a; // prior shape
  real<lower=0> b; // prior scale
  matrix[p,p] V; // prior variance
  vector[p] mu; // prior mean
  real<lower=0> tilde_a; // shape for eta prior
  real<lower=0> tilde_b; // shape for eta prior
}

parameters {
  real<lower=0, upper=1> eta;
  // vector[p] beta;
  // real<lower=0> s2;
}

transformed parameters {
  // Build default parameters
  matrix[p, p+3] pp_hyper = pp_hyper_conj_lm(eta, mu, V, a, b, X0, y0);
  real<lower=0> nu0 = pp_hyper[1, p+2];
  matrix[p,p] inv_Lambda0 = pp_hyper[1:p, 2:(p+1)];
  matrix[p,p] Lambda0 = inverse_spd(inv_Lambda0);
  vector[p] tilde_beta0 = pp_hyper[1:p, 1];
  real<lower=0> H0 = pp_hyper[2, p+3];
  
  matrix[p, p+3] post_par = post_par_conj_lm(tilde_beta0, inv_Lambda0, nu0, H0, 
                                             diag_matrix(rep_vector(1.0, n)), 
                                             X, y);
  vector[p] tilde_beta = post_par[1:p, 1];
  matrix[p,p] inv_Lambda = post_par[1:p, 2:(p+1)];
  matrix[p,p] Lambda = inverse_spd(inv_Lambda);
  real<lower=0> nu = post_par[1, p+2];
  real<lower=0> H = post_par[2, p+3];
}

model {
  // prior
  target += beta_lpdf(eta| tilde_a, tilde_b);
  
  target += lgamma(nu) + 
            0.5 * log_determinant(Lambda0) + 
            nu0 * log(H0) -
            0.5 * log_determinant(Lambda) - 
            lgamma(nu0) -
            nu * log(H);
}

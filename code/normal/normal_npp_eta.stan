functions {
  vector post_par_conj_normal(real m, real v, real a, real b, real w, vector y) {
    int n = rows(y);
    real v_star = 1/(1/w*n + v);
    real m_star = v_star * (m/v + 1/w * sum(y));
    real a_star = a + 0.5*n;
    real b_star = b + 0.5 * (y' * y + 1/v * m^2 - 
                             1/v_star * m_star^2);
    vector[4] out;
    out[1] = m_star;
    out[2] = v_star;
    out[3] = a_star;
    out[4] = b_star;
    
    return out;
  }
  
  vector pp_hyper_conj_normal(real eta, real m, real v, real a, real b, vector y0) {
    int n0 = rows(y0);
    vector[4] out;
    if (eta == 0) {
      out[1] = m;
      out[2] = v;
      out[3] = a;
      out[4] = b;
    } else {
      real w0 = 1/eta;
      real a0 = a + 0.5*n0*(eta-1);
      out = post_par_conj_normal(m, v, a0, b, w0, y0);
    }
    return out;
  }
}

data {
  // int<lower=1> K; // number of historical data
  int<lower=0> n0; // number of observations of each hist data
  int<lower=0> n; // number of current observations 
  vector[n0] y0;
  vector[n] y;
  real<lower=0> a; // prior shape
  real<lower=0> b; // prior scale
  real<lower=0> v; // prior variance
  real m; // prior mean
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
  vector[4] pp_hyper = pp_hyper_conj_normal(eta, m, v, a, b, y0);
  real<lower=0> nu0 = pp_hyper[3];
  real inv_Lambda0 = pp_hyper[2];
  real Lambda0 = 1/inv_Lambda0;
  real tilde_beta0 = pp_hyper[1];
  real<lower=0> H0 = pp_hyper[4];
  
  vector[4] post_par = post_par_conj_normal(tilde_beta0, inv_Lambda0, nu0, H0, 
                                             1, y);
  real tilde_beta = post_par[1];
  real<lower=0> inv_Lambda = post_par[2];
  real<lower=0> Lambda = 1/inv_Lambda;
  real<lower=0> nu = post_par[3];
  real<lower=0> H = post_par[4];
}

model {
  // prior
  target += beta_lpdf(eta| tilde_a, tilde_b);
  
  target += lgamma(nu) + 
            0.5 * log(Lambda0) + 
            nu0 * log(H0) -
            0.5 * log(Lambda) - 
            lgamma(nu0) -
            nu * log(H);
}

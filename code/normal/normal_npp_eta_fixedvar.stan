// Fixed-variance analogue of code/normal/normal_npp_eta.stan.
//
// Implements, directly, the formulas for the fixed-variance normal NPP
// (theta | eta, D, D0 and the eta marginal posterior) as given by:
//
//   theta | eta, D, D0 ~ N(mu_tilde_eta, tau2_tilde_eta),
//     tau2_tilde_eta^-1 = tau2_bullet^-1 + n/sigma2,
//     mu_tilde_eta = tau2_tilde_eta * (mu_bullet/tau2_bullet + n*ybar/sigma2),
//
//   Z(mu, tau2) = (2*pi*tau2)^(1/2) * exp(mu^2 / (2*tau2)),
//
//   pi_NPP(eta | D, D0) propto pi_A(eta) * Z(mu_tilde_eta, tau2_tilde_eta)
//                                          / Z(mu_eta, tau2_eta),
//
// where (mu_eta, tau2_eta) are the same "bullet" update rule applied to
// the ORIGINAL prior (mu_A, tau2_A) using only the eta-power-discounted
// historical data D0 (i.e. mu_eta/tau2_eta play the role of "mu_bullet/
// tau2_bullet" for that first update, before D is seen):
//
//   tau2_eta^-1 = tau2_A^-1 + eta*n0/sigma2,
//   mu_eta = tau2_eta * (mu_A/tau2_A + eta*n0*ybar0/sigma2).
//
// Only eta is a sampled parameter (theta is analytically marginalized out
// of the target, exactly as in normal_npp_eta.stan for the unknown-
// variance/NIG model); theta | eta, D, D0 is then drawn exactly in
// generated quantities at every posterior draw of eta, so (eta, theta)
// together are joint draws from the NPP posterior p(theta, eta | D, D0) --
// eta by NUTS/HMC, theta by an exact closed-form conditional draw.

functions {
  // log Z(mu, tau2) = log[ (2*pi*tau2)^(1/2) * exp(mu^2 / (2*tau2)) ]
  real log_Z(real mu, real tau2) {
    return 0.5 * log(2 * pi() * tau2) + square(mu) / (2 * tau2);
  }
}

data {
  int<lower=0> n0;          // historical sample size
  int<lower=0> n;           // current sample size
  vector[n0] y0;             // historical data D0
  vector[n] y;                // current data D
  real m0;                    // mu_bullet: ORIGINAL prior mean for theta
  real<lower=0> v0;           // tau_bullet^2: ORIGINAL prior variance for theta
  real<lower=0> v;            // known/fixed variance sigma^2, shared by D0 and D
  real<lower=0> tilde_a;      // pi_A(eta) = Beta(tilde_a, tilde_b)
  real<lower=0> tilde_b;
}

transformed data {
  real ybar0 = mean(y0);
  real ybar  = mean(y);
}

parameters {
  real<lower=0, upper=1> eta;
}

transformed parameters {
  // theta | eta ~ N(mu_eta, tau2_eta): power prior for theta after
  // discounting D0 by eta, but before D is seen
  real inv_tau2_eta = 1 / v0 + eta * n0 / v;
  real<lower=0> tau2_eta = 1 / inv_tau2_eta;
  real mu_eta = tau2_eta * (m0 / v0 + eta * n0 * ybar0 / v);

  // theta | eta, D, D0 ~ N(mu_tilde_eta, tau2_tilde_eta): further updated
  // with the current data D
  real inv_tau2_tilde_eta = inv_tau2_eta + n / v;
  real<lower=0> tau2_tilde_eta = 1 / inv_tau2_tilde_eta;
  real mu_tilde_eta = tau2_tilde_eta * (mu_eta / tau2_eta + n * ybar / v);
}

model {
  target += beta_lpdf(eta | tilde_a, tilde_b);
  target += log_Z(mu_tilde_eta, tau2_tilde_eta) - log_Z(mu_eta, tau2_eta);
}

generated quantities {
  // theta, drawn exactly (closed form) at this posterior draw of eta --
  // the "Gibbs" companion draw. Together with eta, gives one joint draw
  // from p(theta, eta | D, D0).
  real theta = normal_rng(mu_tilde_eta, sqrt(tau2_tilde_eta));
}

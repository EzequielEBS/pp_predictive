data {
    int<lower=0> n0;
    vector[n0] y0;
    int<lower=0> n1;
    real m0;
    real<lower=0> v0;
    real<lower=0> a0;
    real<lower=0> b0;
    real<lower=0> eta;
}

parameters {
    real mu;
    real<lower=0> s2;
}

transformed parameters {
    real<lower=0> sigma;
    sigma = sqrt(s2);
}

model {
    for (i in 1:n0) {
        target += eta * normal_lpdf(y0[i] | mu, sigma);
    }

    target += normal_lpdf(mu | m0, sqrt(s2*v0));
    target += inv_gamma_lpdf(s2 | a0, b0);
}

generated quantities {
    vector[n1] y1;
    for (i in 1:n1) {
        y1[i] = normal_rng(mu, sqrt(s2));
    }
}
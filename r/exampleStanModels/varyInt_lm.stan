data {
  int<lower=1> n_cats;
  int<lower=0> n_samples;
  int<lower=1> cats[n_samples];
  real y[n_samples];
}
parameters {
  real<lower=9, upper=11> alpha_mu;
  real<lower=1, upper=3> alpha_tau;
  real alpha[n_cats];
  real<lower=0, upper=1.0> sigma;
}
model {
  for (i in 1:n_cats) {
    alpha[i] ~ normal(alpha_mu, alpha_tau);
  }
  y ~ normal(alpha[cats], sigma);
}

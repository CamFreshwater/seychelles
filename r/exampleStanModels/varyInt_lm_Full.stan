data {
  int<lower=1> n_cats;
  int<lower=0> n_samples;
  int<lower=1> cats[n_samples];
  vector[n_samples] x;
  vector[n_samples] y;
}
parameters {
  real<lower=9, upper=11> alpha_mu;
  real<lower=1, upper=3> alpha_tau;
  real alpha[n_cats];
  real b_x;
  real<lower=0, upper=1.0> sigma;
}
model {
  b_x ~ normal(0, 5);  
  for (i in 1:n_cats) {
    alpha[i] ~ normal(alpha_mu, alpha_tau);
  }

  for (i in 1:n_samples) {
    y[i] ~ normal(alpha[cats[i]], sigma);
  }
}
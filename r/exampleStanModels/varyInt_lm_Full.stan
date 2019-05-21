data {
  int<lower=1> n_cats; // number of categories
  int<lower=0> N; // rows of data
  int<lower=1> J; // n fixed covariates on mean
  int<lower=1> cats[N];
  matrix[N,J] X_ij; // covariate model matrix
  vector[N] y_i; // vector to hold observations
}
parameters {
  real<lower=9, upper=11> alpha_mu;
  real<lower=1, upper=3> alpha_tau;
  real alpha[n_cats];
  vector[J] b_j;
  real<lower=0, upper=1.0> sigma;
}
model {
  b_j ~ normal(0, 5);  
  for (i in 1:n_cats) {
    alpha[i] ~ normal(alpha_mu, alpha_tau);
  }
  for(n in 1:N) {
    y_i[n] ~ normal(alpha 

  // for (i in 1:n_cats) {
  //  alpha[i] ~ normal(alpha_mu, alpha_tau);
  // }
  // y_i ~ normal(X_ij * b_j, sigma);
}

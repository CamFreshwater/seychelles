data {
  int<lower=1> K; // number of fishers
  int<lower=1> J; // number of years
  int<lower=1> N; // rows of data
  int<lower=1, upper=K> fisher_id[N]; // vector of fisher indices
  int<lower=1, upper=J> yr_id[N]; // vector of yr indices
  vector[N] x1; //exp var 1
  vector[N] y; // vector to hold observations
}
parameters {
  vector[K] a_fisher_z; //group level deviate
  vector[J] a_yr_z; //group level deviate
  real b1;
  real mu_a; //global intercept 
  real<lower=0> sigma; //SD of ind. obs
  real<lower=0> sigma_a_fisher; //SD of group intercept
  real<lower=0> sigma_a_yr; //SD of group intercept
}
transformed parameters {
  vector[K] a_fisher; //group-level deviates
  vector[J] a_yr; //group-level deviates
  a_fisher = sigma_a_fisher * a_fisher_z;
  a_yr = sigma_a_yr * a_yr_z;
}
model {
  vector[N] mu; //linear predictor

  //group-level ints
  a_fisher_z ~ normal(0,1);
  a_yr_z ~ normal(0,1);
  
  //priors
  mu_a ~ normal(0, 5);
  b1 ~ normal(0, 5);
  sigma_a_yr ~ cauchy(0, 2.5);
  sigma_a_fisher ~ cauchy(0, 2.5);
  sigma ~ gamma(2, 0.1);

  for (i in 1:N) {
    mu[i] = mu_a + a_fisher[fisher_id[i]] + a_yr[yr_id[i]] + b1*x1[i];
  }

  //likelihood
  y ~ normal(mu, sigma);
}

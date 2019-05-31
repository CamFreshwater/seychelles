data {
  int<lower=1> F; // number of fishers
  int<lower=1> J; // number of years
  int<lower=1> N; // rows of data
  int<lower=1, upper=F> fisher_id[N]; // vector of fisher indices
  int<lower=1, upper=J> yr_id[N]; // vector of yr indices
  vector[N] x1; //exp var 1
  vector[N] x2; //exp var 2
  vector[N] y; // vector to hold observations
}
parameters {
  vector[F] a_fisher_z; //group level deviate
  vector[J] a_yr_z; //group level deviate
  real b1;
  real b2;
  real mu_a; //global intercept 
  real<lower=0> sigma; //SD of ind. obs
  real<lower=0> sigma_a_group; //SD of group intercept
  real<lower=2> nu;
}
transformed parameters {
  vector[F] a_fisher; //group-level deviates
  vector[J] a_year; //group-level deviates
  a_fisher = sigma_a_fisher * a_fisher_z;
  a_year = sigma_a_yr * a_yr_z;
}
model {
  vector[N] mu; //linear predictor

  //group-level ints
  a_fisher_z ~ normal(0,1);
  a_yr_z ~ normal(0,1);
  
  //priors
  mu_a ~ normal(0, 5);
  b1 ~ normal(0, 5);
  b2 ~ normal(0, 5);
  sigma_a_yr ~ student_t(5, 0, 3);
  sigma_a_fisher ~ student_t(5, 0, 3);
  sigma ~ student_t(5, 0, 3);

  for (i in 1:N) {
    mu[i] = mu_a + a_fisher[fisher_id[i]] + a_yr[yr_id[i]] + b1*x1[i] + b2*x2[i];
  }

  //likelihood
  y ~ student_t(nu, mu, sigma);
}

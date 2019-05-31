data {
  int<lower=1> J; // number of groups
  int<lower=1> N; // rows of data
  int<lower=1, upper=J> group_id[N]; // vector of group indices
  vector[N] x1;
  vector[N] x2;
  vector[N] y; // vector to hold observations
}
parameters {
  vector[J] a_group_z; //group level deviate
  real b1;
  real b2;
  real mu_a; //global intercept 
  real<lower=0> sigma; //SD of ind. obs
  real<lower=0> sigma_a_group; //SD of group intercept
  real<lower=2> nu;
}
transformed parameters {
  vector[J] a_group; //group-level deviates
  a_group = sigma_a_group * a_group_z;
}
model {
  vector[N] mu; //linear predictor

  //group-level ints
  a_group_z ~ normal(0,1);
  
  //priors
  nu ~ exponential(0.1);
  mu_a ~ normal(0, 5);
  b1 ~ normal(0, 5);
  b2 ~ normal(0, 5);
  sigma_a_group ~ student_t(5, 0, 3);
  sigma ~ student_t(5, 0, 3);

  for (i in 1:N) {
    mu[i] = mu_a + a_group[group_id[i]] + b1*x1[i] + b2*x2[i];
  }

  //likelihood
  y ~ student_t(nu, mu, sigma);
}

data {
  int<lower=1> J; // number of groups
  int<lower=1> N; // rows of data
  int<lower=1> K; // number of columns in model matrix
  int<lower=1, upper=J> id[N]; // vector of group indices
  matrix[N, K] X_ij; // covariate model matrix
  vector[N] y_i; // vector to hold observations
}
parameters {
  vector[K] gamma; //population-level regression coefficients
  vector<lower=0>[K] tau; //the standard deviation of the regression coefficients
  //implementing Matt's trick
  vector[K] beta_raw[J];
  real<lower=0> sigma; //standard deviation of the individual observations
}
transformed parameters {
  vector[K] beta[J]; //matrix of group-level regression coefficients
  //computing the group-level coefficient, based on non-centered parametrization based on section 22.6 STAN (v2.12) user's guide
  for(j in 1:J){
    beta[j] = gamma + tau .* beta_raw[j];
  }
}
model {
  vector[N] mu; //linear predictor
  //priors
  gamma ~ normal(0,5); //weakly informative priors on the regression coefficients
  tau ~ cauchy(0,2.5); //weakly informative priors, see section 6.9 in STAN user guide
  sigma ~ gamma(2,0.1); //weakly informative priors, see section 6.9 in STAN user guide
  for(j in 1:J){
   beta_raw[j] ~ normal(0,1); //fill the matrix of group-level regression coefficients 
  }
  for(n in 1:N){
    mu[n] = X_ij[n] * beta[id[n]]; //compute the linear predictor using relevant group-level regression coefficients 
  }
  //likelihood
  y_i ~ normal(mu,sigma);
}

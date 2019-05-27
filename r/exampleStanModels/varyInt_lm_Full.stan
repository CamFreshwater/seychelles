data {
  int<lower=1> J; // number of groups
  int<lower=1> N; // rows of data
  int<lower=1> K; // number of columns in model matrix
  int<lower=1, upper=J> id[N]; // vector of group indices
  matrix[N, K] X_ij; // covariate model matrix
  vector[N] y_i; // vector to hold observations
}
parameters {
  vector[K] gamma; // group level reg coefs
  vector<lower=0>[K] tau; // SD of reg chefs  

  vector[K] beta[J]; //matrix of group level reg coefs
  real<lower=0> sigma; // SD of ind. obs
}
model {
  vector[N] mu; // linear predictor
  
  gamma ~ normal(0,5); //weakly informative priors reg coef
  tau ~ cauchy(0,2.5); //weakly informative priors, see section 6.9 in STAN user guide
  sigma ~ gamma(2,0.1); //weakly informative priors, see section 6.9 in STAN user guide
  
  for(j in 1:J){
   beta[j] ~ normal(gamma,tau); //fill the matrix of group-level regression coefficients 
  }
  for(n in 1:N){
    mu[n] = X_ij[n] * beta[id[n]]; //compute the linear predictor using group-level reg coef 
  }

  y_i ~ normal(mu,sigma);
}

library(dplyr)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(12445)
n_cats <- 10
n <- 1000
alpha_mu <- 10
alpha_tau <- 2
beta_x <- 0.4
sigma <- 0.2

df <- data.frame(cat = sample(1:n_cats, n, replace = TRUE))
alpha <- data.frame(cat = 1:n_cats,
                    alpha = rnorm(n_cats, alpha_mu, alpha_tau))

df1 <- left_join(df, alpha, by = "cat") %>%
  mutate(x = rnorm(n, 0, 10)) %>% 
  mutate(y = alpha + (beta_x * x) + rnorm(n, 0, sigma),
         y2 = alpha + rnorm(n, 0, sigma)) %>%
  tbl_df()

dat <- list(n_cats = n_cats,
             n_samples = n,
             cats = df$cat,
             y = df1$y)

# library(lme4)
# lmer(y ~ (1 | cat) + x, df1) %>% coef()
# lmer(y2 ~ (1 | cat), df1) %>% coef()

# model
mod <- stan_model(here::here("r/exampleStanModels/varyInt_lm.stan"))

fitModel <- sampling(mod,
         data = dat,
         iter = 100, chains = 1, cores = 1,
         control = list(adapt_delta = 0.95, max_treedepth = 20))

mod <- '
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
'
model1 <- stan(model_code = mod, data = dat, iter = 100, cores = 1, chains = 1)

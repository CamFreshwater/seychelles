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

X_ij <- model.matrix(~x, data = df1)

dat <- list(n_cats = n_cats,
            n_samples = n,
            cats = df1$cat,
            X_ij = X_ij,
            y_i = df1$y)

library(lme4)
lmer(y ~ (1 | cat) + x, df1) %>% coef()
# lmer(y2 ~ (1 | cat), df1) %>% coef()

# model
mod <- stan_model(here::here("r/exampleStanModels/varyInt_lm_Full.stan"))

fitModel <- sampling(mod,
         data = dat,
         iter = 4000, chains = 4, cores = 1,
         control = list(adapt_delta = 0.95, max_treedepth = 20))

post_n <- plyr::ldply(out, function(x) {
  e <- extract(fitModel)
  sm_summ <- summary(fitModel)$summary
  max_rhat <- max(sm_summ[, "Rhat"])
  min_neff <- min(sm_summ[, "n_eff"])
  data.frame(alpha = e$b_j[, 1], beta = e$b_j[, 2], max_rhat, min_neff)
})
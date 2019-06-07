### Practice fitting varying intercept, varying slope models w/ two groups

library(brms)
library(dplyr)
library(ggplot2) 
library(rstan)
library(rethinking)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

focal <- read.csv(here::here('data', 'scaled_revenue_dataset.csv')) %>% 
  arrange(SZ, year) %>% 
  mutate(yrID = as.numeric(as.factor(year)),
         szID = as.numeric(as.factor(SZ))) %>% 
  select(yrID, szID, div = simpson.diversity, cpue, 
         rev = revenue_per_fisher_perday,
         revKg  = revenue_per_fisher_perday_perkg) %>% 
  mutate(divZ = as.numeric(scale(div)),
         revZ = as.numeric(scale(rev)),
         revKgZ = as.numeric(scale(revKg)))


brms::make_stancode(revZ ~ divZ, data = focal)


revMod <- map2stan(
  alist(
    revZ ~ dnorm(mu, sigma),
    #linear models
    mu <- A + BD*divZ,
    #A <- a + a_sz[szID] + a_yr[yrID],
    A <- a  + a_yr[yrID],
    BD <- bd + bd_yr[yrID],
    #adaptive priors
    c(a_yr, bd_yr)[yrID] ~ dmvnorm2(0, phi_yr, rho_yr),
    # a_sz[szID] ~ dnorm(0, phi_sz),
    #fixed priors
    a ~ dnorm(0, 5),
    bd ~ dnorm(0, 5),
    sigma ~ dcauchy(0, 2),
    phi_yr ~ dcauchy(0, 2),
    # c(phi_sz, phi_yr) ~ dcauchy(0, 2),
    rho_yr ~ dlkjcorr(2)
  ), 
  data = focal,
  iter=2000, warmup=1000, chains=1, cores=1
)

# one intercept model
revMod <- map2stan(
  alist(
    revZ ~ dnorm(mu, sigma),
    #linear models
    mu <- A + bd*divZ,
    A <- a + a_sz[szID],
    #adaptive priors
    a_sz[szID] ~ dnorm(0, phi_sz),
    #fixed priors
    a ~ dnorm(0, 5),
    bd ~ dnorm(0, 5),
    sigma ~ dcauchy(0, 2),
    phi_sz ~ dcauchy(0, 2)
  ), 
  data = focal,
  iter=2000, warmup=1000, chains=1, cores=1
)

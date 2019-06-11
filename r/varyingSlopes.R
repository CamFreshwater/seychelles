### Practice fitting varying intercept, varying slope models w/ two groups

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
         rev = revenue,
         revDay = revenue_per_fisher_perday,
         revDayKg  = revenue_per_fisher_perday_perkg) %>% 
  mutate(divZ = as.numeric(scale(div)),
         revZ = as.numeric(scale(rev)),
         revDayZ = as.numeric(scale(revDay)),
         revDayKgZ = as.numeric(scale(revDayKg)))

## Use rethinking to check what the stan structure should look like for 
# slope/intercept covariance
revDayMod1 <- map2stan(
  alist(
    revDay ~ dnorm(mu, sigma),
    #linear models
    mu <- A + BD*divZ,
    A <- a + a_sz[szID] + a_yr[yrID],
    BD <- bd + bd_yr[yrID],
    #adaptive priors
    c(a_yr, bd_yr)[yrID] ~ dmvnormNC(phi_yr, rho_yr),
    a_sz[szID] ~ dnorm(0, phi_sz),
    #fixed priors
    a ~ dnorm(0, 5),
    bd ~ dnorm(0, 5),
    sigma ~ dcauchy(0, 2),
    c(phi_sz, phi_yr) ~ dcauchy(0, 2),
    rho_yr ~ dlkjcorr(2)
  ), 
  data = focal,
  iter=3000, warmup=1000, chains=4, cores=4
)
revDayMod2 <- map2stan(
  alist(
    revDayZ ~ dnorm(mu, sigma),
    #linear models
    mu <- A + BD*divZ,
    A <- a + a_sz[szID] + a_yr[yrID],
    BD <- bd + bd_yr[yrID],
    #adaptive priors
    c(a_yr, bd_yr)[yrID] ~ dmvnormNC(phi_yr, rho_yr),
    a_sz[szID] ~ dnorm(0, phi_sz),
    #fixed priors
    a ~ dnorm(0, 5),
    bd ~ dnorm(0, 5),
    sigma ~ dcauchy(0, 2),
    c(phi_sz, phi_yr) ~ dcauchy(0, 2),
    rho_yr ~ dlkjcorr(2)
  ), 
  data = focal,
  iter=3000, warmup=1000, chains=4, cores=4
)

stancode(revDayMod1)
precis(revDayMod1, depth = 2)


## Fit equivalent model with stan
normMod <- stan_model(here::here("r", "stanModels", 
                                 "varyIntSlope_normal.stan"))

# Prep input data 
N <- nrow(focal) #sample size
nYr <- length(unique(focal$yrID)) #number of yrs 
nSz <- length(unique(focal$szID)) #number of sz groups
yrID <- focal$yrID
szID <- focal$szID

# test normal model
m_hier <- sampling(normMod, data=list(N=N, J=nYr, K=nSz, fisher_id = szID,
                                      yr_id = yrID, x1=focal[, "divZ"],
                                      y = focal[, "revZ"]),
                   iter = 3000, chains = 4, cores = 4,
                   control = list(adapt_delta = 0.9, max_treedepth = 20))
shinystan::launch_shinystan(m_hier)
## converges well

## Fit student-t model with stan
studTMod <- stan_model(here::here("r", "stanModels", 
                                 "varyIntSlope_studT_rand.stan"))


## Fit student-t models to each response in sequence
# This doesn't seem to work, not sure why. Fit separately instead
# respSeq <- list(focal[ ,"rev"], focal[ ,"revDay"], focal[ ,"revDayKg"])
# modOutList <- lapply(respSeq, function(x) {
#   sampling(studTMod, data=list(N=N, J=nYr, K=nSz, fisher_id = szID, 
#                                yr_id = yrID, x1=focal[, "divZ"], 
#                                y = x),
#            iter = 3000, chains = 4, cores = 4, 
#            control = list(adapt_delta = 0.9, max_treedepth = 20))
#   })
# names(modOutList) <- c("rev", "revDay", "revDayKg")

# revMod <- sampling(studTMod, data=list(N=N, J=nYr, K=nSz, fisher_id = szID, 
#                              yr_id = yrID, x1=focal[, "divZ"], 
#                              y = focal[ ,"rev"]),
#                    iter = 3000, chains = 4, cores = 4,
#                    control = list(adapt_delta = 0.9, max_treedepth = 20))
revDayMod <- sampling(studTMod, data=list(N=N, J=nYr, K=nSz, fisher_id = szID,
                                       yr_id = yrID, x1=focal[, "divZ"],
                                       y = focal[ ,"revDay"]),
                   iter = 3000, chains = 4, cores = 4,
                   control = list(adapt_delta = 0.9, max_treedepth = 20))
revDayModRaw <- sampling(studTMod, data=list(N=N, J=nYr, K=nSz, fisher_id = szID, 
                                          yr_id = yrID, x1=focal[, "div"], 
                                          y = focal[ ,"revDay"]),
                      iter = 3000, chains = 4, cores = 4,
                      control = list(adapt_delta = 0.9, max_treedepth = 20))
revDayModZ <- sampling(studTMod, data=list(N=N, J=nYr, K=nSz, fisher_id = szID, 
                                          yr_id = yrID, x1=focal[, "divZ"], 
                                          y = focal[ ,"revDayZ"]),
                      iter = 3000, chains = 4, cores = 4,
                      control = list(adapt_delta = 0.9, max_treedepth = 20))
# revDayKgMod <- sampling(studTMod, data=list(N=N, J=nYr, K=nSz, fisher_id = szID, 
#                                           yr_id = yrID, x1=focal[, "divZ"], 
#                                           y = focal[ ,"revDayKg"]),
#                       iter = 3000, chains = 4, cores = 4,
#                       control = list(adapt_delta = 0.9, max_treedepth = 20))

# modOutList <- list(revMod, revDayMod, revDayKgMod)
# names(modOutList) <- c("rev", "revDay", "revDayKg")
# respSeq <- list("rev", "revDay", "revDayKg")

modOutList <- list(revDayMod, revDayModZ)
names(modOutList) <- c("revDay", "revDayZ")
respSeq <- list("revDay", "revDayZ")

# shinystan::launch_shinystan(revMod)
shinystan::launch_shinystan(revDayMod)
shinystan::launch_shinystan(revDayModZ)
shinystan::launch_shinystan(revDayModRaw)
# shinystan::launch_shinystan(revDayKgMod)
#all look good


# Look at posterior predictions
ppFunc <- function(respName) {
  y_data <- focal[, respName]
  mod <- modOutList[[respName]]
  y_rep <- as.matrix(mod, pars = "y_rep")
  bayesplot::ppc_dens_overlay(y_data, y_rep[1:200, ]) + 
    labs(title = respName) +
    xlim(min(y_data), max(y_data))
}
postPlots <- lapply(respSeq, ppFunc)

postPlots
#some of these look like they could benefit from a skewed distribution 
#(last one fits best)

## Look at coefficient estimates for one model
theta <- extract(modOutList[["revDayKg"]])
d <- data.frame(mu_alpha = theta$mu_a, betaDiv = theta$b1, sigma = theta$sigma,
                harvSig = theta$sigma_a_fisher, yrSigAlpha = theta$phi_yr[,1],
                yrSigBeta = theta$phi_yr[,2], nu = theta$nu) %>% 
  tidyr::gather(key = "par", value = "est")

ggplot(d, aes(as.factor(par), est)) +
  geom_violin() + 
  labs(title = "rev")

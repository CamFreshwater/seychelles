## Fit varying intercept models w/ normal or student-t variance parameters to 
# Seychelles fisheries revenue data

library(dplyr) 
library(ggplot2) 
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

focal <- read.csv(here::here('data', 'scaled_revenue_dataset.csv')) %>% 
  arrange(SZ, year) %>% 
  mutate(yrID = as.numeric(as.factor(year)),
         szID = as.numeric(as.factor(SZ))) %>% 
  select(yrID, szID, div = simpson.diversity, cpue, 
         rev = revenue_per_fisher_perday) %>% 
  mutate(divZ = as.numeric(scale(div)),
         revZ = as.numeric(scale(rev)))
  
# Models
normMod <- stan_model(here::here("r", "varyInt_normal_oneBeta.stan"))
studTMod <- stan_model(here::here("r", "varyInt_studT_oneBeta.stan"))

# Prep input data 
N <- nrow(focal) #sample size
nYr <- length(unique(focal$yrID)) #number of yrs 
nSz <- length(unique(focal$szID)) #number of sz groups
yrID <- focal$yrID
szID <- focal$szID

m_hier <- sampling(normMod, data=list(N=N, J=nYr, K=nSz, fisher_id = szID,
                                   yr_id = yrID, x1=focal[, "divZ"], 
                                   y = focal[, "revZ"]),
                   iter = 4000, chains = 4, cores = 4,
                   control = list(adapt_delta = 0.9, max_treedepth = 20))
shinystan::launch_shinystan(m_hier)


m_hierT <- sampling(studTMod, data=list(N=N, J=nYr, K=nSz, fisher_id = szID,
                                      yr_id = yrID, x1=focal[, "divZ"],
                                      y = focal[, "revZ"]),
                   iter = 4000, chains = 4, cores = 4,
                   control = list(adapt_delta = 0.9, max_treedepth = 20))
shinystan::launch_shinystan(m_hierT)

print(m_hier)
print(m_hierT)

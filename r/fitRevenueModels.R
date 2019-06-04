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
         rev = revenue_per_fisher_perday,
         revKg  = revenue_per_fisher_perday_perkg) %>% 
  mutate(divZ = as.numeric(scale(div)),
         revZ = as.numeric(scale(rev)),
         revKgZ = as.numeric(scale(revKg)))
  
# Models
normMod <- stan_model(here::here("r", "varyInt_normal_oneBeta_rand.stan"))
studTMod <- stan_model(here::here("r", "varyInt_studT_oneBeta_rand.stan"))

# Prep input data 
N <- nrow(focal) #sample size
nYr <- length(unique(focal$yrID)) #number of yrs 
nSz <- length(unique(focal$szID)) #number of sz groups
yrID <- focal$yrID
szID <- focal$szID

# revenue by day models
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

# revenue by day and mass models
m_hier2 <- sampling(normMod, data=list(N=N, J=nYr, K=nSz, fisher_id = szID,
                                      yr_id = yrID, x1=focal[, "divZ"], 
                                      y = focal[, "revKgZ"]),
                   iter = 4000, chains = 4, cores = 4,
                   control = list(adapt_delta = 0.9, max_treedepth = 20))
shinystan::launch_shinystan(m_hier2)

m_hierT2 <- sampling(studTMod, data=list(N=N, J=nYr, K=nSz, fisher_id = szID,
                                        yr_id = yrID, x1=focal[, "divZ"],
                                        y = focal[, "revKgZ"]),
                    iter = 4000, chains = 4, cores = 4,
                    control = list(adapt_delta = 0.9, max_treedepth = 20))
shinystan::launch_shinystan(m_hierT2)

modOutList <- list(m_hier, m_hierT, m_hier2, m_hierT2)
names(modOutList) <- c("revByDayNormal", "revByDayStudT", "revByDayKGNormal", 
                       "revByDayKGStudT")
saveRDS(modOutList, here::here("data", "modFits", "parEstMCMC.rds"))

#function to extract fixed effects
extFunc <- function(modOut, studT = FALSE) {
  e <- extract(modOut)
  if (studT == FALSE) {
    d <- data.frame(mu_alpha = e$mu_a, betaDiv = e$b1, sigma = e$sigma,
                    harvSig = e$sigma_a_fisher, yrSig = e$sigma_a_yr) %>% 
      mutate(model = "norm")
  } else {
    d <- data.frame(mu_alpha = e$mu_a, betaDiv = e$b1, sigma = e$sigma,
                    harvSig = e$sigma_a_fisher, yrSig = e$sigma_a_yr, 
                    nu = e$nu) %>% 
      mutate(model = "studT")
  }
  d %>% 
    tidyr::gather(key = "par", value = "est", -model)
}

normEst1 <- extFunc(m_hier, studT = FALSE)
studTEst1 <- extFunc(m_hierT, studT = TRUE)
post1 <- rbind(normEst1, studTEst1) %>% 
  mutate(response = "revByDay")

normEst2 <- extFunc(m_hier2, studT = FALSE)
studTEst2 <- extFunc(m_hierT2, studT = TRUE)
postOut <- rbind(normEst2, studTEst2) %>% 
  mutate(response = "revByDayKg") %>% 
  rbind(., post1) %>% 
  mutate(par = factor(par, levels = c("mu_alpha", "betaDiv", "nu", "sigma",
                                      "harvSig", "yrSig")))

## parameter estimates seem similar, but the fact that nu is small (with high
# confidence), suggests that there is strong evidence for heavy tails; similar 
# patterns with both response variables as well
png(here("figs", "parEstimates.png"), height = 6, 
    width = 6, units = "in", res = 300)
ggplot(postOut %>% filter(response == "revByDay"), 
       aes(as.factor(par), est)) +
  geom_violin() + 
  facet_wrap(~model) + 
  labs(title = "revByDay")
ggplot(postOut %>% filter(response == "revByDayKg"), 
       aes(as.factor(par), est)) +
  geom_violin() + 
  facet_wrap(~model) + 
  labs(title = "revByDayKg")
dev.off()


## look at posterior predictive checks
y_data <- focal[, "revZ"]
y_rep <- as.matrix(m_hierT, pars = "y_rep")
dim(y_rep)

bayesplot::ppc_dens_overlay(y_data, y_rep[1:200, ]) + xlim(-2, 2)

#given similar par estimates focus on revByDay, studentT models for now
mcmcOut <- extract(modOutList[["revByDayStudT"]])
str(mcmcOut)
range(focal$divZ)
preds <- cbind(mcmcOut$mu_a, mcmcOut$b1)

X_new <- model.matrix(~x, data = data.frame(x = seq(-2, 2, by = 0.2)))
pred_X <- apply(preds, 1, function(estBeta) X_new %*% estBeta) %>% 
  apply(., 1, quantile, probs = c(0.025, 0.5, 0.975))

plot(focal[, "revZ"] ~ focal[, "divZ"], pch=1, xlab="Diversity", 
     ylab="Revenue")
lines(seq(-2,2, by=0.2), pred_X[1,], lty=2, col="red")
lines(seq(-2,2, by=0.2), pred_X[2,], lty=1, lwd=3, col="blue")
lines(seq(-2,2, by=0.2), pred_X[3,], lty=2, col="red")


plot(log(focal[, "revZ"]) ~ focal[, "divZ"], pch=1, xlab="Diversity", 
     ylab="Revenue")



# coefficient plot
#now we could look at the variation in the regression coefficients between the 
#groups doing caterpillar plots
yr_coeff <- apply(mcmc_hier$beta, c(2,3), quantile, probs = c(0.025, 0.5, 0.975))
df_ind_coeff <- data.frame(Coeff = rep(c("(Int)"), each=10),
                           LI = c(ind_coeff[1, , 1]), 
                           Median = c(ind_coeff[2, , 1], ind_coeff[2, , 2],
                                      ind_coeff[2,,3]),
                           HI=c(ind_coeff[3, , 1], ind_coeff[3, , 2], 
                                ind_coeff[3, , 3]))
gr <- paste("Gr", 1:10)
df_ind_coeff$Group <- factor(gr, levels=gr)
#we may also add the population-level median estimate
pop_lvl <- data.frame(Coeff = c("(Int)", "X1", "X2"), 
                      Median = apply(mcmc_hier$gamma, 2, quantile, probs = 0.5))
ggplot(df_ind_coeff, aes(x = Group, y = Median)) + 
  geom_point() +
  geom_linerange(aes(ymin = LI, ymax = HI)) +
  coord_flip() +
  facet_grid(.~Coeff) +
  geom_hline(data = pop_lvl, aes(yintercept = Median), color="blue", 
             linetype = "dashed") +
  labs(y = "Regression parameters")
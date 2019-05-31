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
            # X_ij = X_ij,
            y_i = df1$y2)

library(lme4)
lmer(y ~ (1 | cat) + x, df1) %>% coef()
# lmer(y2 ~ (1 | cat), df1) %>% coef()

# model
mod <- stan_model(here::here("r/exampleStanModels/varyInt_lm.stan"))

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

#-----

# A slightly different but equivalent model with varying slopes as well
mod2 <- stan_model(here::here("r/exampleStanModels/varyInt_lm_Full.stan"))

#simulate some data
set.seed(20161110)
N<-100 #sample size
J<-10 #number of plant species
id<-rep(1:J,each=10) #index of plant species
K<-3 #number of regression coefficients
#population-level regression coefficient
gamma<-c(2,-1,3)
#standard deviation of the group-level coefficient
tau<-c(0.3,2,1)
#standard deviation of individual observations
sigma<-1
#group-level regression coefficients
beta<-mapply(function(g,t) rnorm(J,g,t),g=gamma,t=tau) 
#the model matrix
X<-model.matrix(~x1 + x2,data = data.frame(x1 = runif(N,-2,2),
                                     x2 = runif(N,-2,2)))
y<-vector(length = N)
for(n in 1:N){
  #simulate response data
  y[n]<-rnorm(1, X[n,] %*% beta[id[n], ], sigma)
}

m_hier <- sampling(mod2, data=list(N=N, J=J, K=K, id=id, X_ij=X, y_i=y),
                   iter = 4000, chains = 4, cores = 1,
                   control = list(adapt_delta = 0.95, max_treedepth = 20))

## results in divergent iterations - reparameterize
mod3 <- stan_model(here::here("r/exampleStanModels/varyInt_lm_Full_repar.stan"))

m_hier2 <- sampling(mod3, data=list(N=N, J=J, K=K, id=id, X_ij=X, y_i=y),
                   iter = 4000, chains = 4, cores = 1,
                   control = list(adapt_delta = 0.95, max_treedepth = 20))
#estimates of hyper parameters
print(m_hier2)

## As above but with ML
df <- cbind(id, X) %>% 
  as.data.frame() %>% 
  rename(int = "(Intercept)") %>% 
  cbind(y, .)

ml2 <- lmer(y ~ (x1 + x2 + 1) | id, df)

ml2 %>% 
  summary()
apply(mlOut$id, 2, mean)

## Revisit coefficient specific estimates
mcmc_hier <- extract(m_hier)
str(mcmc_hier)

#plot average response to explanatory variables
X_new <- model.matrix(~x+y,data = data.frame(x = seq(-2, 2, by = 0.2),
                                             y = 0))
#get predicted values for each MCMC sample
pred_x1 <- apply(mcmc_hier$gamma, 1, function(estBeta) X_new %*% estBeta)
#now get median and 95% credible intervals
pred_x1 <- apply(pred_x1, 1, quantile, probs=c(0.025, 0.5, 0.975))
#same stuff for the second explanatory variables
X_new <- model.matrix(~x+y,data = data.frame(x = 0, 
                                             y= seq(-2, 2, by=0.2)))
pred_x2 <- apply(mcmc_hier$gamma, 1, function(estBeta) X_new %*% estBeta)
pred_x2 <- apply(pred_x2, 1, quantile, probs=c(0.025, 0.5, 0.975))

# coefficient plot
#now we could look at the variation in the regression coefficients between the groups doing caterpillar plots
ind_coeff <- apply(mcmc_hier$beta, c(2,3), quantile, probs = c(0.025, 0.5, 0.975))
df_ind_coeff <- data.frame(Coeff = rep(c("(Int)", "X1", "X2"), each=10),
                           LI = c(ind_coeff[1, , 1],ind_coeff[1, , 2], 
                                  ind_coeff[1, ,3]), 
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



#-----

# As above but with fixed slopes
mod4 <- stan_model(here::here("r/exampleStanModels", 
                              "varyInt_lm_Full_repar_noCoefMatrix.stan"))
mod5 <- stan_model(here::here("r/exampleStanModels", 
                              "varyInt_lm_studT.stan"))

#simulate some data
set.seed(20161110)
N<-100 #sample size
J<-10 #number of plant species
id<-rep(1:J,each=10) #index of plant species
#population-level regression coefficient
mu_alpha <- 2 
#standard deviation of the group-level coefficient
sigma_a <- 0.5
b1 <- -1
b2 <- 3
  
#standard deviation of individual observations
sigma <- 1
#group-level regression coefficients
beta <- mapply(function(g,t) rnorm(J,g,t), g=mu_alpha, t=sigma_a) %>% 
  cbind(rep(b1, J),
        rep(b2, J))

#the model matrix
X <- model.matrix(~x1 + x2, data = data.frame(x1 = runif(N, -2, 2),
                                            x2 = runif(N, -2, 2)))
y <- vector(length = N)
for(n in 1:N){
  #simulate response data
  y[n]<-rnorm(1, X[n, ] %*% beta[id[n], ], sigma)
}

m_hier <- sampling(mod4, data=list(N=N, J=J, group_id=id, x1=X[,"x1"], 
                                   x2=X[,"x2"], y=y),
                   iter = 4000, chains = 4, cores = 1,
                   control = list(adapt_delta = 0.9, max_treedepth = 20))
m_hierT <- sampling(mod5, data=list(N=N, J=J, group_id=id, x1=X[,"x1"], 
                                   x2=X[,"x2"], y=y),
                   iter = 4000, chains = 1, cores = 1,
                   control = list(adapt_delta = 0.9, max_treedepth = 20))

#estimates of hyper parameters
print(m_hier)
print(m_hierT)


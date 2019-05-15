library(here) 
library(dplyr) 
library(ggplot2) 
library(funk) 
library(rethinking)


focal<-read.csv('scaled_revenue_dataset.csv')

qplot(simpson.diversity,revenue_per_fisher_perday_perkg, data = focal)
qplot(simpson.diversity,revenue_per_fisher_perday, data = focal)

### Is daily revenue, per fisher predicted by catch diversity and catch size?
m2 <- map2stan(
	alist(
	    revenue_per_fisher_perday ~ dnorm( mu , sigma ) ,
	    mu <- a + ar[SZ] + ar2[year] + 
	    			   bD*simpson.diversity +
	    			   bE*cpue ,
	    a ~ dnorm(0, 5),
	    c(ar)[SZ] ~  dnorm(0,sigmar1), 
	    c(ar2)[year] ~  dnorm(0,sigmar2),
  	    c(bD, bE) ~ dnorm(0, 5),
	    c(sigma, sigmar1, sigmar2) ~ dcauchy( 0 , 2)
), data=focal, iter=3000, chains=3)
precis(m2)



### Is daily revenue, per fisher and corrected for catch size, predicted by catch diversity?
m3 <- map2stan(
	alist(
	    revenue_per_fisher_perday_perkg ~ dnorm( mu , sigma ) ,
	    mu <- a + ar[SZ] + ar2[year] + 
	    			   bD*simpson.diversity 	
	    			   ,
	    a ~ dnorm(0, 5),
	    c(ar)[SZ] ~  dnorm(0,sigmar1), 
	    c(ar2)[year] ~  dnorm(0,sigmar2),
  	    c(bD) ~ dnorm(0, 5),
	    c(sigma, sigmar1, sigmar2) ~ dcauchy( 0 , 2)
), data=focal, iter=3000, chains=3)
precis(m3)




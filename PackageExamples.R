set.seed(1)
######################################
##  R file to run examples for the
## cptmcmc R package
######################################
library(cptmcmc)

# The package tries to imitate our eyes.
# Finds the first instance of a changepoint
# by looking for the first time-index when the 
# running mean estimate crosses over the overall mean (method = 1)
# or
# the running average crosses over the average beyond
# the index of the running average (method = 2)(default) 


############
## Univariate MH for a N(0,1)
# good jump-size with h = 1 and bad starting value
chain <- MH_uni(n = 1e3, start = 10, h = 1, mean = 0, var = 1)

# this gives a plot by default
burned_chain <- cptmcmc(chain)


# bad jump-size with h = .1 and bad starting value
chain <- MH_uni(n = 1e3, start = 10, h = .1, mean = 0, var = 1)
burned_chain <- cptmcmc(chain)


## 5-dim normal target with independent covariance matrix
# in the multivariate case, the univariate rule is applied
# to each component and the largest burn-in time is chosen
chain <- MH_multi(n = 1e4, k = 5, h = rep(.01, 5), mean = rep(0,5), sigma = diag(1,5), start = c(10, -10, 10, -10, 10))

# this produced  two plots
# (1) the trace plot with burn-in time for the component that has the largest burn-in time
# (2) a boxplot of all the burn-in time for all components
# (yes, the color of the boxplot is yellow. Kids!).
burned_chain <- cptmcmc(chain)

# plotting the individual trace plots with the
# determined burn-in time
par(mfrow = c(3,2))
for(i in 1:5)
{
	plot(chain[,i], type = 'l')
	abline(v = max(burned_chain$breakpt), col = "red")
}


## real data example
# multinomial regression on the Nethvote
# dataset in library MCMCpack. 22-dim chain.
# ignore the warnings()
chain <- multinomial_regression(n = 1e4)
burned_chain <- cptmcmc(chain)

# change par if you have large screen to c(5,5)
# notice that each individual burn-in times would've been
# smaller, but choosing the maximum time over here helps us
par(mfrow = c(3,4))
for(i in 1:dim(chain)[2])
{
	plot(chain[,i], type = 'l')
	abline(v = max(burned_chain$breakpt), col = "red")
}




#### Examples to sample chains

#' @title Multinomial Regression examples
#' @description Generate multinomial regression data
#' @name multinomial_regression
#' @param n : No. of observations required
#'
#' @return
#' Returns a matrix with dimensions n*22.
#'
#' @examples
#'\dontrun{
#' multinomial_regression(1e4)}

#' @rdname multinomial_regression
#' @export
multinomial_regression <- function(n)
{

  running_multinomial <- function(rep, m, n,p) {
    data("Nethvote")
    a <- floor(sqrt(n))
    b <- a
    # stores Markov chains
    X <- array(0, dim = c(m, n, p))
    random <- matrix(sample(1:(rep * m), rep * m), nrow = rep)

    for (j in 1:rep) {
        MLE <- MCMCmnl(vote ~ choicevar(distD66, "sqdist", "D66") +
                     choicevar(distPvdA, "sqdist", "PvdA") +
                     choicevar(distVVD, "sqdist", "VVD") +
                     choicevar(distCDA, "sqdist", "CDA") +
                     relig + class + income + educ + age + urban,
                   baseline = "D66", mcmc.method = "RWM", B0 = 0,
                   verbose = FALSE, mcmc = 1, thin = 1, tune = 0.5, burnin = 0, seed = j,
                   data = Nethvote)

      starting <- matrix(0, nrow = m, ncol = p)
      for (i in 1:p) {
        starting[, i] <- seq(MLE[1, i] - 2, MLE[1, i] + 2, length = m)
     }

    # sample Markov chains
    for (i in 1:m) {

      X[i,,] <- MCMCmnl(vote ~ choicevar(distD66, "sqdist", "D66") +
                          choicevar(distPvdA, "sqdist", "PvdA") +
                          choicevar(distVVD, "sqdist", "VVD") +
                          choicevar(distCDA, "sqdist", "CDA") +
                          relig + class + income + educ + age + urban,
                        baseline = "D66", beta.start = starting[i,], mcmc.method = "RWM", B0 = 0,
                        verbose = 0, mcmc = n, thin = 1, tune = .5, burnin = 0, seed = i * j,
                        data = Nethvote)
     }

   }
    return(X)
  }

  rep <- 1
  m <- 2
  p <- 22
  x <- running_multinomial(rep, m, n,p)
  return(x[1,,])
}

#' @title MH for multivariate normal
#' @description Generate Markov Chain for target = Nk(mean, sigma), proposal = Nk(x,diag(h))
#' @name MH_multi
#' @param n :no. of observations
#' @param start : vector of starting values
#' @param k : no of dimensions i.e k-variate normal(mean,sigma) is the target
#' @param h : vector of step size for MH, i.e proposal will have covariance matrix diag(h)
#' @param mean : true target mean
#' @param sigma : true target covariance
#'
#' @return
#' Returns a markov chain matrix with dimensions n*k
#'
#' @examples
#'\dontrun{
#'MH_multi(n=1e4, start=c(5,5,5), k=3, h =c(0.1,0.1,0.01), mean=c(0,0,0), sigma=diag(3))
#'}




### target = Nk(mean, sigma), proposal = Nk(x,diag(h))
#' @rdname MH_multi
#' @export
MH_multi <- function(n, start, k, h, mean, sigma)
{
  chain = matrix(0,n,k)
  chain[1,] = start

  for(i in 2:n)
  {
    accept = 0
    proposal = rmvnorm(1, mean = chain[i-1,], sigma=diag(h))
    alpha = exp(log(dmvnorm(proposal,mean,sigma)) - log(dmvnorm(chain[i-1,],mean,sigma)))
    u = runif(1)
    if(u < alpha)
      chain[i,] = proposal
    else
    {
      chain[i,] = chain[i-1,]
      accept = accept + 1
    }
  }
  return(chain)
}

#' @title MH for univariate normal
#' @description Generate Markov Chain for target = N(mean, var), proposal = N(x,h)
#' @name MH_multi
#' @param n :no. of observations
#' @param start : starting value
#' @param h : var for proposal
#' @param mean : true target mean
#' @param var : true target covariance
#'
#' @return
#' Returns a markov chain vector
#'
#' @examples
#'\dontrun{
#'MH_uni(n=1e4, start=5, h=0.1, mean=0, var=1)
#'}



### target = N(mean, var), proposal = N(x,h)
#' @rdname MH_uni
#' @export
MH_uni <- function(n, start, h, mean, var)
{
  chain = numeric(n)
  chain[1] = start

  for(i in 2:n)
  {
    accept = 0
    proposal = rnorm(1,chain[i-1],h)
    alpha = exp(log(dnorm(proposal,mean,sqrt(var))) - log(dnorm(chain[i-1],mean,sqrt(var))))
    u = runif(1)
    if(u < alpha)
      chain[i] = proposal
    else
    {
      chain[i] = chain[i-1]
      accept = accept + 1
    }
  }
  return(chain)
}

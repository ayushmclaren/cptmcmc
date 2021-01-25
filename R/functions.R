###################################################################
#### Function for component wise (univariate) detection using Local Max

cpt_uni <- function(chain)
{
  n = length(chain)
  # defining the cusum process
  z = (cumsum(chain) - (1:n)*mean(chain))/sqrt(n)
  localmax_idx = 0

  for(i in 1:n)
  {
    if(abs(z[i])>abs(z[i+1]))
    {
      localmax_idx = i+1
      break
    }
  }
  #cat('\nLocalMax method:\nThe estimated break pt. is at:', localmax_idx, sep = "")
  return(localmax_idx)
}

###################################################################
#### Function for component wise (multivariate) detection using Local Max
cpt_multi <- function(chain)
{
  breakpts = numeric(ncol(chain))
  breakpts = apply(chain,2,cpt_uni)
  return(breakpts)
}

###################################################################
#### Function for component wise (univariate) detection using Modified Local Max method
##cpt_uni_modf <- function(chain)
##{
##  n = length(chain)
##  z = numeric(n)
##  z[1] = 0
##  localmax_idx = 0

##  for(i in 1:(n-1))
##  {
##      z[i+1] = z[i] + (chain[i] - mean(chain[(i+1):n]))/sqrt(n)

##      if(abs(z[i])>abs(z[i+1]))
##      {
##        localmax_idx = i+1
##        break
##      }
##  }
  #cat('\nModified LocalMax method:\nThe estimated break pt. is at:', localmax_idx, sep = "")
##  return(localmax_idx)
##}


###################################################################
#### Function for component wise (multivariate) detection using Modified Local Max
#cpt_multi_modf <- function(chain)
#{
 # breakpts = numeric(ncol(chain))
#  breakpts = apply(chain,2,cpt_uni_modf)
 # return(breakpts)
#}

## usethis namespace: start
#' @useDynLib cptmcmc, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
cpt_uni_modf <- function(chain) {
  .Call(`_cptmcmc_cpt_uni_modf`, chain)
}

## usethis namespace: start
#' @useDynLib cptmcmc, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
cpt_multi_modf <- function(chain) {
  .Call(`_cptmcmc_cpt_multi_modf`, chain)
}

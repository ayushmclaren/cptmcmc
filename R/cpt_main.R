######################################3
### method = 1 is local max, method =2 is modified local max from functions.R
#### for multiple chains run in a loop, each column is treated like a component

#' @title Changepoint for MCMC
#' @description The proposed changepoint algorithms as in the document.
#' @name cptmcmc
#' @param chain : matrix or vector of the Markov chain. For matrix each column is treated as a component.
#' @param method : '1' for local-max method in section 3.1 and 3.2, '2' for modifies as in section 3.3
#' @param logpost : F by default, should be T if passing logposterior vector of the chain
#' @param plot : T by default
#' @param thresh : (0,1). Threshold of maximum viable breakpt fraction.
#'
#' @return
#' Returns a data frame with all the estimatewd breakpts and the truncated chain.
#'
#' @examples
#'\dontrun{
#'cptmcmc(chain_multi,method=1,plot=TRUE)
#'cptmcmc(chain_log,method=2,logpost=TRUE)
#'}


#' @rdname cptmcmc
#' @export
cptmcmc <- function(chain, method = 2, logpost = FALSE, plot = TRUE, thresh = 0.50)
{
  on.exit(par(mfrow=c(1,1)))
  if(method == 1)
  {
    if(logpost)
    {
      breakpt = cpt_uni(chain)
      cat('\nLocalMax method:\nThe estimated break pt. is at:', breakpt,"\n", sep = "")
      if(breakpt > length(chain)*thresh)
        cat("\nThe change pt. detected is at >", thresh*100,
            "%. We recommend re-running the chain with different parameters.\n",sep="")
      if(plot)
      {
        plot.ts(chain, main="Log-Posterior Chain", ylab="chain", xlab="Time")
        abline(v=breakpt, col="red")
      }
      return(list(breakpt = breakpt, chain_truncated = chain[breakpt:length(chain)]))
    }

    if(is.vector(chain))
    {
      breakpt = cpt_uni(chain)
      cat('\nLocalMax method:\nThe estimated break pt. is at:', breakpt,"\n", sep = "")
      if(breakpt > length(chain)*thresh)
        cat("\nThe change pt. detected is at >", thresh*100,
            "%. We recommend re-running the chain with different parameters.\n",sep="")
      if(plot)
      {
        plot.ts(chain, main="Markov Chain", ylab="chain", xlab="Time")
        abline(v=breakpt, col="red")
      }
      return(list(breakpt = breakpt, chain_truncated = chain[breakpt:length(chain)]))
    }

    else
    {
      breaks = cpt_multi(chain)
      component = which.max(breaks)
      breakpt = breaks[which.max(breaks)]
      cat('\nLocalMax method:\nThe estimated break pt. is at:', breakpt,
          "\nThe component with Max break pt. is:", component, "\n", sep = "")
      if(breakpt > nrow(chain)*thresh)
        cat("\nThe change pt. detected is at >", thresh*100,
            "%. We recommend re-running the chain with different parameters.\n",sep="")
      if(plot)
      {
        par(mfrow=c(1,2))
        plot(chain[,component], main=paste0("Markov Chain Component ",component),
                ylab="chain", xlab="Time",type='l')
        abline(v=breakpt, col="red")
        boxplot(breaks, main="Estimated Change Points", xlab = "Estimates",
                col = "orange", border = "brown", horizontal = TRUE)
      }
      return(list(breakpt = breaks, chain_truncated = chain[breakpt:nrow(chain),]))
    }

  }

  if(method == 2)
  {
    if(logpost)
    {
      breakpt = cpt_uni_modf(chain)
      cat('\nModified LocalMax method:\nThe estimated break pt. is at:', breakpt,"\n", sep = "")
      if(breakpt > length(chain)*thresh)
        cat("\nThe change pt. detected is at >", thresh*100,
            "%. We recommend re-running the chain with different parameters.\n",sep="")
      if(plot)
      {
        plot.ts(chain, main="Log-Posterior Chain", ylab="chain", xlab="Time")
        abline(v=breakpt, col="red")
      }
      return(list(breakpt = breakpt, chain_truncated = chain[breakpt:length(chain)]))
    }

    if(is.vector(chain))
    {
      breakpt = cpt_uni_modf(chain)
      cat('\nModified LocalMax method:\nThe estimated break pt. is at:', breakpt,"\n", sep = "")
      if(breakpt > length(chain)*thresh)
        cat("\nThe change pt. detected is at >", thresh*100,
            "%. We recommend re-running the chain with different parameters.\n",sep="")
      if(plot)
      {
        plot.ts(chain, main="Markov Chain", ylab="chain", xlab="Time")
        abline(v=breakpt, col="red")
      }
      return(list(breakpt = breakpt, chain_truncated = chain[breakpt:length(chain)]))
    }

    else
    {
      breaks = cpt_multi_modf(chain)
      component = which.max(breaks)
      breakpt = breaks[which.max(breaks)]
      cat('\nModified LocalMax method:\nThe estimated break pt. is at:', breakpt,
          "\nThe component with Max break pt. is:", component, "\n", sep = "")
      if(breakpt > length(chain)*thresh)
        cat("\nThe change pt. detected is at >", thresh*100,
            "%. We recommend re-running the chain with different parameters.\n", sep="")
      if(plot)
      {
        par(mfrow=c(1,2))
        plot(chain[,component], main= paste0("Markov Chain Component ",component),
                ylab="chain", xlab="Time",type='l')
        abline(v=breakpt, col="red")
        boxplot(breaks, main="Estimated Change Points", xlab = "Estimates",
                col = "orange", border = "brown", horizontal = TRUE)
      }
      return(list(breakpt = breaks, chain_truncated = chain[breakpt:nrow(chain),]))
    }

  }
  stop("Specified method should be either '1' or '2'.\n")
}

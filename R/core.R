#' The main mapper function
#' @description
#' Maps a pair of mean and c-statistic value to the parameters of a specified distribution for risk
#'
#' @param target A vector of size 2. The first element is mean and the second element is c-statistic.
#' @param type One of "beta", "logitnorm", "probitnorm". Loosy matching is enabled (so "b" will be mapped to "beta").
#' @return An object of class mcmapper. The "value" component returns the parameter. Any warning or error from the integration or gradient ascent will also be returned in the "info" component.
#' @examples
#' mcmap(c(0.1, 0.75), "beta")
#'@export
mcmap <- function(target, type=c("beta", "logitnorm", "probitnorm"))
{
  type <- match.arg(type)
  if(!is.vector(target) | length(target)!=2)
  {
    stop("Error: target should be a vector oftwo elements (m, c)")
  }
  if(is.null(names(target)))
  {
    m <- target[1]
    c <- target[2]
  }
  else
  {
    m <- target['m']
    c <- target['c']
  }

  if(m<=0 | m>=1 | c<=0.5 | c>=1)
  {
    stop("Error: inadmissible value(s) for m and/or c.")
  }

  if(abs(m-0.5)>0.499  | abs(c-0.75)>0.24)
  {
    message("Extreme values for parameters are requested. The algorithms might not be reliable. Check the results (e.g., via simulation).")
  }

  out <- list(type=type)

  out$value <- switch(
    type,

    "beta"=mcmap_beta(c(m,c)),

    "logitnorm"=
      {
        if(abs(m-0.5)>0.49  | abs(c-0.75)>0.19)
        {
          message("Invoking 1D version of logitnorm mapper given the extrema parameter values.")
          mcmap_logitnorm(c(m,c), method="meansolve")
        }
        else
        {
          mcmap_logitnorm(c(m,c))
        }
      },

    "probitnorm"=mcmap_probitnorm(c(m,c)),

    {stop("Type is not recognized.")}
  )

  class(out)<-"mcmapper_output"

  out
}





#' A generic mapper function
#' @description
#' Maps a pair of mean and c-statistic value to the parameters of an unspecified distribution that is indexed by two parameters
#'
#' @param target A vector of size 2. The first element is mean and the second element is c-statistic.
#' @param CDF Cumulative distribution function of an unnspecified distribution. The CDF must be indexed by two parameters.
#' @param integrate_controls (optional): parameters to be passed to integrate()
#' @param optim_controls (optional): parameters to be passed to optim()
#' @return A vector of size two that contains the distribution parameters
#' @examples
#' mcmap_generic(c(0.1, 0.75), pbeta)
#'@export
mcmap_generic <- function(target, CDF, integrate_controls=list(), optim_controls=list())
{
  m <- target[1]
  c <- target[2]
  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  if(is.null(integrate_controls$lower)) integrate_controls$lower<-0
  if(is.null(integrate_controls$upper)) integrate_controls$upper<-1
  integrate_controls$f <- CDF

  integrate_controls2 <- integrate_controls
  integrate_controls2$f <- function(x, par1, par2) {CDF(x,par1,par2)^2}

  f <- function(x)
  {
    f1 <- do.call(stats::integrate, args=c(integrate_controls, x[1], x[2]))$value
    f2 <- do.call(stats::integrate, args=c(integrate_controls2,x[1], x[2]))$value
    (f1-F1)^2+(f2-F2)^2
  }

  if(is.null(optim_controls$par)) optim_controls$par<-c(0.5,1)
  optim_controls$fn <- f
  res <- do.call(stats::optim, args=optim_controls) #, control=list(trace=100))

  if(res$convergence==0)
    unname(res$par)
  else
    NULL
}












##'@export
# plot.mcmapper_output <- function(mcmapper_output, CDF=F, bins=1000, ...)
# {
#   x <- (0:bins)/bins
#
#   strFun <- paste0(ifelse(CDF,"p","d"), mcmapper_output$type)
#
#   tmp <- as.list(c(NA,unname(mcmapper_output$value)))
#   tmp[[1]] <- x
#
#   y <- do.call(strFun,args=tmp)
#
#   plot(x,y,...)
# }

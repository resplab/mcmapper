#' Mapper function for beta distribution
#' @description
#' Maps a pair of mean and c-statistic value to the parameters of a beta distribution
#'
#' @param target A vector of size 2. The first element is mean and the second element is c-statistic.
#' @param method Not implemented for this funciton yet; leave as empty string.
#' @param integrate_controls (optional): parameters to be passed to integrate()
#' @param optim_controls (optional): parameters to be passed to optim()
#' @return A vector of size two that contains the distribution parameters
#' @examples
#' mcmap_beta(c(0.1, 0.75))
#' @export
mcmap_beta <- function(target, method="", integrate_controls=list(), optim_controls=list())
{
  if(method=="")
  {
    return(mcmap_beta_default(target, integrate_controls, optim_controls))
  }
  else
  {
    stop("The requested method is not supplied.")
  }
}






mcmap_beta_default <- function(target=c(m=0.25,c=0.75), integrate_controls=list(), optim_controls=list())
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- mcmap_beta_default(c(1-m,c), integrate_controls, optim_controls)
    return(c(tmp[2],tmp[1]))
  }

  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  if(is.null(integrate_controls$lower)) integrate_controls$lower<-0
  if(is.null(integrate_controls$upper)) integrate_controls$upper<-1
  integrate_controls$f <- function(x, alpha) {stats::pbeta(x, alpha, alpha*(1-m)/m)^2}

  f <- function(x)
  {
    f2 <- do.call(stats::integrate, args=c(integrate_controls, x))$value
    (f2-F2)^2
  }

  if(is.null(optim_controls$par)) optim_controls$par <- 1 #sigma
  if(is.null(optim_controls$method)) optim_controls$method <- "Brent"
  if(is.null(optim_controls$lower)) optim_controls$lower <- 0.0001
  if(is.null(optim_controls$upper)) optim_controls$upper <- 10000
  optim_controls$fn <- f
  res <- do.call(stats::optim, args=optim_controls)

  if(res$convergence==0)
    c(alpha=unname(res$par[1]), beta=unname(res$par[1]*(1-m)/m))
  else
    NULL
}








#Log-transformed alpha and beta
mcmap_beta2 <- function(target=c(m=0.25,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- mcmap_beta(c(1-m,c), init)
    return(c(tmp[2],tmp[1]))
  }

  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    L <- 0
    U <- 1
    f2 <- stats::integrate(function(x, ln_alpha) {stats::pbeta(x, exp(ln_alpha), exp(ln_alpha)*(1-m)/m)^2}, L, U, ln_alpha=x)$value

    (f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(0)
  }
  res <- stats::optim(init, f, method="Brent", lower=-c(10), upper=c(10)) #, control=list(trace=100))

  if(res$convergence==0)
    c(alpha=exp(res$par[1]), beta=exp(res$par[1])*(1-m)/m)
  else
    NULL
}

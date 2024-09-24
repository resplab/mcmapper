logit <- function(x) {log(x/(1-x))}
expit <- function(x) {1/(1+exp(-x))}


#' Functions related to logit-normal distribution.
#' @description
#' Functions related to logit-normal distribution.
#' @name logitnorm
#' @param n Number of draws requested (for rlogitnorm)
#' @param x For density, CDF, and quantile functions
#' @param mu Mean of the logit-transformed variable
#' @param sigma SD of the logit-transformed variable
#' @return Depends on the function
#' @export
rlogitnorm <- function(n, mu, sigma) {1/(1+exp(-(stats::rnorm(n,mu,sigma))))}
#' @rdname logitnorm
#'@export
dlogitnorm <- function(x, mu, sigma) {1/sigma/sqrt(2*pi)/x/(1-x)*exp(-(logit(x)-mu)^2/(2*sigma^2))}
#' @rdname logitnorm
#'@export
plogitnorm <- function(x, mu, sigma) {stats::pnorm(logit(x),mu,sigma)}
#' @rdname logitnorm
#'@export
qlogitnorm <- function(x, mu, sigma) {1/(1+exp(-(stats::qnorm(x,mu,sigma))))}



elogitnorm_base <- function(mu, sigma, N=100)
{
  A <- B <- C <- 0
  n <- 1:N
  A <- exp(-sigma^2*n^2/2)*sinh(n*mu)*tanh(n*sigma^2/2)
  B <- exp(-(2*n-1)^2*pi^2/2/sigma^2)*sin((2*n-1)*pi*mu/sigma^2)/sinh((2*n-1)*pi^2/sigma^2)
  C <- exp(-sigma^2*n^2/2)*cosh(n*mu)

  bads <- c(which(is.nan(A)), which(is.nan(B)), which(is.nan(C)))
  if(length(bads)>0)
  {
    lg <- min(bads)-1
    A <- A[1:lg]
    B <- B[1:lg]
    C <- C[1:lg]
  }

  1/2+(sum(A)+2*pi/sigma^2*sum(B))/(1+2*sum(C))
}


elogitnorm_recursive <- function(mu, sigma)
{
  if(mu>0)
  {
    return(1-elogitnorm_recursive(-mu,sigma))
  }
  if(abs(-mu-sigma^2)<abs(mu))
  {
    return(exp(mu+sigma^2/2)*elogitnorm_recursive(-mu-sigma^2,sigma))
  }

  elogitnorm_base(mu, sigma)
}


#standard
elogitnorm_nonrecursive <- function(mu, sigma)
{
  if(mu==0) return(0.5)
  if(mu>0)
  {
    return(1-elogitnorm_nonrecursive(-mu,sigma)) #lol this one is just once recursive
  }

  steps <- floor(-mu/sigma^2)

  if(steps==0) {return(elogitnorm_base(mu,sigma))}

  if(steps>1000) {warning("Too many steps when evaluating (mu,sigma): ",mu,sigma)}

  S <- 0;
  P <- 1
  g <- exp(mu+sigma^2/2)
  I <- 1
  for(i in 1:steps)
  {
    P <- P*g
    S <- S+I*P
    g <- g*exp(sigma^2)
    I <- -I
  }
  S <- S + I*P*elogitnorm_base(mu+steps*sigma^2,sigma)

  S
}



elogitnorm <- function(mu,sigma)
{
  #message(mu,",",sigma)
  steps <- floor(-mu/sigma^2)
  if(steps>1000)
  {
    return(1-stats::integrate(plogitnorm,0,1,mu=mu,sigma=sigma)$value)
  }
  else
  {
    return(elogitnorm_nonrecursive(mu,sigma))
  }
}





find_logitnorm_modes <- function(mu, sigma)
{
  f <- function(x) {log(x/(1-x))-(2*x-1)*sigma^2-mu}

  delta <- 1-4/(2*sigma^2)

  if(delta<0)
  {
    res <- stats::uniroot(f,interval=c(0,1))$root
  }
  else
  {
    y <- c(1-sqrt(delta),1+sqrt(delta))/2
    res <- c()
    y <- c(0,y[which(y>0 & y<1)],1)

    for(i in 2:length(y))
    {
      if(sign(f(y[i-1]))!=sign(f(y[i])))
      {
        res <- c(res, stats::uniroot(f,interval=c(y[i-1],y[i]))$root)
      }
    }
  }

  res
}




#Needs to be vectorized on sigma as integrate() will ask for it
find_logitnorm_mu <- function(m,sigma)
{
  if(m==0.5) return(0)
  if(m>0.5) return(-find_logitnorm_mu(1-m,sigma))

  out <- rep(NA,length(sigma))

  for(i in 1:length(sigma))
  {
    #the actual solution is always between logit(m),0. and remember: m<0.5 so mu<0
    L <- logit(m)-sigma[i]
    repeat
    {
      if(elogitnorm(L,sigma[i])<m) break;
      L <- L-sigma[i]
    }

    interval <- c(L,logit(m))
    out[i] <- stats::uniroot(function(x,s) {elogitnorm(x,s)-m}, interval, sigma[i])$root
  }
  out
}




########################################mcmappers


#' Mapper function for logit-normal distribution
#' @description
#' Maps a pair of mean and c-statistic value to the parameters of a logit-normal distribution
#'
#' @param target A vector of size 2. The first element is mean and the second element is c-statistic.
#' @param method Either empty string, which invoked the default method; or "meansolve" which uses two 1-dimensional optimization approach.
#' @param integrate_controls (optional): parameters to be passed to integrate()
#' @param optim_controls (optional): parameters to be passed to optim()
#' @return A vector of size two that contains the distribution parameters
#' @examples
#' mcmap_logitnorm(c(0.1, 0.75))
#' @export
mcmap_logitnorm <- function(target=c(m=0.25,c=0.75), method="", integrate_controls=list(), optim_controls=list())
{
  if(method=="")
  {
    return(mcmap_logitnorm_default(target, integrate_controls, optim_controls))
  }

  if(method=="meansolve")
  {
    return(mcmap_logitnorm_meansolve_optim(target=target, integrate_controls, optim_controls))
  }

  stop("The requested method is not supplied.")
}






mcmap_logitnorm_default <- function(target=c(m=0.25,c=0.75), integrate_controls=list(), optim_controls=list())
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- mcmap_logitnorm_default(c(1-m,c), integrate_controls, optim_controls)
    if(is.null(tmp)) return(NULL);
    return(c(-tmp[1],tmp[2]))
  }

  F1 <- 1-m
  F2 <- 1-(2*c*m-(2*c-1)*m^2)

  if(is.null(integrate_controls$lower)) integrate_controls$lower<-0
  if(is.null(integrate_controls$upper)) integrate_controls$upper<-1
  integrate_controls$f <- plogitnorm

  integrate_controls2 <- integrate_controls
  integrate_controls2$f <- function(x, mu, sigma) {plogitnorm(x,mu,sigma)^2}

  f <- function(x)
  {
    f1 <- do.call(stats::integrate, args=c(integrate_controls, x[1], exp(x[2])))$value
    f2 <- do.call(stats::integrate, args=c(integrate_controls2,x[1], exp(x[2])))$value
    (f1-F1)^2+(f2-F2)^2
  }

  if(is.null(optim_controls$par)) optim_controls$par<-c(logit(m),0)
  if(is.null(optim_controls$method))
  {
    optim_controls$method <- "Nelder-Mead"
    optim_controls$control=list(maxit=1000)
  }

  optim_controls$fn <- f
  res <- do.call(stats::optim, args=optim_controls) #, control=list(trace=100))

  if(res$convergence==0)
    c(mu=unname(res$par[1]), sigma=unname(exp(res$par[2])))
  else
    NULL
}







#This one reduces the problem into two root-finding ones
mcmap_logitnorm_meansolve_optim <- function(target=c(m=0.25,c=0.75), integrate_controls=list(), optim_controls=list())
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- mcmap_logitnorm_meansolve_optim(c(1-m,c), integrate_controls, optim_controls)
    return(c(-tmp[1],tmp[2]))
  }

  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  if(is.null(integrate_controls$lower)) integrate_controls$lower<-0
  if(is.null(integrate_controls$upper)) integrate_controls$upper<-1
  integrate_controls$f <- function(x, sigma) {plogitnorm(x, find_logitnorm_mu(m, sigma), sigma)^2}

  f <- function(x)
  {
    f2 <- do.call(stats::integrate, args=c(integrate_controls, x))$value
    (f2-F2)^2
  }

  if(is.null(optim_controls$par)) optim_controls$par <- 0.5 #sigma
  if(is.null(optim_controls$method)) optim_controls$method <- "Brent"
  if(is.null(optim_controls$lower)) optim_controls$lower <- 0.00001
  if(is.null(optim_controls$upper)) optim_controls$upper <- 10
  optim_controls$fn <- f
  res <- do.call(stats::optim, args=optim_controls)
  sigma <- res$par[1]

  if(res$convergence==0)
    c(mu=find_logitnorm_mu(m,sigma), sigma=sigma)
  else
    NULL

}




#
# mcmap_logitnorm_meansolve_uniroot <- function(target=c(m=0.25,c=0.75), integrate_controls=list(), optim_controls=list())
# {
#   m <- target[1]
#   c <- target[2]
#
#   if(m>0.5)
#   {
#     tmp <- mcmap_logitnorm_meansolve_uniroot(c(1-m,c), integrate_controls, optim_controls)
#     return(c(-tmp[1],tmp[2]))
#   }
#
#   F1 <- 1-m
#   F2 <- 1- (2*c*m-(2*c-1)*m^2)
#
#   if(is.null(integrate_controls$lower)) integrate_controls$lower<-0
#   if(is.null(integrate_controls$upper)) integrate_controls$upper<-1
#   integrate_controls$f <- function(x, sigma) {plogitnorm(x, find_logitnorm_mu(m, sigma), sigma)^2}
#
#   f <- function(x)
#   {
#     f2 <- do.call(stats::integrate, args=c(integrate_controls, x))$value
#     (f2-F2)
#   }
#
#   if(is.null(optim_controls$interval))
#   {
#     interval <- c(0.001,1)
#     repeat
#     {
#       if(f(interval[1])>0) break;
#       interval[1] <- interval[1]/2
#     }
#     repeat
#     {
#       if(f(interval[2])<0) break;
#       interval[1] <- interval[2]
#       interval[2] <- interval[2]+1
#     }
#
#     optim_controls$interval <- interval
#   }
#
#   #optim_controls$interval <- c(0,1)
#   #if(is.null(optim_controls$par)) optim_controls$par <- 1 #sigma
#   optim_controls$f <- f
#   res <- do.call(stats::uniroot, args=optim_controls)
#
#   sigma <- res$root
#
#   c(mu=find_logitnorm_mu(m,sigma), sigma=sigma)
# }
#

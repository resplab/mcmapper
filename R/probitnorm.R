
#'@export
dprobitnorm <- function(x, mu, sigma) {dnorm(qnorm(x),mu,sigma)/dnorm(qnorm(x))}
#'@export
pprobitnorm <- function(x, mu, sigma) {pnorm(qnorm(x),mu,sigma)}
#'@export
rprobitnorm <- function(n,mu,sigma){pnorm(rnorm(n,mu,sigma))}




#'@export
mcmap_probitnorm <- function(target=c(m=0.25,c=0.75), method="", integrate_controls=list(), optim_controls=list())
{
  if(method=="")
  {
    return(mcmap_probitnorm_default(target, integrate_controls, optim_controls))
  }
  else
  {
    stop("The requested method is not supplied.")
  }
}



mcmap_probitnorm_default <- function(target=c(m=0.25,c=0.75), integrate_controls=list(), optim_controls=list())
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- mcmap_probitnorm(c(1-m,c), integrate_controls, optim_controls)
    return(c(-tmp[1],tmp[2]))
  }

  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  if(is.null(integrate_controls$lower)) integrate_controls$lower<-0
  if(is.null(integrate_controls$upper)) integrate_controls$upper<-1
  integrate_controls$f <- function(x, sigma) {pprobitnorm(x,qnorm(m)*sqrt(1+sigma^2),sigma)^2}

  f <- function(x)
  {
    f2 <- do.call(integrate, args=c(integrate_controls, x))$value
    (f2-F2)^2
  }

  if(is.null(optim_controls$par)) optim_controls$par <- 1 #sigma
  if(is.null(optim_controls$method)) optim_controls$method <- "Brent"
  if(is.null(optim_controls$lower)) optim_controls$lower <- 0.0001
  if(is.null(optim_controls$upper)) optim_controls$upper <- 10
  optim_controls$fn <- f
  res <- do.call(optim, args=optim_controls)

  if(res$convergence==0)
    c(mu=unname(qnorm(m)*sqrt(1+res$par[1]^2)), sigma=unname(res$par[1]))
  else
    NULL
}






#Removes sigma and solves for mu. Given weird ranges on mu it is obsolete in favor of mcmap_probitnorm
mcmap_probitnorm2 <- function(target=c(m=0.25,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]

  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    L <- 0
    U <- 1
    f2 <- integrate(function(x, mu) {pprobitnorm(x, mu, sqrt((mu/qnorm(m))^2-1))^2}, L, U, mu=x)$value

    (f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(sign(m-0.5))
  }
  if(m<0.5)
  {
    L <- -10
    U <- qnorm(m)
  }
  else
  {
    L <- qnorm(m)
    U <- 10
  }

  res <- optim((L+U)/2, f, method="Brent", lower=L, upper=U) #, control=list(trace=100))

  if(res$convergence==0)
    c(mu=res$par[1], sigma=sqrt((res$par[1]/qnorm(m))^2-1))
  else
    NULL
}



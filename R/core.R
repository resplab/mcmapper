
logit <- function(x) {log(x/(1-x))}
expit <- function(x) {1/(1+exp(-x))}

#'@export
dlogitnormal <- function(x, mu, sigma) {1/sigma/sqrt(2*pi)/x/(1-x)*exp(-(logit(x)-mu)^2/(2*sigma^2))}
#'@export
plogitnormal <- function(x, mu, sigma) {pnorm(logit(x),mu,sigma)}

#'@export
dprobitnormal <- function(x, mu, sigma) {dnorm(qnorm(x),mu,sigma)/dnorm(qnorm(x))}
#'@export
pprobitnormal <- function(x, mu, sigma) {pnorm(qnorm(x),mu,sigma)}



#'@export
solve_logitnormal <- function(target=c(m=0.5,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- solve_logitnormal(c(1-m,c), init)
    return(c(-tmp[1],tmp[2]))
  }

  F1 <- 1-m
  F2 <- 1-(2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    L <- 0#expit(x[1]-6*x[2])
    U <- 1#expit(x[1]+6*x[2])
    f1 <- integrate(plogitnormal, L, U , mu=x[1],sigma=(x[2]))$value
    f2 <- integrate(function(x, mu, sigma) {plogitnormal(x,mu,sigma)^2} , L, U ,mu=x[1], sigma=(x[2]))$value

    (f1-F1)^2+(f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(logit(m), (0.2))
  }
  res <- optim(init, f, method="L-BFGS-B", lower=c(-15,0), upper=c(10,10)) #, control=list(trace=100))

  if(res$convergence==0)
    c(mu=res$par[1], sigma=(res$par[2]))
  else
    NULL
}


#'@export
solve_generic <- function(target=c(m=0.5,c=0.75), CDF, init)
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
    f1 <- integrate(CDF, L, U , x[1], x[2])$value
    f2 <- integrate(function(x, par1, par2) {CDF(x,par1,par2)^2} , L, U ,par1=x[1], pa2=x[2])$value

    (f1-F1)^2+(f2-F2)^2
  }

  res <- optim(init, f) #, control=list(trace=100))

  if(res$convergence==0)
    res$par
  else
    NULL
}





#Removes sigma and solves for mu. Given weird ranges on mu it is obsolete in favor of V3.
solve_probitnormal2 <- function(target=c(m=0.5,c=0.75), init=NULL)
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
    f2 <- integrate(function(x, mu) {pprobitnormal(x, mu, sqrt((mu/qnorm(m))^2-1))^2}, L, U, mu=x)$value

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



#'@export
solve_probitnormal <- function(target=c(m=0.5,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- solve_probitnormal(c(1-m,c), init)
    return(c(-tmp[1],tmp[2]))
  }

  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    L <- 0
    U <- 1
    f2 <- integrate(function(x, sigma) {pprobitnormal(x,qnorm(m)*sqrt(1+sigma^2),sigma)^2}, L, U, sigma=x)$value

    (f2-F2)^2
  }

  L=0
  U=10

  if(is.null(init))
  {
    init <- 1
  }

  res <- optim(init, f, method="Brent", lower=L, upper=U) #, control=list(trace=100))

  if(res$convergence==0)
    c(mu=qnorm(m)*sqrt(1+res$par[1]^2), sigma=res$par[1])
  else
    NULL
}








#'@export
solve_beta <- function(target=c(m=0.5,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- solve_beta(c(1-m,c), init)
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
    f2 <- integrate(function(x, alpha) {pbeta(x, alpha, alpha*(1-m)/m)^2}, L, U, alpha=x)$value

    (f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(1)
  }
  res <- optim(init, f, method="Brent", lower=c(0.0001), upper=c(10000)) #, control=list(trace=100))

  if(res$convergence==0)
    c(alpha=res$par[1], beta=res$par[1]*(1-m)/m)
  else
    NULL
}





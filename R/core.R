
logit <- function(x) {log(x/(1-x))}
expit <- function(x) {1/(1+exp(-x))}

#'@export
rlogitnorm <- function(n, mu, sigma) {1/(1+exp(-(rnorm(n,mu,sigma))))}
#'@export
dlogitnorm <- function(x, mu, sigma) {1/sigma/sqrt(2*pi)/x/(1-x)*exp(-(logit(x)-mu)^2/(2*sigma^2))}
#'@export
plogitnorm <- function(x, mu, sigma) {pnorm(logit(x),mu,sigma)}
#'@export
qlogitnorm <- function(x, mu, sigma) {1/(1+exp(-(qnorm(x,mu,sigma))))}


#'@export
dprobitnorm <- function(x, mu, sigma) {dnorm(qnorm(x),mu,sigma)/dnorm(qnorm(x))}
#'@export
pprobitnorm <- function(x, mu, sigma) {pnorm(qnorm(x),mu,sigma)}


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


#'@export
elogitnorm <- function(mu,sigma)
{
  #message(mu,",",sigma)
  steps <- floor(-mu/sigma^2)
  if(steps>1000)
  {
    return(1-integrate(plogitnorm,0,1,mu=mu,sigma=sigma)$value)
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
    res <- uniroot(f,interval=c(0,1))$root
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
        res <- c(res, uniroot(f,interval=c(y[i-1],y[i]))$root)
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
    out[i] <- uniroot(function(x,s) {elogitnorm(x,s)-m}, interval, sigma[i])$root
  }
  out
}


#This is based on the general algorithm, only that the limits are fine-tunes for logitnormal
#Works with ln(sigma) to avoid constraints
#'@export
mcmap_logitnorm <- function(target=c(m=0.25,c=0.75), init=NULL, optim_method="Nelder-Mead")
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- mcmap_logitnorm(c(1-m,c), init)
    if(is.null(tmp)) return(NULL);
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
    f1 <- integrate(plogitnorm, L, U , mu=x[1],sigma=exp(x[2]))$value
    f2 <- integrate(function(x, mu, sigma) {plogitnorm(x,mu,sigma)^2} , L, U ,mu=x[1], sigma=exp(x[2]))$value

    (f1-F1)^2+(f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(logit(m),0)
  }

  if(optim_method=="L-BFGS-B")
  {
    res <- optim(init, f, method=optim_method, lower=c(-15,-10), upper=c(10,10))
  }
  else
  {
    res <- optim(init, f, method=optim_method, control=list(maxit=1000))
  }

  if(res$convergence==0)
    c(mu=res$par[1], sigma=exp(res$par[2]))
  else
    NULL
}



#More advanced based on dividing the integral if required
#'@export
mcmap_logitnorm2 <- function(target=c(m=0.25,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- mcmap_logitnorm(c(1-m,c), init)
    return(c(-tmp[1],tmp[2]))
  }

  F1 <- 1-m
  F2 <- 1-(2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(x)
  {
    #message(paste(x,collapse=","))
    modes <- find_logitnorm_modes(x[1],x[2])
    if(length(modes)==1)
    {
      if(dlogitnorm(modes,x[1],x[2])>10) #Sharp mode, break it around 0.01 and 0.99
      {
        points <- c(0, qlogitnorm(c(0.01,0.99),x[1],x[2]),1)
        message("Breaking!")
      }
      else
      {
        points <- c(0,1)
      }
    }
    else #Bimodal should always be broken! For now around the local minimum in the middle
    {
      plow <- pnorm(modes[2],x[1],x[2])
      points <- c(0,qlogitnorm(c(0.01,0.99)*plow,x[1],x[2]), qlogitnorm(plow+c(0.01,0.99)*(1-plow),x[1],x[2]),1)
      message("Breaking bimodal:",paste(points,collapse=","))
    }

    f1<- f2 <- 0
    for(i in 2:length(points))
    {
      f1 <- f1+integrate(plogitnorm, points[i-1], points[i] , mu=x[1],sigma=(x[2]))$value
      f2 <- f2+integrate(function(x, mu, sigma) {plogitnorm(x,mu,sigma)^2} , points[i-1], points[i] ,mu=x[1], sigma=(x[2]))$value
    }

    (f1-F1)^2+(f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(logit(m), (0.2))
  }
  res <- optim(init, f, method="L-BFGS-B", lower=c(-15,0), upper=c(10,10))

  if(res$convergence==0)
    c(mu=res$par[1], sigma=(res$par[2]))
  else
    NULL
}




#This one reduces the prblem into two root-finding ones
#'@export
mcmap_logitnorm3 <- function(target=c(m=0.25,c=0.75), interval=NULL)
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- mcmap_logitnorm3(c(1-m,c), interval)
    return(c(-tmp[1],tmp[2]))
  }

  F1 <- 1-m
  F2 <- 1- (2*c*m-(2*c-1)*m^2)

  #message(paste0("F1:",F1," , F2:",F2))

  f <- function(sigma)
  {
    f2 <- integrate(function(y) {plogitnorm(y, find_logitnorm_mu(m, sigma), sigma)^2}, 0, 1)$value

    f2-F2
  }

  if(is.null(interval))
  {
    interval <- c(0.1,1)
    repeat
    {
      if(f(interval[1])>0) break;
      interval[1] <- interval[1]/2
    }
    repeat
    {
      if(f(interval[2])<0) break;
      interval[1] <- interval[2]
      interval[2] <- interval[2]+1
    }
  }


  res <- uniroot(f, interval=interval)

  sigma <- res$root
  c(mu=find_logitnorm_mu(m,sigma), sigma=sigma)
}





#'@export
mcmap_generic <- function(target=c(m=0.25,c=0.75), CDF, init)
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
    f2 <- integrate(function(x, par1, par2) {CDF(x,par1,par2)^2} , L, U ,par1=x[1], par2=x[2])$value

    (f1-F1)^2+(f2-F2)^2
  }

  res <- optim(init, f) #, control=list(trace=100))

  if(res$convergence==0)
    res$par
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



#'@export
mcmap_probitnorm <- function(target=c(m=0.25,c=0.75), init=NULL)
{
  m <- target[1]
  c <- target[2]

  if(m>0.5)
  {
    tmp <- mcmap_probitnorm(c(1-m,c), init)
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
    f2 <- integrate(function(x, sigma) {pprobitnorm(x,qnorm(m)*sqrt(1+sigma^2),sigma)^2}, L, U, sigma=x)$value

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
mcmap_beta <- function(target=c(m=0.25,c=0.75), init=NULL)
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



#Log-transformed alpha and beta
#'@export
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
    f2 <- integrate(function(x, ln_alpha) {pbeta(x, exp(ln_alpha), exp(ln_alpha)*(1-m)/m)^2}, L, U, ln_alpha=x)$value

    (f2-F2)^2
  }

  if(is.null(init))
  {
    init <- c(0)
  }
  res <- optim(init, f, method="Brent", lower=-c(10), upper=c(10)) #, control=list(trace=100))

  if(res$convergence==0)
    c(alpha=exp(res$par[1]), beta=exp(res$par[1])*(1-m)/m)
  else
    NULL
}


#'@export
mcmap <- function(target=c(m=0.75, c=0.80), type=c("beta", "logitnorm", "probitnorm"), init=NULL)
{
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
    stop("Error: inadmissible value for m and/or c.")
  }

  if(abs(m-0.5)>0.499  | abs(c-0.75)>0.24)
  {
    warning("Extreme values for parameters are requested. The algorithms might not be reliable. Check the results (e.g., via simulation).")
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
        mcmap_logitnorm3(c(m,c))
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



#'@export
plot.mcmapper_output <- function(mcmapper_output, CDF=T, bins=1000, ...)
{
  x <- (0:bins)/bins

  strFun <- paste0(ifelse(CDF,"p","d"), mcmapper_output$type)

  tmp <- as.list(c(NA,res$value))
  tmp[[1]] <- x

  y <- do.call(strFun,args=tmp)

  plot(x,y,...)
}

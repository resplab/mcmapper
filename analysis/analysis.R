simulate <- function(args, type=c('beta','logitnormal','probitnormal'), n_sim=10^6)
{
  if(type=="logitnormal")
  {
    pi <- 1/(1+exp(-rnorm(n_sim,args[1],args[2])))
    Y <- rbinom(n_sim,1,pi)
  }

  if(type=="probitnormal")
  {
    pi <- pnorm(rnorm(n_sim,args[1],args[2]))
    Y <- rbinom(n_sim,1,pi)
  }

  if(type=="beta")
  {
    pi <- rbeta(n_sim,args[1],args[2])
    Y <- rbinom(n_sim,1,pi)
  }

  require(pROC)

  c(m=mean(pi), c=pROC::roc(Y~pi,)$auc)
}




m <- 0.086
c <- 0.81

res_logitnormal <- solve_logitnormal(c(m,c))
res_probitnormal <- solve_probitnormal(c(m,c))
res_beta <- solve_beta(c(m,c))

x <- (0:(10000-1))/10000
plot(x, dlogitnormal(x,res_logitnormal[1],res_logitnormal[2]), type='l')
lines(x, dprobitnormal(x,res_probitnormal[1],res_probitnormal[2]), type='l', col='blue')
lines(x, dbeta(x,res_beta[1],res_beta[2]), type='l', col='red')

plot(x, plogitnormal(x,res_logitnormal[1],res_logitnormal[2]), type='l')
lines(x, pprobitnormal(x,res_probitnormal[1],res_probitnormal[2]), type='l', col='blue')
lines(x, pbeta(x,res_beta[1],res_beta[2]), type='l', col='red')



populate_table <- function()
{
  out <- data.frame(m=double(), c=double(), type=character(), parm1=double(), parm2=double())

  ms <- (1:99)/100
  cs <- (55:95)/100

  N <- length(ms)*length(cs)*3

  out[N,1] <- 0

  index <- 1
  for(m in ms)
    for(c in cs)
    {
      cat(c(m,c),"|");

      res <- solve_logitnormal(c(m,c)); if(is.null(res)) {res<-c(NA,NA); message("Bad")}
      out[index,] <- list(m,c,"logitnormal",res[1],res[2])
      index <- index+1

      res <- solve_probitnormal(c(m,c)); if(is.null(res)) {res<-c(NA,NA); message("Bad")}
      out[index,] <- list(m,c,"probitnormal",res[1],res[2])
      index <- index+1

      res <- solve_beta(c(m,c)); if(is.null(res)) {res<-c(NA,NA); message("Bad")}
      out[index,] <- list(m,c,"beta",res[1],res[2])
      index <- index+1
    }

  out
}




plot_dist <- function(args, type=c('beta','logitnormal','probitnormal'), CDF=T, n_bins=10^3)
{
  x <- seq(from=0, to=1, length.out=n_bins)
  tmp <- as.list(c(NA,unname(args)))
  tmp[[1]] <- x
  if(CDF) pfx<-"p" else pfx<-"d";
  y <- do.call(args=tmp, what=paste0(pfx,type))

  plot(x,y,type='l')
}

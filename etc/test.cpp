#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
double Celogitnormal(double mu, double sigma, int N)
{
  Rprintf("%f,%f\n",mu,sigma);

  const double pi=3.141592653589793238462643383279502884;

  if(mu>0)
  {
    return(1-Celogitnormal(-mu,sigma,N));
  }
  if(abs(-mu-sigma*sigma)<abs(mu))
  {
    return(exp(mu+sigma*sigma/2)*Celogitnormal(-mu-sigma*sigma,sigma,N));
  }

  double A=0, B=0, C=0;
  for(int n=1;n<=N;n++)
  {
    A += exp(-sigma*sigma*n*n/2)*sinh(n*mu)*tanh(n*sigma*sigma/2);
    B += exp(-(2*n-1)*(2*n-1)*pi*pi/(2*sigma*sigma))*sin((2*n-1)*pi*mu/(sigma*sigma))/sinh((2*n-1)*pi*pi/(sigma*sigma));
    C += exp(-sigma*sigma*n*n/2)*cosh(n*mu);
    //Rprintf("%f,%f,%f; %f\n",A,B,C,-pi);
  }

  return (double)0.5+(A+2*pi/(sigma*sigma)*B)/(1+2*C);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

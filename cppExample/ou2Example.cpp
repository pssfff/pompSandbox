// -*- C++ -*-

//#include <R.h> // From pomp template
//#include <Rmath.h> // From pomp template
//#include <Rdefines.h> // From pomp template
#include <Rcpp.h> // From Rcpp template
using namespace Rcpp;



/*** R

procSim_R <- discrete.time.sim(
        step.fun=function (x, t, params, ...) {
          eps <- rnorm(n=2,mean=0,sd=1) # noise terms
          xnew <- c(
                    x1=params["alpha.1"]*x["x1"]+params["alpha.3"]*x["x2"]+
                        params["sigma.1"]*eps[1],
                    x2=params["alpha.2"]*x["x1"]+params["alpha.4"]*x["x2"]+
                        params["sigma.2"]*eps[1]+params["sigma.3"]*eps[2]
                    )
          names(xnew) <- c("x1","x2")
xnew }
        )
        
*/
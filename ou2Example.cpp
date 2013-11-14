#include <Rcpp.h> 
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector convolveCpp(NumericVector a, NumericVector b) {
    int na = a.size(), nb = b.size();
    int nab = na + nb - 1;
    NumericVector xab(nab);
for (int i = 0; i < na; i++)
for (int j = 0; j < nb; j++)
            xab[i + j] += a[i] * b[j];
return xab; }

/*** R
## this is the rprocess function written in R
procSim_R <- discrete.time.sim(
        step.fun=function (x, t, params, ...) {
                eps <- rnorm(n=2,mean=0,sd=1) # noise terms
                xnew <- c(x1=params["alpha.1"]*x["x1"]+params["alpha.3"]*x["x2"]+ params["sigma.1"]*eps[1],
                          x2=params["alpha.2"]*x["x1"]+params["alpha.4"]*x["x2"]+params["sigma.2"]*eps[1]+params["sigma.3"]*eps[2])
                names(xnew) <- c("x1","x2")
                xnew 
        }
)

## this command shows the examples in C code 
## file.show(file=system.file("examples/ou2.c",package="pomp"))

*/
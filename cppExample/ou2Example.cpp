// -*- C++ -*-

// according to http://www.lindonslog.com/programming/r/rcpp/, the following
// two includes are unnecessary when using Rcpp.h
// #include <Rdefines.h> // From pomp template
// #include <R.h> // From pomp template

#include <Rmath.h> // From pomp template
#include <Rcpp.h> // From Rcpp template
using namespace Rcpp;


// simple 2D Ornstein-Uhlenbeck process simulation
// [[Rcpp::export]]
static void sim_ou2 (double *x,
                     double alpha1, double alpha2, double alpha3, double alpha4, 
		     double sigma1, double sigma2, double sigma3)
{
  double eps[2], xnew[2];

  if (!(R_FINITE(x[0]))) return;
  if (!(R_FINITE(x[1]))) return;
  if (!(R_FINITE(alpha1))) return;
  if (!(R_FINITE(alpha2))) return;
  if (!(R_FINITE(alpha3))) return;
  if (!(R_FINITE(alpha4))) return;
  if (!(R_FINITE(sigma1))) return;
  if (!(R_FINITE(sigma2))) return;
  if (!(R_FINITE(sigma3))) return;

  eps[0] = rnorm(0,1);
  eps[1] = rnorm(0,1);

  xnew[0] = alpha1*x[0]+alpha3*x[1]+sigma1*eps[0];
  xnew[1] = alpha2*x[0]+alpha4*x[1]+sigma2*eps[0]+sigma3*eps[1];

  x[0] = xnew[0];
  x[1] = xnew[1];
}

#define ALPHA1     (p[parindex[0]])
#define ALPHA2     (p[parindex[1]])
#define ALPHA3     (p[parindex[2]])
#define ALPHA4     (p[parindex[3]])
#define SIGMA1     (p[parindex[4]])
#define SIGMA2     (p[parindex[5]])
#define SIGMA3     (p[parindex[6]])
#define TAU        (p[parindex[7]])

#define X1    (x[stateindex[0]])
#define X2    (x[stateindex[1]])
#define Y1    (y[obsindex[0]])
#define Y2    (y[obsindex[1]])

// onestep simulator for use in 'discrete.time.sim' plug-in
// [[Rcpp::export]]
void ou2_step (double *x, const double *p,
               const int *stateindex, const int *parindex, const int *covindex,
	       int ncovars, const double *covars,
	       double t, double dt) 
{
  sim_ou2(x,ALPHA1,ALPHA2,ALPHA3,ALPHA4,SIGMA1,SIGMA2,SIGMA3);
}

#undef ALPHA1
#undef ALPHA2
#undef ALPHA3
#undef ALPHA4
#undef SIGMA1
#undef SIGMA2
#undef SIGMA3
#undef TAU

#undef X1
#undef X2
#undef Y1
#undef Y2




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
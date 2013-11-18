#include <Rcpp.h>

using namespace Rcpp;

// testing the define syntax at the top to abstract a variable call
#define p (k)


// [[Rcpp::export]]
int timesTwo(int x) {
   return x * 2;
}

// [[Rcpp::export]]
int timesK(int x, int k) {
   return x * p;
}

## foo.R
setwd("~/Documents/code_versioned/pompSandbox/cppExample")

library(Rcpp)

sourceCpp("foo.cpp")
timesTwo(10)
timesK(10, 3)

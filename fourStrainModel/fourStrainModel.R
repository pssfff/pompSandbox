
setwd("~/Documents/code_versioned/pompSandbox/fourStrainModel/")
require(pomp)
source("fourStrainFuncs.R")


############################
## create base TSIR model ##
############################

## make spline basis
tbasis <- seq(0, 100, by=1/26)
basis <- periodic.bspline.basis(tbasis,nbasis=3,degree=2,period=1,names="b%d")
basisCoefs <- c(1.5,-.1,1) ## chosen to give range between ~.2 and ~1.2
betas <- basis %*% basisCoefs
#plot(betas[1:52], type="l")


## create the TSIR pomp object 
tsirModel <- pomp(
        data=data.frame(
                time=seq(0, 100,by=1/26),
                y1=NA, y2=NA, y3=NA, y4=NA
        ),
        times="time",
        tcovar=tbasis,
        covar=basis,
        t0=0,
        rprocess=discrete.time.sim(
                step.fun=fourStrainProcSim,
                delta.t=1/26
        ),
        rmeasure=fourStrainMeasSim,
        dmeasure=fourStrainMeasDens,
        # initial condition parameters 
        ic.pars=c(paste0("I", 1:4, ".0"),
                  paste0("S", 1:4, ".0"),
                  paste0("C", 1:4, ".0")
        ),
        # names of the compartments
        comp.names=c(paste0("I", 1:4),
                     paste0("S", 1:4),
                     paste0("C", 1:4)
        ),
        ## from estimation scale to natural scale
        parameter.transform=function(params,...){
                exp.idx <- c("N", "mu", "iota", "lambda",
                             "S1.0", "I1.0", "C1.0",
                             "S2.0", "I2.0", "C2.0",
                             "S3.0", "I3.0", "C3.0",
                             "S4.0", "I4.0", "C4.0")
                params[exp.idx] <- exp(params[exp.idx])
                expit.idx <- paste0("rho", 1:4)
                params[expit.idx] <- exp(params[expit.idx])/(1+exp(params[expit.idx]))
                return(params)
        },
        ## from natural scale to estimation scale
        parameter.inv.transform=function(params,...){
                log.idx <- c("N", "mu", "iota", "lambda",
                             "S1.0", "I1.0", "C1.0",
                             "S2.0", "I2.0", "C2.0",
                             "S3.0", "I3.0", "C3.0",
                             "S4.0", "I4.0", "C4.0")
                params[log.idx] <- log(params[log.idx])
                logit.idx <- paste0("rho", 1:4)
                params[logit.idx] <- log(params[logit.idx]/(1-params[logit.idx]))
                return(params)
        }
) 


##########################################
## refining basic TSIR model --> model2 ##
##########################################

## first, modify the seasonal basis
## chosen to give range between ~.2 and ~1.2
basisCoefs <- c(1.1,.95,1) 
betas <- basis %*% basisCoefs
#plot(betas[1:52], type="l")

paramsModel2 <- c(N=50000,
                  mu=1/50,
                  rho1=.5, rho2=.1, rho3=.1, rho4=.1,
                  b1=basisCoefs[1], 
                  b2=basisCoefs[2],  
                  b3=basisCoefs[3], 
                  alpha1=1,alpha2=1,
                  iota=1/700,
                  lambda=1/50,
                  S1.0=10000,I1.0=100,C1.0=100,
                  S2.0=10000,I2.0=100,C2.0=100,
                  S3.0=10000,I3.0=100,C3.0=100,
                  S4.0=10000,I4.0=100,C4.0=100
)

## simulate data from TSIR model and select just a subset of the data for fitting
tsirModel2 <- simulate(
        tsirModel,
        params=paramsModel2,
        seed=677573454L
) 

tsirModel2_Short <- window(tsirModel2, start=90, end=100)
#plot(tsirModel2, variables=c("I1","I2", "I3", "I4"))
#plot(tsirModel2, variables=c("S1","S2", "S3", "S4"))




## start with the truth, adjusted for windowing   
index0 <- max(which(tsirModel2@times<90))## variables for timezero
theta.truth <- paramsModel2     
time0names <- c("S1.0", "I1.0", "C1.0", "S2.0", "I2.0", "C2.0", "S3.0", "I3.0", "C3.0", "S4.0", "I4.0", "C4.0")
theta.truth[time0names] <- states(tsirModel2)[,index0]

pf.truth <- pfilter(tsirModel2_Short, params=theta.truth, 
                    Np=4000, max.fail=261, tol=1e-15,
                    pred.mean=TRUE, filter.mean=TRUE)

## now a small lie (1% from the truth)
theta.lie.small <- theta.truth
theta.lie.small[-1] <-  theta.lie.small[-1] + theta.truth[-1]*.01

pf.lie.small <- pfilter(tsirModel2_Short, params=theta.lie.small, Np=4000, max.fail=261, tol=1e-15,
                        pred.mean=TRUE, filter.mean=TRUE)

## now a big lie (10% from the truth)
theta.lie.big <- theta.truth
theta.lie.big[-1] <-  theta.lie.big[-1] + theta.truth[-1]*.1

pf.lie.big <- pfilter(tsirModel2_Short, params=theta.lie.big, Np=4000, max.fail=261, tol=1e-15,
                      pred.mean=TRUE, filter.mean=TRUE)

save.image("pfiltersWithInitialConditions.rda")


## Comparing data with one step ahead predictions..
plot.resids(pf.truth, standardize=TRUE)
plot.resids(pf.truth, standardize=FALSE)
plot.means(pf.truth)

plot.resids(pf.lie.small, standardize=TRUE)
plot.resids(pf.lie.small, standardize=FALSE)
plot.means(pf.lie.small)

#plot.resids(pf.lie.big, standardize=TRUE) ## gives error?
#plot.resids(pf.lie.big, standardize=FALSE)
#plot.means(pf.lie.big)


print(paste('The likelihood of truth is: ',logLik(pf.truth), sep=''))
print(paste('The likelihood of a small lie is: ',logLik(pf.lie.small), sep=''))
print(paste('The likelihood of a big lie is: ',logLik(pf.lie.big), sep=''))

## rewrite process simulator in C 
## use mif to create an estimate



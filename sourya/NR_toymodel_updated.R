

require(pomp)


############################
## create base TSIR model ##
############################

fourStrainProcSim <- function(x, t, params, delta.t, covars, ...) {
        ## unpack the params vector:
        b <- params[c("b1", "b2", "b3")]
        beta <- b %*% covars#params["beta"]
        alpha1 <- params["alpha1"]
        alpha2 <- params["alpha2"]
        lambda <- params["lambda"]
        mu <- params["mu"]
        N <- params["N"]
        iota <- params["iota"]
        ## copy x for easy access
        newx <- x
        namesI <- c("I1", "I2", "I3", "I4")
        namesS <- c("S1", "S2", "S3", "S4")
        namesC <- c("C1", "C2", "C3", "C4")
        ## transitions for each strain
        newx[namesI] <- rpois(4, beta*((x[namesI]/N+iota)^alpha1)*(x[namesS]^alpha2))
        ## sums for newly cross-protecteds for each strain 
        sumMatrix <- 1 - diag(4) ## for strain i, sum all but i 
        NCP <- newx[namesI]%*%sumMatrix
        ## Poisson approximation to exponential losses, those moving from CP back to S
        NS <- rpois(4, lambda*x[namesC]) 
        ## udpate counts
        newx[namesS] <- pmax(x[namesS] + N*mu - sum(newx[namesI]) + NS, 0)
        newx[namesC] <- x[namesC] - NS + NCP
        return(newx[c(namesS, namesI, namesC)])
}

fourStrainMeasSim <- function(x, t, params, covars, ...) {
        namesR <- paste0("rho", 1:4)
        namesI <- paste0("I", 1:4)
        rhos <- params[namesR]
        ys <- rbinom(4, size=x[namesI], prob=rhos)
        unname(ys)
}


## make spline basis
tbasis <- seq(0, 100, by=1/26)
basis <- periodic.bspline.basis(tbasis,nbasis=3,degree=2,period=1,names="b%d")
basisCoefs <- c(1.5,-.1,1) ## chosen to give range between ~.2 and ~1.2
betas <- basis %*% basisCoefs
#plot(betas[1:52], type="l")


## creating the TSIR pomp object 
pomp(
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
        # initial condition parameters 
        ic.pars=c(paste0("I", 1:4, ".0"),
                  paste0("S", 1:4, ".0"),
                  paste0("C", 1:4, ".0")
        ),
        # names of the compartments
        comp.names=c(paste0("I", 1:4),
                     paste0("S", 1:4),
                     paste0("C", 1:4)
        )
) -> tsir


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

## create 
tsirModel2 <- simulate(
        tsir,
        params=paramsModel2,
        seed=677573454L
) 
#plot(tsirModel2, variables=c("I1","I2", "I3", "I4"))
#plot(tsirModel2, variables=c("S1","S2", "S3", "S4"))
tsirModel2 <- window(tsirModel2, start=90, end=100)



## updating to remove NAs
fourStrainMeasDens <- function(y, x, t, params, covars, log, ...) {
        namesR <- paste0("rho", 1:4)
        namesI <- paste0("I", 1:4)
        rhos <- params[namesR]
        tol <- 10e-15
        f <- ifelse(is.na(dbinom(y, size=x[namesI], prob=rhos, log=log)), 
        		tol, 
        		dbinom(y, size=x[namesI], prob=rhos, log=log) + tol)
        return(unname(prod(f)))
}

tsirModel2 <- pomp(
             tsirModel2,
             dmeasure=fourStrainMeasDens
             )



tsirModel2 <- pomp(
                tsirModel2,
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


## functions to plot residuals from pfilters
plot.resids <- function(pf, standardize=FALSE) {
        par(mfrow=c(2,2))
        par(mar=c(4,4,1,1))
        if(standardize){
                resid1 <- (pf@data['y1',] - pred.mean(pf)['I1',]*pf@params['rho1'])/sqrt(pred.mean(pf)['I1',]*pf@params['rho1'])
                resid2 <- (pf@data['y2',] - pred.mean(pf)['I2',]*pf@params['rho2'])/sqrt(pred.mean(pf)['I2',]*pf@params['rho2'])
                resid3 <- (pf@data['y3',] - pred.mean(pf)['I3',]*pf@params['rho3'])/sqrt(pred.mean(pf)['I3',]*pf@params['rho3'])
                resid4 <- (pf@data['y4',] - pred.mean(pf)['I4',]*pf@params['rho4'])/sqrt(pred.mean(pf)['I4',]*pf@params['rho4'])
        } else {
                resid1 <- pf@data['y1',] - pred.mean(pf)['I1',]*pf@params['rho1']
                resid2 <- pf@data['y2',] - pred.mean(pf)['I2',]*pf@params['rho2']
                resid3 <- pf@data['y3',] - pred.mean(pf)['I3',]*pf@params['rho3']
                resid4 <- pf@data['y4',] - pred.mean(pf)['I4',]*pf@params['rho4']                
        }
        ylim <- range(resid1, resid2, resid3, resid4)
        plot(resid1, type='o', col=2, ylab='y1', xlab='weeks', ylim=ylim)
        plot(resid2, type='o', col=2, ylab='y2', xlab='weeks', ylim=ylim)
        plot(resid3, type='o', col=2, ylab='y3', xlab='weeks', ylim=ylim)
        plot(resid4,type='o', col=2, ylab='y4', xlab='weeks', ylim=ylim)
        
}

plot.means <- function(pf) {
        par(mfrow=c(2,2))
        par(mar=c(4,4,1,1))
        mean1 <- pred.mean(pf)['I1',]*pf@params['rho1']
        mean2 <- pred.mean(pf)['I2',]*pf@params['rho2']
        mean3 <- pred.mean(pf)['I3',]*pf@params['rho3']
        mean4 <- pred.mean(pf)['I4',]*pf@params['rho4']                
        
        ylim <- range(mean1, mean2, mean3, mean4, pf@data[c("y1", "y2", "y3", "y4"),])
        plot(mean1, type='o', col=2, ylab='y1', xlab='weeks', ylim=ylim)
        lines(pf@data['y1',], type='l')
        plot(mean2, type='o', col=2, ylab='y2', xlab='weeks', ylim=ylim)
        lines(pf@data['y2',], type='l')
        plot(mean3, type='o', col=2, ylab='y3', xlab='weeks', ylim=ylim)
        lines(pf@data['y3',], type='l')
        plot(mean4,type='o', col=2, ylab='y4', xlab='weeks', ylim=ylim)
        lines(pf@data['y4',], type='l')
}

## start with the truth       
theta.truth <- paramsModel2     
            
pf.truth <- pfilter(tsirModel2, params=theta.truth, 
                    Np=4000, max.fail=261, tol=1e-15,
                    pred.mean=TRUE, filter.mean=TRUE)

## now a small lie (1% from the truth)
theta.lie.small <- theta.truth
theta.lie.small[-1] <-  theta.lie.small[-1] + theta.truth[-1]*.01

pf.lie.small <- pfilter(tsirModel2, params=theta.lie.small, Np=4000, max.fail=261, tol=1e-15,
                        pred.mean=TRUE, filter.mean=TRUE)

## now a big lie (10% from the truth)
theta.lie.big <- theta.truth
theta.lie.big[-1] <-  theta.lie.big[-1] + theta.truth[-1]*.1

pf.lie.big <- pfilter(tsirModel2, params=theta.lie.big, Np=4000, max.fail=261, tol=1e-15,
                      pred.mean=TRUE, filter.mean=TRUE)

## Comparing data with one step ahead predictions..
plot.resids(pf.truth, standardize=TRUE)
plot.resids(pf.truth, standardize=FALSE)
plot.means(pf.truth)

plot.resids(pf.lie.small, standardize=TRUE)
plot.resids(pf.lie.small, standardize=FALSE)
plot.means(pf.lie.small)

plot.resids(pf.lie.big, standardize=TRUE) ## gives error?
plot.resids(pf.lie.big, standardize=FALSE)
plot.means(pf.lie.big)


print(paste('The likelihood of truth is: ',logLik(pf.truth), sep=''))
print(paste('The likelihood of a small lie is: ',logLik(pf.lie.small), sep=''))
print(paste('The likelihood of a big lie is: ',logLik(pf.lie.big), sep=''))

## reset initial conditions
## rewrite process simulator in C 
## use mif to create an estimate



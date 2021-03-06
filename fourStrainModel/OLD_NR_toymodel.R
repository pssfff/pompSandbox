

require(pomp)
proc.sim <- function(x, t, params, delta.t, ...) {
        ## unpack the params vector:
        ##  biweek
        biweek <- params["biweek"]
        beta.string <- paste0("beta", biweek)
        ## beta1, ...., beta26
        beta_t <- params[beta.string]
        ## alpha1 and alpha2
        alpha1 <- params["alpha1"]
        alpha2 <- params["alpha2"]
        ## CP params: delta and k
        delta <- params["delta"]
        k <- params["k"]
        ## state at time t:
        X <- x["X"]
        ## compute the state at time t+delta.t
        log.lambda <- beta*log(X) + alpha
        xnew <- c(X=unname( rpois(n=1, lambda=exp(log.lambda)) ))
        return(xnew)
}



toy.proc.sim <- function(x, t, params, delta.t, ...) {
        ## unpack the params vector:
        beta <- params["beta"]
        alpha1 <- params["alpha1"]
        alpha2 <- params["alpha2"]
        lambda <- params["lambda"]
        mu <- params["mu"]
        N <- params["N"]
        iota <- params["iota"]
        ## copy x for easy access
        newx <- x
        ## transitions for each strain
        newx["I1"] <- rpois(1, beta*((x["I1"]/N+iota)^alpha1)*(x["S1"]^alpha2))
        newx["I2"] <- rpois(1, beta*((x["I2"]/N+iota)^alpha1)*(x["S2"]^alpha2))
        ## Poisson approximation to exponential losses
        C1.loss <- rpois(1, lambda*x["C1"]) 
        C2.loss <- rpois(1, lambda*x["C2"]) 
        ## udpate counts
        newx["S1"] <- x["S1"] + N*mu - newx["I1"] - newx["I2"] + C1.loss
        newx["S2"] <- x["S2"] + N*mu - newx["I2"] - newx["I1"] + C2.loss
        newx["C1"] <- x["C1"] - C1.loss + newx["I2"]
        newx["C2"] <- x["C2"] - C2.loss + newx["I1"]
        return(newx[c("S1", "I1", "C1", "S2", "I2", "C2")])
}

toy.meas.sim <- function(x, t, params, ...) {
        rho1 <- params["rho1"]
        rho2 <- params["rho2"]
        y1 <- rbinom(1, size=x["I1"], prob=rho1)
        y2 <- rbinom(1, size=x["I2"], prob=rho2)
        unname(c(y1, y2))
}



simulate(
        pomp(
                data=data.frame(
                        time=seq(0, 5000,by=1),
                        y1=NA,
                        y2=NA
                ),
                times="time",
                t0=0,
                rprocess=discrete.time.sim(
                        step.fun=toy.proc.sim,
                        delta.t=1
                ),
                rmeasure=toy.meas.sim,
                # initial condition parameters 
                ic.pars=c("S1.0","I1.0","C1.0", 
                          "S2.0","I2.0","C2.0"), 
                # names of the compartments
                comp.names=c("S1","I1","C1", "S2","I2","C2") 
        ),
        params=c(
                N=50000,
                mu=1/500,
                rho1=.5, rho2=.1,
                beta=1.4,alpha1=1,alpha2=1,
                iota=1/10000,
                lambda=1/100,
                S1.0=40000,I1.0=100,C1.0=100,
                S2.0=40000,I2.0=100,C2.0=100
        ),
        seed=677573454L
) -> tsir
#plot(tsir, variables=c("y1", "y2", "C1","C2", "I1","I2", "S1", "S2"))



tbasis <- seq(0, 100, by=1/26)
basis <- periodic.bspline.basis(tbasis,nbasis=3,degree=2,period=1,names="b%d")
basisCoefs <- c(1.5,-.1,1) ## chosen to give range between ~.2 and ~1.2
betas <- basis %*% basisCoefs
#plot(betas[1:52], type="l")



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



pomp(
        data=data.frame(
                time=seq(0, 10,by=1/26),
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



paramsZeroes <- c(N=50000,
                  mu=1/500,
                  rho1=.5, rho2=.1, rho3=.1, rho4=.1,
                  b1=basisCoefs[1], 
                  b2=basisCoefs[2],  
                  b3=basisCoefs[3], 
                  #beta=1.1,
                  alpha1=1,alpha2=1,
                  iota=1/200000,
                  lambda=1/50,
                  S1.0=10000,I1.0=100,C1.0=100,
                  S2.0=10000,I2.0=100,C2.0=100,
                  S3.0=10000,I3.0=100,C3.0=100,
                  S4.0=10000,I4.0=100,C4.0=100
)
simulate(
        tsir,
        params=paramsZeroes,
        seed=677573454L
) -> tsirZeroes
#plot(tsirZeroes, variables=c("I1","I2", "I3", "I4"))
#plot(tsirZeroes, variables=c("S1","S2", "S3", "S4"))



basisCoefs <- c(1.1,.95,1) ## chosen to give range between ~.2 and ~1.2
betas <- basis %*% basisCoefs
#plot(betas[1:52], type="l")
paramsMore <- c(N=50000,
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
simulate(
        tsir,
        params=paramsMore,
        seed=677573454L
) -> tsirMore
#plot(tsirMore, variables=c("I1","I2", "I3", "I4"))
#plot(tsirMore, variables=c("S1","S2", "S3", "S4"))
#tsirMore <- window(tsirMore, start=30, end=100)



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
tsirMore <- pomp(
             tsirMore,
             dmeasure=fourStrainMeasDens
             )



tsirMore <- pomp(
                tsirMore,
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


## start with the truth       
theta.truth <- paramsMore     
            
pf.truth <- pfilter(tsirMore, params=theta.truth, Np=2000, max.fail=261, tol=1e-15,
               pred.mean=TRUE, filter.mean=TRUE)

## Comparing data with one step ahead predictions..
pdf(file='predmean_truth.pdf')
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
plot(pf.truth@data['y1',] - pred.mean(pf.truth)['I1',]*pf.truth@params['rho1'], type='o', col=2, ylab='y1', xlab='weeks')
plot(pf.truth@data['y2',] - pred.mean(pf.truth)['I2',]*pf.truth@params['rho2'], type='o', col=2, ylab='y2', xlab='weeks')
plot(pf.truth@data['y3',] - pred.mean(pf.truth)['I3',]*pf.truth@params['rho3'], type='o', col=2, ylab='y3', xlab='weeks')
plot(pf.truth@data['y4',] - pred.mean(pf.truth)['I4',]*pf.truth@params['rho4'], type='o', col=2, ylab='y4', xlab='weeks')
dev.off()

print(paste('The likelihood of truth is: ',logLik(pf.truth), sep=''))

## now a small lie (1% from the truth)
theta.lie.small <- theta.truth
theta.lie.small[-1] <-  theta.lie.small[-1] + theta.truth[-1]*.01

pf.lie.small <- pfilter(tsirMore, params=theta.lie.small, Np=2000, max.fail=261, tol=1e-15,
               pred.mean=TRUE, filter.mean=TRUE)
print(paste('The likelihood of a small lie is: ',logLik(pf.lie.small), sep=''))
pdf(file='predmean_small_lie.pdf')
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
plot(pf.lie.small@data['y1',] - pred.mean(pf.lie.small)['I1',]*pf.lie.small@params['rho1'], type='o', col=2, ylab='y1', xlab='weeks')
plot(pf.lie.small@data['y2',] - pred.mean(pf.lie.small)['I2',]*pf.lie.small@params['rho2'], type='o', col=2, ylab='y2', xlab='weeks')
plot(pf.lie.small@data['y3',] - pred.mean(pf.lie.small)['I3',]*pf.lie.small@params['rho3'], type='o', col=2, ylab='y3', xlab='weeks')
plot(pf.lie.small@data['y4',] - pred.mean(pf.lie.small)['I4',]*pf.lie.small@params['rho4'], type='o', col=2, ylab='y4', xlab='weeks')
dev.off()


## now a big lie (10% from the truth)
theta.lie.big <- theta.truth
theta.lie.big[-1] <-  theta.lie.big[-1] + theta.truth[-1]*.1

pf.lie.big <- pfilter(tsirMore, params=theta.lie.big, Np=2000, max.fail=261, tol=1e-15,
               pred.mean=TRUE, filter.mean=TRUE)
print(paste('The likelihood of a big lie is: ',logLik(pf.lie.big), sep=''))

pdf(file='predmean_big_lie.pdf')
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
plot(pf.lie.big@data['y1',] - pred.mean(pf.lie.big)['I1',]*pf.lie.big@params['rho1'], type='o', col=2, ylab='y1', xlab='weeks')
plot(pf.lie.big@data['y2',] - pred.mean(pf.lie.big)['I2',]*pf.lie.big@params['rho2'], type='o', col=2, ylab='y2', xlab='weeks')
plot(pf.lie.big@data['y3',] - pred.mean(pf.lie.big)['I3',]*pf.lie.big@params['rho3'], type='o', col=2, ylab='y3', xlab='weeks')
plot(pf.lie.big@data['y4',] - pred.mean(pf.lie.big)['I4',]*pf.lie.big@params['rho4'], type='o', col=2, ylab='y4', xlab='weeks')
dev.off()



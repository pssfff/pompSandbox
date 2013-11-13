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


fourStrainMeasDens <- function(y, x, t, params, covars, log, ...) {
        namesR <- paste0("rho", 1:4)
        namesI <- paste0("I", 1:4)
        rhos <- params[namesR]
        tol <- 10e-15
        ## removing NAs
        f <- ifelse(is.na(dbinom(y, size=x[namesI], prob=rhos, log=log)), 
                    tol, 
                    dbinom(y, size=x[namesI], prob=rhos, log=log) + tol)
        return(unname(prod(f)))
}


## convenience plotting functions to plot residuals and means from pfilters
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
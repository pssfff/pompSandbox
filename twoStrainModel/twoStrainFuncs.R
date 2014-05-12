############################
## proc sim functions     ##
############################

twoStrainProcSim <- function(x, t, params, delta.t, covars, ...) {
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
        namesI <- c("I1", "I2")
        namesS <- c("S1", "S2")
        namesC <- c("C1", "C2")
        ## transitions for each strain
        newx[namesI] <- rpois(2, beta*((x[namesI]/N+iota)^alpha1)*(x[namesS]^alpha2))
        ## sums for newly cross-protecteds for each strain 
        sumMatrix <- 1 - diag(2) ## for strain i, sum all but i 
        NCP <- newx[namesI]%*%sumMatrix
        ## Poisson approximation to exponential losses, those moving from CP back to S
        NS <- rpois(2, lambda*x[namesC]) 
        ## udpate counts
        newx[namesS] <- pmax(x[namesS] + N*mu - sum(newx[namesI]) + NS, 0)
        newx[namesC] <- x[namesC] - NS + NCP
        return(newx[c(namesS, namesI, namesC)])
}

## skeleton functions in C
step.fn <- '
        // define beta, the transmission rate
        double beta;        		
        beta = b1*seas1+b2*seas2+b3*seas3;

        // new infections
        I1 = rpois(beta*pow(I1/N+iota, alpha1)*pow(S1, alpha2));
        I2 = rpois(beta*pow(I2/N+iota, alpha1)*pow(S2, alpha2));

        // Poisson approximation to exponential losses, transitions from CP to S
        int C1_loss = rpois(lambda*C1);
        int C2_loss = rpois(lambda*C2);

        // counts of newly cross-protected for each strain
        int C1_add = I2 ;
        int C2_add = I1 ;

        // udpate susceptible counts
        int allI = I1 + I2;
        S1 += N*mu - allI + C1_loss;
        if ( S1 <=0 ) S1 = 0;
        if ( S1 >=N ) S1 = N;
        S2 += N*mu - allI + C2_loss;
        if ( S2 <=0 ) S2 = 0;
        if ( S2 >=N ) S2 = N;

        // update cross-protected counts 
        C1 += C1_add - C1_loss;
        if ( C1 >=N ) C1 = N;
        C2 += C2_add - C2_loss;
        if ( C2 >=N ) C2 = N;
'
skel <- '
return;
'


#######################
## measure functions ##
#######################

twoStrainMeasSim <- function(x, t, params, covars, ...) {
        namesR <- paste0("rho", 1:2)
        namesI <- paste0("I", 1:2)
        rhos <- params[namesR]
        ys <- rbinom(2, size=x[namesI], prob=rhos)
        unname(ys)
}


twoStrainMeasDens <- function(y, x, t, params, covars, log, ...) {
        namesR <- paste0("rho", 1:2)
        namesI <- paste0("I", 1:2)
        rhos <- params[namesR]
        tol <- 10e-15
        ## removing NAs
        f <- ifelse(is.na(dbinom(y, size=x[namesI], prob=rhos, log=log)), 
                    tol, 
                    dbinom(y, size=x[namesI], prob=rhos, log=log) + tol)
        return(unname(prod(f)))
}

## measure functions for C
rmeas <- "
        cases1 = rbinom(I1, rho1);
        cases2 = rbinom(I2, rho2);
"

dmeas <- "
        double lik1 = dbinom(cases1, I1, rho1, give_log);
        double lik2 = dbinom(cases2, I2, rho2, give_log);
        double tol = 10E-15;
        lik = lik1*lik2;
        // need to add tolerance here?
"


###############################
## parameter transformations ##
###############################

## in R
partransR <- function(params,...){
        exp.idx <- c("N", "mu", "iota", "lambda",
                     "S1.0", "I1.0", "C1.0",
                     "S2.0", "I2.0", "C2.0")
        params[exp.idx] <- exp(params[exp.idx])
        expit.idx <- paste0("rho", 1:2)
        params[expit.idx] <- exp(params[expit.idx])/(1+exp(params[expit.idx]))
        return(params)
}

paruntransR  <- function(params,...){
        log.idx <- c("N", "mu", "iota", "lambda",
                     "S1.0", "I1.0", "C1.0",
                     "S2.0", "I2.0", "C2.0")
        params[log.idx] <- log(params[log.idx])
        logit.idx <- paste0("rho", 1:2)
        params[logit.idx] <- log(params[logit.idx]/(1-params[logit.idx]))
        return(params)
}

## in C
partrans <- "
        TN = exp(N);
        Tmu = exp(mu);
        Talpha1 = exp(alpha1);
        Talpha2 = exp(alpha2);
        Tiota = exp(iota);
        Tlambda = exp(lambda);
        // rhos        
        Trho1 = expit(rho1);
        Trho2 = expit(rho2);
        // strain 1
        TS1_0 = exp(S1_0);
        TI1_0 = exp(I1_0);
        TC1_0 = exp(C1_0);
        // strain 2
        TS2_0 = exp(S2_0);
        TI2_0 = exp(I2_0);
        TC2_0 = exp(C2_0);
"
paruntrans <- "
        TN = log(N);
        Tmu = log(mu);
        Talpha1 = log(alpha1);
        Talpha2 = log(alpha2);
        Tiota = log(iota);
        Tlambda = log(lambda);
        // rhos        
        Trho1 = logit(rho1);
        Trho2 = logit(rho2);
        // strain 1
        TS1_0 = log(S1_0);
        TI1_0 = log(I1_0);
        TC1_0 = log(C1_0);
        // strain 2
        TS2_0 = log(S2_0);
        TI2_0 = log(I2_0);
        TC2_0 = log(C2_0);
"

########################
## Create memory term ##
########################

createMemory <- function(pompDat, k, lambda=NULL) {
        ## assume that dat is a two-strain pomp object
        ## k is a duration
        ## lambda is an exponential parameter
        ## first, we will create a matrix to help with the memory term summations.
        ## this matrix will have 0s and 1s and will be multiplied by the data vectors
        
        dat <- pompDat@data
        ## vector of 1s of length k
        if(is.null(lambda)) {
                vec1 <- rep(1, k)
        } else {
                ## rough cumulative probabilities of still being susceptible
                vec1 <- 1-pexp((k-1):0, rate=1/lambda)
        }
        ## vector of 0s of length ncol(dat)-k+1
        vec0 <- rep(0, ncol(dat)-k+1)
        ## the matrix
        memMatrixRows <- ncol(dat)-k
        firstReps <- rep(c(vec1, vec0), times=memMatrixRows-1)
        last1s <- c(vec1,0)
        memMatrix <- matrix(c(firstReps, last1s), byrow=TRUE, nrow=memMatrixRows)
        M1 <- memMatrix%*%dat[2,]
        M2 <- memMatrix%*%dat[1,]
        
        dataIdx <- (k+1):ncol(dat)
        dat1 <- data.frame(y=dat[1,dataIdx], 
                           yAR = dat[1, dataIdx-1],
                           M=M1, 
                           time=pompDat@times[dataIdx],
                           strain=1)
        dat2 <- data.frame(y=dat[2,dataIdx], 
                           yAR = dat[1, dataIdx-1],
                           M=M2,
                           time=pompDat@times[dataIdx],
                           strain=2)
        data <- rbind(dat1, dat2)
        return(data)
}


######################################
## miscellaneous plotting functions ##
######################################

## convenience plotting functions to plot residuals and means from pfilters
plot.resids <- function(pf, standardize=FALSE) {
        par(mfrow=c(2,1))
        par(mar=c(4,4,1,1))
        if(standardize){
                resid1 <- (pf@data['cases1',] - pred.mean(pf)['I1',]*pf@params['rho1'])/sqrt(pred.mean(pf)['I1',]*pf@params['rho1'])
                resid2 <- (pf@data['cases2',] - pred.mean(pf)['I2',]*pf@params['rho2'])/sqrt(pred.mean(pf)['I2',]*pf@params['rho2'])
        } else {
                resid1 <- pf@data['cases1',] - pred.mean(pf)['I1',]*pf@params['rho1']
                resid2 <- pf@data['cases2',] - pred.mean(pf)['I2',]*pf@params['rho2']
        }
        ylim <- range(resid1, resid2)
        plot(resid1, type='o', col=2, ylab='cases1', xlab='weeks', ylim=ylim)
        plot(resid2, type='o', col=2, ylab='cases2', xlab='weeks', ylim=ylim)
        
}

plot.means <- function(pf, ...) {
        par(mfrow=c(2,1))
        par(mar=c(4,4,1,1))
        mean1 <- pred.mean(pf)['I1',...]*pf@params['rho1']
        mean2 <- pred.mean(pf)['I2',...]*pf@params['rho2']
        
        ylim <- range(mean1, mean2, pf@data[c("cases1", "cases2"),])
        plot(mean1, type='o', col=2, ylab='y1', xlab='weeks', ylim=ylim)
        lines(pf@data['cases1',], type='l')
        plot(mean2, type='o', col=2, ylab='y2', xlab='weeks', ylim=ylim)
        lines(pf@data['cases2',], type='l')
}
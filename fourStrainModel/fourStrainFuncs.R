############################
## proc sim functions     ##
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

## skeleton functions in C
step.fn <- '
        // define beta, the transmission rate
        double beta;        		
        beta = b1*seas1+b2*seas2+b3*seas3;

        // new infections
        I1 = rpois(beta*pow(I1/N+iota, alpha1)*pow(S1, alpha2));
        I2 = rpois(beta*pow(I2/N+iota, alpha1)*pow(S2, alpha2));
        I3 = rpois(beta*pow(I3/N+iota, alpha1)*pow(S3, alpha2));
        I4 = rpois(beta*pow(I4/N+iota, alpha1)*pow(S4, alpha2));

        // Poisson approximation to exponential losses, transitions from CP to S
        int C1_loss = rpois(lambda*C1);
        int C2_loss = rpois(lambda*C2);
        int C3_loss = rpois(lambda*C3);
        int C4_loss = rpois(lambda*C4);

        // counts of newly cross-protected for each strain
        int C1_add = I2 + I3 + I4;
        int C2_add = I1 + I3 + I4;
        int C3_add = I1 + I2 + I4;
        int C4_add = I1 + I2 + I3;

        // udpate susceptible counts
        int allI = I1 + I2 + I3 + I4;
        S1 += N*mu - allI + C1_loss;
        if ( S1 <=0 ) S1 = 0;
        if ( S1 >=N ) S1 = N;
        S2 += N*mu - allI + C2_loss;
        if ( S2 <=0 ) S2 = 0;
        if ( S2 >=N ) S2 = N;
        S3 += N*mu - allI + C3_loss;
        if ( S3 <=0 ) S3 = 0;
        if ( S3 >=N ) S3 = N;
        S4 += N*mu - allI + C4_loss;
        if ( S4 <=0 ) S4 = 0;
        if ( S4 >=N ) S4 = N;

        // update cross-protected counts 
        C1 += C1_add - C1_loss;
        if ( C1 >=N ) C1 = N;
        C2 += C2_add - C2_loss;
        if ( C2 >=N ) C2 = N;
        C3 += C3_add - C3_loss;
        if ( C3 >=N ) C3 = N;
        C4 += C4_add - C4_loss;
        if ( C4 >=N ) C4 = N;
'
skel <- '
return;
'


#######################
## measure functions ##
#######################

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

## measure functions for C
rmeas <- "
        cases1 = rbinom(I1, rho1);
        cases2 = rbinom(I2, rho2);
        cases3 = rbinom(I3, rho3);
        cases4 = rbinom(I4, rho4);
"

dmeas <- "
        double lik1 = dbinom(cases1, I1, rho1, give_log);
        double lik2 = dbinom(cases2, I2, rho2, give_log);
        double lik3 = dbinom(cases3, I3, rho3, give_log);
        double lik4 = dbinom(cases4, I4, rho4, give_log);
        double tol = 10E-15;
        lik = lik1*lik2*lik3*lik4;
        // need to add tolerance here?
"


###############################
## parameter transformations ##
###############################

## in R
partransR <- function(params,...){
        exp.idx <- c("N", "mu", "iota", "lambda",
                     "S1.0", "I1.0", "C1.0",
                     "S2.0", "I2.0", "C2.0",
                     "S3.0", "I3.0", "C3.0",
                     "S4.0", "I4.0", "C4.0")
        params[exp.idx] <- exp(params[exp.idx])
        expit.idx <- paste0("rho", 1:4)
        params[expit.idx] <- exp(params[expit.idx])/(1+exp(params[expit.idx]))
        return(params)
}

paruntransR  <- function(params,...){
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
        Trho3 = expit(rho3);
        Trho4 = expit(rho4);
        // strain 1
        TS1_0 = exp(S1_0);
        TI1_0 = exp(I1_0);
        TC1_0 = exp(C1_0);
        // strain 2
        TS2_0 = exp(S2_0);
        TI2_0 = exp(I2_0);
        TC2_0 = exp(C2_0);
        // strain 3
        TS3_0 = exp(S3_0);
        TI3_0 = exp(I3_0);
        TC3_0 = exp(C3_0);
        // strain 4
        TS4_0 = exp(S4_0);
        TI4_0 = exp(I4_0);
        TC4_0 = exp(C4_0);
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
        Trho3 = logit(rho3);
        Trho4 = logit(rho4);
        // strain 1
        TS1_0 = log(S1_0);
        TI1_0 = log(I1_0);
        TC1_0 = log(C1_0);
        // strain 2
        TS2_0 = log(S2_0);
        TI2_0 = log(I2_0);
        TC2_0 = log(C2_0);
        // strain 3
        TS3_0 = log(S3_0);
        TI3_0 = log(I3_0);
        TC3_0 = log(C3_0);
        // strain 4
        TS4_0 = log(S4_0);
        TI4_0 = log(I4_0);
        TC4_0 = log(C4_0);
"

######################################
## miscellaneous plotting functions ##
######################################

## convenience plotting functions to plot residuals and means from pfilters
plot.resids <- function(pf, standardize=FALSE) {
        par(mfrow=c(2,2))
        par(mar=c(4,4,1,1))
        if(standardize){
                resid1 <- (pf@data['cases1',] - pred.mean(pf)['I1',]*pf@params['rho1'])/sqrt(pred.mean(pf)['I1',]*pf@params['rho1'])
                resid2 <- (pf@data['cases2',] - pred.mean(pf)['I2',]*pf@params['rho2'])/sqrt(pred.mean(pf)['I2',]*pf@params['rho2'])
                resid3 <- (pf@data['cases3',] - pred.mean(pf)['I3',]*pf@params['rho3'])/sqrt(pred.mean(pf)['I3',]*pf@params['rho3'])
                resid4 <- (pf@data['cases4',] - pred.mean(pf)['I4',]*pf@params['rho4'])/sqrt(pred.mean(pf)['I4',]*pf@params['rho4'])
        } else {
                resid1 <- pf@data['cases1',] - pred.mean(pf)['I1',]*pf@params['rho1']
                resid2 <- pf@data['cases2',] - pred.mean(pf)['I2',]*pf@params['rho2']
                resid3 <- pf@data['cases3',] - pred.mean(pf)['I3',]*pf@params['rho3']
                resid4 <- pf@data['cases4',] - pred.mean(pf)['I4',]*pf@params['rho4']                
        }
        ylim <- range(resid1, resid2, resid3, resid4)
        plot(resid1, type='o', col=2, ylab='cases1', xlab='weeks', ylim=ylim)
        plot(resid2, type='o', col=2, ylab='cases2', xlab='weeks', ylim=ylim)
        plot(resid3, type='o', col=2, ylab='cases3', xlab='weeks', ylim=ylim)
        plot(resid4,type='o', col=2, ylab='cases4', xlab='weeks', ylim=ylim)
        
}

plot.means <- function(pf, ...) {
        par(mfrow=c(2,2))
        par(mar=c(4,4,1,1))
        mean1 <- pred.mean(pf)['I1',...]*pf@params['rho1']
        mean2 <- pred.mean(pf)['I2',...]*pf@params['rho2']
        mean3 <- pred.mean(pf)['I3',...]*pf@params['rho3']
        mean4 <- pred.mean(pf)['I4',...]*pf@params['rho4']                
        
        ylim <- range(mean1, mean2, mean3, mean4, pf@data[c("cases1", "cases2", "cases3", "cases4"),])
        plot(mean1, type='o', col=2, ylab='y1', xlab='weeks', ylim=ylim)
        lines(pf@data['cases1',], type='l')
        plot(mean2, type='o', col=2, ylab='y2', xlab='weeks', ylim=ylim)
        lines(pf@data['cases2',], type='l')
        plot(mean3, type='o', col=2, ylab='y3', xlab='weeks', ylim=ylim)
        lines(pf@data['cases3',], type='l')
        plot(mean4,type='o', col=2, ylab='y4', xlab='weeks', ylim=ylim)
        lines(pf@data['cases4',], type='l')
}
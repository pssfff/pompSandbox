## two strain pomp model in R
## Nicholas Reich
## November 2013

setwd("~/Documents/code_versioned/pompSandbox/twoStrainModel")

require(pomp)

####################
## step functions ##
####################

## step function in R

toy.proc.sim <- function(x, t, params, delta.t, ...) {
        ## unpack the params vector:
        beta1 <- params["beta1"]
        alpha1 <- params["alpha1"]
        alpha2 <- params["alpha2"]
        lambda <- params["lambda"]
        mu <- params["mu"]
        N <- params["N"]
        iota <- params["iota"]
        ## copy x for easy access
        newx <- x
        ## transitions for each strain
        newx["I1"] <- rpois(1, beta1*((x["I1"]/N+iota)^alpha1)*(x["S1"]^alpha2))
        newx["I2"] <- rpois(1, beta1*((x["I2"]/N+iota)^alpha1)*(x["S2"]^alpha2))
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

## skeleton functions in C
## SIR process model with extra-demographic stochasticity
step.fn <- '
return;
'
skel <- '
return;
'

#######################
## measure functions ##
#######################

## measure function in R
toy.meas.sim <- function(x, t, params, ...) {
        rho1 <- params["rho1"]
        rho2 <- params["rho2"]
        cases1 <- rbinom(1, size=x["I1"], prob=rho1)
        cases2 <- rbinom(1, size=x["I2"], prob=rho2)
        unname(c(cases1, cases2))
}

## measure functions for C
rmeas <- "
        cases1 = rbinom(I1, rho1);
        cases2 = rbinom(I2, rho2);
"

dmeas <- "
        double lik1 = dbinom(cases1, I1, rho1, give_log);
        double lik2 = dbinom(cases2, I2, rho2, give_log);
        lik = lik1*lik2;
"

###############################
## parameter transformations ##
###############################
 
## in R
## from estimation scale to natural scale
partransR <- function(params,...){
        exp.idx <- c("N", "mu", "beta1", "alpha1", "alpha2",
                     "iota", "lambda",
                     "S1.0", "I1.0", "C1.0",
                     "S2.0", "I2.0", "C2.0")
        params[exp.idx] <- exp(params[exp.idx])
        expit.idx <- c("rho1", "rho2")
        params[expit.idx] <- exp(params[expit.idx])/(1+exp(params[expit.idx]))
        return(params)
}
## from natural scale to estimation scale
paruntransR <- function(params,...){
        log.idx <- c("N", "mu", "beta1", "alpha1", "alpha2",
                     "iota", "lambda",
                     "S1.0", "I1.0", "C1.0",
                     "S2.0", "I2.0", "C2.0")
        params[log.idx] <- log(params[log.idx])
        logit.idx <- c("rho1", "rho2")
        params[logit.idx] <- log(params[logit.idx]/(1-params[logit.idx]))
        return(params)
}

## in C
partrans <- "
        TN = exp(N);
        Tmu = exp(mu);
        Tbeta1 = exp(beta1);
        Talpha1 = exp(alpha1);
        Talpha2 = exp(alpha2);
        Tiota = exp(iota);
        Tlambda = exp(lambda);
        TS1_0 = exp(S1_0);
        TI1_0 = exp(I1_0);
        TC1_0 = exp(C1_0);
        TS2_0 = exp(S2_0);
        TI2_0 = exp(I2_0);
        TC2_0 = exp(C2_0);
        Trho1 = expit(rho1);
        Trho2 = expit(rho2);
"
paruntrans <- "
        TN = log(N);
        Tmu = log(mu);
        Tbeta1 = log(beta1);
        Talpha1 = log(alpha1);
        Talpha2 = log(alpha2);
        Tiota = log(iota);
        Tlambda = log(lambda);
        TS1_0 = log(S1_0);
        TI1_0 = log(I1_0);
        TC1_0 = log(C1_0);
        TS2_0 = log(S2_0);
        TI2_0 = log(I2_0);
        TC2_0 = log(C2_0);
        Trho1 = logit(rho1);
        Trho2 = logit(rho2);
"

######################
## build the models ##
######################

## in R
tsirR <- pomp(
        data=data.frame(
                time=seq(0, 5000,by=1),
                cases1=NA,
                cases2=NA
        ),
        times="time",
        t0=0,
        rprocess=discrete.time.sim(
                step.fun=toy.proc.sim,
                delta.t=1
        ),
        parameter.transform = partransR,
        parameter.inv.transform = paruntransR,
        rmeasure=toy.meas.sim,
        # initial condition parameters 
        ic.pars=c("S1.0","I1.0","C1.0", 
                  "S2.0","I2.0","C2.0"), 
        # names of the compartments
        comp.names=c("S1","I1","C1", "S2","I2","C2"),
        # names of the parameters
        paramnames=c("N", "mu", "rho1", "rho2", 
                     "beta1", "alpha1", "alpha2",
                     "iota", "lambda", 
                     "S1.0", "I1.0", "C1.0",
                     "S2.0", "I2.0", "C2.0"),
        statenames=c("S1","I1","C1", "S2","I2","C2")
)

## in C
pompBuilder(
        name="TSIR",
        data=data.frame(
                time=seq(0, 5000,by=1),
                cases1=NA,
                cases2=NA
        ),
        times="time",
        t0=0,
        dmeasure=dmeas,
        rmeasure=rmeas,
        step.fn=step.fn,
        step.fn.delta.t=1,
        skeleton.type="vectorfield",
        skeleton=skel,
        tcovar="time",
        parameter.transform=partrans,
        parameter.inv.transform=paruntrans,
        ic.pars=c("S1.0","I1.0","C1.0", 
                  "S2.0","I2.0","C2.0"), 
        # names of the compartments
        comp.names=c("S1","I1","C1", "S2","I2","C2"),
        # names of the parameters
        paramnames=c("N", "mu", "rho1", "rho2", 
                     "beta1", "alpha1", "alpha2",
                     "iota", "lambda", 
                     "S1.0", "I1.0", "C1.0",
                     "S2.0", "I2.0", "C2.0"),
        statenames=c("S1","I1","C1", "S2","I2","C2")
) -> tsirC

#####################
## simulate models ##
#####################

## in R
simulate(tsirR, 
        params=c(
                N=50000,
                mu=1/500,
                rho1=.5, rho2=.1,
                beta1=1.4,alpha1=1,alpha2=1,
                iota=1/10000,
                lambda=1/100,
                S1.0=40000,I1.0=100,C1.0=100,
                S2.0=40000,I2.0=100,C2.0=100
        ),
        seed=677573454L
) -> tsirR
plot(tsirR, variables=c("cases1", "cases2", "C1","C2", "I1","I2", "S1", "S2"))

## in C
simulate(tsirC, 
         params=c(
                 N=50000,
                 mu=1/500,
                 rho1=.5, rho2=.1,
                 beta1=1.4,alpha1=1,alpha2=1,
                 iota=1/10000,
                 lambda=1/100,
                 S1.0=40000,I1.0=100,C1.0=100,
                 S2.0=40000,I2.0=100,C2.0=100
         ),
         seed=677573454L
) -> tsirC
plot(tsirC, variables=c("cases1", "cases2", "C1","C2", "I1","I2", "S1", "S2"))


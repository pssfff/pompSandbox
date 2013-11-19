## two strain pomp model in R
## Nicholas Reich
## November 2013


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
plot(tsir, variables=c("y1", "y2", "C1","C2", "I1","I2", "S1", "S2"))

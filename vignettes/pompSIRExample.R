sir.proc.sim <- function (x, t, params, delta.t, ...) {
        ## unpack the parameters
        N <- params["N"]             # population size
        gamma <- params["gamma"]     # recovery rate
        mu <- params["mu"]           # birth rate = death rate
        beta <- params["beta"]       # contact rate
        foi <- beta*x["I"]/N         # the force of infection
        trans <- c(
                rpois(n=1,lambda=mu*N*delta.t), # births are assumed to be Poisson 
                reulermultinom(n=1,size=x["S"],rate=c(foi,mu),dt=delta.t), # exits from S 
                reulermultinom(n=1,size=x["I"],rate=c(gamma,mu),dt=delta.t), # exits from I
                reulermultinom(n=1,size=x["R"],rate=c(mu),dt=delta.t) # exits from R
        )
        ## now connect the compartments
        x[c("S","I","R","cases")]+c(
                trans[1]-trans[2]-trans[3],
                trans[2]-trans[4]-trans[5],
                trans[4]-trans[6],
                trans[4]
        ) }

simulate(
        pomp(
                data=data.frame(
                        time=seq(1/52,15,by=1/52),
                        reports=NA
                ),
                times="time",
                t0=0,
                rprocess=euler.sim(
                        step.fun=sir.proc.sim,
                        delta.t=1/52/20
                ),
                measurement.model=reports~binom(size=cases,prob=rho),
                initializer=function(params, t0, ic.pars, comp.names, ...){
                        x0 <- c(S=0,I=0,R=0,cases=0)
                        N <- params["N"]
                        fracs <- params[ic.pars]
                        x0[comp.names] <- round(N*fracs/sum(fracs))
                        x0
                },
                zeronames=c("cases"), # 'cases' is an accumulator variable 
                ic.pars=c("S0","I0","R0"), # initial condition parameters 
                comp.names=c("S","I","R") # names of the compartments
        ),
        params=c(
                N=50000,
                beta=60,gamma=8,mu=1/50,
                rho=0.6,
                S0=8/60,I0=0.002,R0=1-8/60-0.001
        ),
        seed=677573454L
) -> sir
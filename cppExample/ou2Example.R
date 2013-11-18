require(pomp)
require(Rcpp)

data(ou2)
ou2.dat <- as.data.frame(ou2)

sourceCpp("ou2Example.cpp")


#######################################
## Example for non-vectorized R code ##
#######################################

ou2.Rplug <-  pomp(
        data=ou2.dat[c("time","y1","y2")],
        times="time",
        t0=0,
        rprocess=procSim_R
)

theta <- c(
        x1.0=-3, x2.0=4,
        tau=1,
        alpha.1=0.8, alpha.2=-0.5, alpha.3=0.3, alpha.4=0.9,
        sigma.1=3, sigma.2=-0.5, sigma.3=2
)

tic <- Sys.time()
simdat.Rplug <- simulate(ou2.Rplug, 
                         params=theta, 
                         states=T, 
                         nsim=1000)
toc <- Sys.time()
(etime.Rplug <- toc-tic)

#######################################
## Example for non-vectorized C code ##
#######################################
ou2.Cplug <- pomp(
        ou2.Rplug,
        rprocess=discrete.time.sim("ou2_step"),
        paramnames=c(
                "alpha.1","alpha.2","alpha.3","alpha.4",
                "sigma.1","sigma.2","sigma.3",
                "tau"
        ),
        statenames=c("x1","x2"),
        obsnames=c("y1","y2")
)

tic <- Sys.time()
simdat.Cplug <- simulate(ou2.Cplug,params=theta,states=T,nsim=100000)
toc <- Sys.time()
(etime.Cplug <- toc-tic)

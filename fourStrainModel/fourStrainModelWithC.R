## a four-strain cross-protection model with C code
## Nicholas Reich
## November 2013

setwd("~/Documents/code_versioned/pompSandbox/fourStrainModel/")
require(pomp)
source("fourStrainFuncs.R")

############################
## create base TSIR model ##
############################

## make spline basis
t.end <- 100
step <- 1/26
tbasis <- seq(0, t.end, by=step)
basis <- cbind(time=tbasis,
               as.data.frame(periodic.bspline.basis(tbasis,nbasis=3,degree=2,period=1,names="seas%d")))

## define state names for easy assignment
statenames <- c(paste0("I", 1:4), paste0("S", 1:4), paste0("C", 1:4))
ic.names <- paste0(statenames, ".0")
                
## create the TSIR pomp object 
tsirR <- pomp(
        data=data.frame(
                time=tbasis,
                cases1=NA, cases2=NA, cases3=NA, cases4=NA
        ),
        times="time",
        tcovar="time",
        covar=basis,
        t0=0,
        rprocess=discrete.time.sim(
                step.fun=fourStrainProcSim,
                delta.t=step
        ),
        rmeasure=fourStrainMeasSim,
        dmeasure=fourStrainMeasDens,
        ## from estimation scale to natural scale
        parameter.transform=partransR,
        ## from natural scale to estimation scale
        parameter.inv.transform=paruntransR,
        paramnames=c("N", "mu", 
                     "rho1", "rho2", "rho3", "rho4", 
                     "b1", "b2", "b3", 
                     "alpha1", "alpha2",
                     "iota", "lambda",
                     ic.names),
        # initial condition parameters 
        ic.pars=ic.names,
        # names of the compartments
        comp.names=statenames,
        statenames=statenames
) 

## in C
pompBuilder(
        name="TSIR_four_strain",
        data=data.frame(
                time=seq(0, t.end, by=step),
                cases1=NA,
                cases2=NA,
                cases3=NA,
                cases4=NA
        ),
        times="time",
        t0=0,
        dmeasure=dmeas,
        rmeasure=rmeas,
        step.fn=step.fn,
        step.fn.delta.t=step, ## not really treating this as a continuous system
        skeleton.type="vectorfield",
        skeleton=skel,
        tcovar="time",
        covar=basis,
        parameter.transform=partrans,
        parameter.inv.transform=paruntrans,
        ## parameter and state stuff
        paramnames=c("N", "mu", 
                     "rho1", "rho2", "rho3", "rho4", 
                     "b1", "b2", "b3", 
                     "alpha1", "alpha2",
                     "iota", "lambda",
                     ic.names),
        # initial condition parameters 
        ic.pars=ic.names,
        # names of the compartments
        comp.names=statenames,
        statenames=statenames
) -> tsirC

##########################################
## refining basic TSIR model --> model2 ##
##########################################

## first, modify the seasonal basis
## chosen to give range between ~.2 and ~1.2
basisCoefs <- c(1.1,.95,1) # basisCoefs <- c(1.5,-.1,1) 
betas <- as.matrix(basis[,2:4]) %*% basisCoefs
#plot(betas[1:52], type="l")

paramsFourStrain <- c(N=50000,
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


#####################
## simulate models ##
#####################

## in R
tic <- Sys.time()
simulate(tsirR, 
         params=paramsFourStrain,
         seed=677573454L
) -> tsirR
toc <- Sys.time()
(tictoc1 <- toc-tic)
plot(tsirR, variables=c("cases1", "cases2", "C1","C2", "I1","I2", "S1", "S2"))

## in C
tic <- Sys.time()
simulate(tsirC, 
         params=paramsFourStrain,
         seed=677573454L
) -> tsirC
toc <- Sys.time()
(tictoc2 <- toc-tic)
plot(tsirC, variables=c("cases1", "cases2", "C1","C2", "I1","I2", "S1", "S2"))
plot(tsirC, variables=c("C1","C2", "C3", "C4", "I1","I2", "I3", "I4"))

(as.numeric(tictoc1)/as.numeric(tictoc2))


###########################
## run a particle filter ##
###########################

## start with the truth, with initial conditions adjusted for windowing   
dur <- 10 ## duration of time series
index0 <- max(which(tsirC@times<t.end-dur))## variables for timezero
theta.truth <- paramsFourStrain     
theta.truth[ic.names] <- states(tsirC)[statenames,index0]


## in C
tsirC_short <- window(tsirC, start=t.end-dur, end=t.end)
tic <- Sys.time()
pfC <- pfilter(tsirC_short, params=theta.truth, Np=100, max.fail=length(tsirC_short@times)+1)
toc <- Sys.time()
(tictoc.pfC <- toc-tic)
print(round(logLik(pfC),1))

## in R
tsirR_short <- window(tsirR, start=t.end-dur, end=t.end)
tic <- Sys.time()
pfR <- pfilter(tsirR_short, params=theta.truth, Np=100, max.fail=length(tsirC_short@times)+1)
toc <- Sys.time()
(tictoc.pfR <- toc-tic)
print(round(logLik(pfR),1))

(as.numeric(tictoc.pfR)/as.numeric(tictoc.pfC))

#######################
## miffing the data  ##
#######################

nmif <- 10
## parallelizing
require(doMC)
registerDoMC(10)
estpars <- c("lambda")
tic <- Sys.time()
mf <- foreach(i=1:nmif) %dopar% {
        theta.guess <- theta.truth
        theta.guess[estpars] <- rlnorm(
                n=length(estpars),
                meanlog=log(theta.guess[estpars]),
                sdlog=.1
        )
        mif(
                tsirC_short,
                Nmif=100,
                start=theta.guess,
                transform=TRUE,
                pars=estpars,
                rw.sd=c(lambda=0.02),
                Np=2000,
                var.factor=4,
                ic.lag=10,
                cooling.factor=0.999,
                max.fail=length(tsirC_short@times)+1
        ) 
}
toc <- Sys.time()
(toc-tic)
compare.mif(mf)

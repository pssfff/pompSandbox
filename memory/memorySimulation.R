## Memory cross-protection analysis simulation
## Nicholas Reich
## May 2014

## preliminaries
setwd('~/Documents/code_versioned/pompSandbox/memory/')
source('~/Documents/code_versioned/pompSandbox/memory/fourStrainMemory.R')
path <- "~/Documents/work/dengueData/simdata_5_local/"
key <- read.csv("~/Documents/work/research/dengueCrossProtection/data/simulatedData_fromMI/key_5.csv")
colnames(key)[1] <- "filename"

## parameters
k <- 300
lambda <- 26*4
file_idx <- 1:12000 #c(701:750, 4681:4730, 2911:2934, 2936:2961)
n_files <- length(file_idx)

## storage
results <- matrix(NA, ncol=9, nrow=n_files)
colnames(results) <- c(colnames(key), "k", "lambda", "loglik", "z", "beta_M")
results[,1] <- file_idx
results[,2:4] <- as.matrix(key[file_idx,2:4])

## create loop for each filename
require(doMC)
registerDoMC(18)
tic <- Sys.time()
parResults <- foreach(i = 1:n_files, .combine=rbind) %dopar% {
        fn <- paste0(path, key[file_idx[i],"filename"])
        raw_data <- read.csv(fn, row.names=1)
        tsir_data <- raw_data[,2:5]
        runCrossProtectMemoryAnalysis(data=tsir_data, k=k, 
                                      subset=500:1041,
                                      max_lambda=lambda, 
                                      plot=FALSE, verbose=FALSE)
}
toc <- Sys.time()
toc-tic
results[,5:9] <- parResults

write.csv(results, file="simResults.csv", quote=FALSE, row.names=FALSE)

#results <- read.csv("simResults.csv")

results <- as.data.frame(results)
results$trueLambda <- 26/results$delta
results$memoryIndicator <- results$lambda>13# abs(results$beta_M)>.02
#results
qplot(factor(trueLambda), lambda, data=results, geom="violin")
qplot(factor(trueLambda), lambda, data=results, geom="boxplot", facets=.~rho)
qplot(factor(trueLambda), z, data=results, geom="violin")
qplot(factor(trueLambda), z, data=results, geom="boxplot")
qplot(factor(trueLambda), sign(z)*lambda, data=results, geom="boxplot")
qplot(factor(trueLambda), beta_M, data=results, geom="boxplot")
qplot(factor(trueLambda), beta_M, data=results, geom="violin")
qplot(factor(memoryIndicator), data=results, geom="bar", fill=factor(trueLambda))
qplot(exp(beta_M), lambda, data=results, alpha=I(.2), color=factor(trueLambda), facets=chi~trueLambda) 

qplot(exp(beta_M), lambda, data=results, alpha=I(.2), color=factor(trueLambda), facets=rho+chi~trueLambda) 

qplot(exp(beta_M), lambda, data=results, alpha=I(.2), color=factor(trueLambda), facets=rho~trueLambda) 


qplot(exp(beta_M), lambda, data=results, alpha=I(.4), color=factor(trueLambda), facets=.~trueLambda) + xlim(exp(-.1), exp(.1)) 


with(results, table(trueLambda, memoryIndicator))/2000

## tree exploration
treeFit <- rpart(CP ~ lambda + exp(beta_M), data=results)
plot(treeFit)
text(treeFit, use.n=TRUE)

## run Dengue analysis
load("/Users/nick/Dropbox/work/research/dengueCrossProtection/data/bkk.dengue.cases.new.rda")
dengue_data <- bkk.dengue.cases[,3:6]
dengue_data[which(dengue_data<0, arr.ind=TRUE)] <- 0 ## remove negative values
tmp <- runCrossProtectMemoryAnalysis(data=dengue_data, k=300, 
                                     max_lambda=26*6, 
                                     plot=TRUE, verbose=FALSE)
exp(tmp["beta_M"])

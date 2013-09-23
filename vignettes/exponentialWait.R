

a <- round(100*dnorm(seq(-5, 5, .01)))
c <- rep(NA, length(a))

lambda <- 1/1000

c[1] <- 0
for(i in 2:length(a)) {
        outgoing <- rpois(1, c[i-1]*lambda)
        c[i] <- c[i-1] + a[i] - outgoing
}

par(mfrow=c(2,1))
plot(a, type="l", lty=2)
abline(v=500)
plot(c, type="l")
abline(v=c(500, 500+1/lambda)) ## line at middle and at middle + exponential mean
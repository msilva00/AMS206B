rm(list=ls(all=TRUE))
set.seed(454)

source("fn-density-1.R")

den <- read.csv("pdensity.txt", header = T, sep = "")
N <- nrow(den)

pdf("dat.pdf")
plot(den[, 2], den[, 3])

par(mfrow=c(2,2))
for(i in 1:10)
{
    ind <- (1:N)[den[,1]==i]
    plot(den[ind, 2], den[ind, 3])
}
dev.off()

### preparing data
mydat <- NULL
mydat$X <- as.matrix(cbind(den$density, (den$density)^2))
mydat$y <- den$yield
mydat$I <- 10
mydat$J <- 8
mydat$plot <- den$plot
mydat$IJ <- mydat$I*mydat$J
mydat$sum_x2 <-  tapply(mydat$X[,2], mydat$plot, sum)
mydat$sum_x4 <-  tapply((mydat$X[,2])^2, mydat$plot, sum)

#### fixed hyperparameters
hyper <- NULL
hyper$a_sig <- 0.1
hyper$b_sig <- 0.1

hyper$a_t0 <- 0.1
hyper$b_t0 <- 0.1
hyper$v2_0 <- 100

hyper$a_t1 <- 0.1
hyper$b_t1 <- 0.1
hyper$v2_1 <- 100

hyper$a_t2 <- 0.1
hyper$b_t2 <- 0.1
hyper$v2_2 <- 100

### initialize the chain - used OLS estimates
cur_sam <- NULL
cur_sam$sig2 <- (0.9822)^2
cur_sam$b0 <- rep(2.86875, mydat$I)
cur_sam$b1 <- rep(1.85485, mydat$I)
cur_sam$b2 <- rep(-0.15925, mydat$I)
cur_sam$mu0 <- 0
cur_sam$mu1 <- 0
cur_sam$mu2 <- 0
cur_sam$tau2_0 <- 1
cur_sam$tau2_1 <- 1
cur_sam$tau2_2 <- 1

### MCMC parameters
N_burn <- 3000
N_sam <- 5000
MCMC_sam <- NULL
MCMC_sam$sig2 <- MCMC_sam$mu1 <- MCMC_sam$mu2 <- MCMC_sam$tau2_0 <- MCMC_sam$tau2_1 <- MCMC_sam$tau2_2 <- MCMC_sam$mu0 <- rep(0, N_sam)
MCMC_sam$b0 <- MCMC_sam$b1 <- MCMC_sam$b2 <- array(NA, dim=c(N_sam, mydat$I))


### for burn-in
for(i_iter in 1:N_burn)
{
    cur_sam$sig2 <- fn.update.sig2(mydat, cur_sam$b0, cur_sam$b1, cur_sam$b2, hyper$a_sig, hyper$b_sig)
    
    cur_sam$b0 <- fn.update.b0(mydat, cur_sam$sig2, cur_sam$tau2_0, cur_sam$mu0, cur_sam$b1, cur_sam$b2)
    cur_sam$mu0 <- fn.update.mu0(mydat$I, cur_sam$tau2_0, hyper$v2_0, cur_sam$b0)
    cur_sam$tau2_0 <- fn.update.tau2.0(hyper$a_t0, hyper$b_t0, cur_sam$b0, cur_sam$mu0, mydat$I)
    
    cur_sam$b1 <- fn.update.b1(mydat, cur_sam$sig2, cur_sam$tau2_1, cur_sam$mu1, cur_sam$b0, cur_sam$b2)
    cur_sam$mu1 <- fn.update.mu0(mydat$I, cur_sam$tau2_1, hyper$v2_1, cur_sam$b1)
    cur_sam$tau2_1 <- fn.update.tau2.0(hyper$a_t1, hyper$b_t1, cur_sam$b1, cur_sam$mu1, mydat$I)
    
    cur_sam$b2 <- fn.update.b2(mydat, cur_sam$sig2, cur_sam$tau2_2, cur_sam$mu2, cur_sam$b0, cur_sam$b1)
    cur_sam$mu2 <- fn.update.mu0(mydat$I, cur_sam$tau2_2, hyper$v2_2, cur_sam$b2)
    cur_sam$tau2_2 <- fn.update.tau2.0(hyper$a_t2, hyper$b_t2, cur_sam$b2, cur_sam$mu2, mydat$I)
} ##for(i_iter in 1:N_burn)


### after burn-in
for(i_iter in 1:N_sam)
{
    cur_sam$sig2 <- fn.update.sig2(mydat, cur_sam$b0, cur_sam$b1, cur_sam$b2, hyper$a_sig, hyper$b_sig)
    
    cur_sam$b0 <- fn.update.b0(mydat, cur_sam$sig2, cur_sam$tau2_0, cur_sam$mu0, cur_sam$b1, cur_sam$b2)
    cur_sam$mu0 <- fn.update.mu0(mydat$I, cur_sam$tau2_0, hyper$v2_0, cur_sam$b0)
    cur_sam$tau2_0 <- fn.update.tau2.0(hyper$a_t0, hyper$b_t0, cur_sam$b0, cur_sam$mu0, mydat$I)
    
    
    cur_sam$b1 <- fn.update.b1(mydat, cur_sam$sig2, cur_sam$tau2_1, cur_sam$mu1, cur_sam$b0, cur_sam$b2)
    cur_sam$mu1 <- fn.update.mu0(mydat$I, cur_sam$tau2_1, hyper$v2_1, cur_sam$b1)
    cur_sam$tau2_1 <- fn.update.tau2.0(hyper$a_t1, hyper$b_t1, cur_sam$b1, cur_sam$mu1, mydat$I)
    
    cur_sam$b2 <- fn.update.b2(mydat, cur_sam$sig2, cur_sam$tau2_2, cur_sam$mu2, cur_sam$b0, cur_sam$b1)
    cur_sam$mu2 <- fn.update.mu0(mydat$I, cur_sam$tau2_2, hyper$v2_2, cur_sam$b2)
    cur_sam$tau2_2 <- fn.update.tau2.0(hyper$a_t2, hyper$b_t2, cur_sam$b2, cur_sam$mu2, mydat$I)
    
    
    MCMC_sam$sig2[i_iter] <- cur_sam$sig2
    MCMC_sam$b0[i_iter, ] <- cur_sam$b0
    MCMC_sam$tau2_0[i_iter] <- cur_sam$tau2_0
    MCMC_sam$mu0[i_iter] <- cur_sam$mu0
    
    MCMC_sam$b1[i_iter, ] <- cur_sam$b1
    MCMC_sam$tau2_1[i_iter] <- cur_sam$tau2_1
    MCMC_sam$mu1[i_iter] <- cur_sam$mu1
    
    
    MCMC_sam$b2[i_iter, ] <- cur_sam$b2
    MCMC_sam$tau2_2[i_iter] <- cur_sam$tau2_2
    MCMC_sam$mu2[i_iter] <- cur_sam$mu2
} ## for(i_iter in 1:N_sam)



pdf("para-1.pdf")
par(mfrow=c(2,2))
hist(MCMC_sam$sig2)

par(mfrow=c(2,2))
hist(MCMC_sam$mu0)
hist(MCMC_sam$tau2_0)

for(i in 1:mydat$I)
{
    hist(MCMC_sam$b0[,i], main=i)
}

par(mfrow=c(2,2))
hist(MCMC_sam$mu1)
hist(MCMC_sam$tau2_1)

for(i in 1:mydat$I)
{
    hist(MCMC_sam$b1[,i], main=i)
}

par(mfrow=c(2,2))
hist(MCMC_sam$mu2)
hist(MCMC_sam$tau2_2)

for(i in 1:mydat$I)
{
    hist(MCMC_sam$b2[,i], main=i)
}


dev.off()

x <- 2

pdf("pred-2-1.pdf")
par(mfrow=c(2,2))
for(i in 1:mydat$I)
{
    y1 <- rnorm(N_sam, MCMC_sam$b0[,i] + MCMC_sam$b1[,i]*x + MCMC_sam$b2[,i]*(x^2), sqrt(MCMC_sam$sig2))
    hist(y1)
    qq <- mydat$y[(mydat$X[,1]==x)&(mydat$plot==i)]
    points(qq, c(0, length(qq)), cex=4, pch=4, col=2, lwd=3)
}
dev.off()





b0_hat <- apply(MCMC_sam$b0, 2, mean)
b0_1 <- apply(MCMC_sam$b0, 2, quantile, 0.025)
b0_2 <- apply(MCMC_sam$b0, 2, quantile, 0.975)

cbind(b0_1, b0_hat, b0_2)


tmp <- lm(mydat$y ~ 1 + mydat$X)

tmp <- lm(den$yield ~ 1 + den$density + I((den$density)^2))



library(lme4)

tmp1 <- lme(mydat$y ~ mydat$X, random=~1|mydat$plot)



tmp1 <- lmer(mydat$y ~ (1|mydat$plot) + mydat$X)


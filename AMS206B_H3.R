alpha = 1/2
beta = 1/2
a = alpha + x
b = beta + n - x
#### Problem 6a ####
x = 1
n =10
a 
lower = qbeta(0.025, a, b)
upper = qbeta(0.975, a, b)
print(c(lower, upper))

#### Problem 6b ####
# calculate the MAP
theta_hat <- (a-1)/(a+b-2)
# evaluate the first and second derivaties at the MAP
h <- (a-1)*log(theta_hat) + (b-1)*log(1-theta_hat) 
h.2 <- -(a-1)/theta_hat^2 - (b-1)/(1-theta_hat)^2
#calculate the constant
Const <- exp(h)*sqrt(2*pi/(-h.2))/beta(a, b)
# compute the interval
lower<-theta_hat + qnorm((1 - 0.95/Const)/2, 0, 1)/sqrt(-h.2)
upper<-theta_hat + qnorm((1 - 0.95/Const)/2, 0, 1, lower.tail=FALSE)/sqrt(-h.2)
print(c(lower,upper))

#### Problem 6c ####
th <- rbeta(10000, a, b) 
quantile(th, probs=c(0.025, 0.975))

#### Problem 6d ####
n = c(10, 100, 1000)
x = c(1, 10, 100)
a = alpha + x
b = beta + n - x
lower = qbeta(0.025, a, b)
upper = qbeta(0.975, a, b)
exact = cbind(lower, upper)
print(exact)
# calculate the MAP
theta_hat <- (a-1)/(a+b-2)
# evaluate the first and second derivaties at the MAP
h <- (a-1)*log(theta_hat) + (b-1)*log(1-theta_hat) 
h.2 <- -(a-1)/theta_hat^2 - (b-1)/(1-theta_hat)^2
#calculate the constant
Const <- exp(h)*sqrt(2*pi/(-h.2))/beta(a, b)
# compute the interval
lower<-theta_hat + qnorm((1 - 0.95/Const)/2, 0, 1)/sqrt(-h.2)
upper<-theta_hat + qnorm((1 - 0.95/Const)/2, 0, 1, lower.tail=FALSE)/sqrt(-h.2)
laplace = (cbind(lower,upper))
print(laplace)


th1 <- rbeta(10000, a[1], b[1]) 
q1 = quantile(th1, probs=c(0.025, 0.975))
th2 <- rbeta(10000, a[2], b[2]) 
q2 = quantile(th2, probs=c(0.025, 0.975))
th3 <- rbeta(10000, a[3], b[3]) 
q3 = quantile(th3, probs=c(0.025, 0.975))

rbind(q1,q2,q3)

#### Problem 8e ####
# set the number of observations and true value of parameters
n <- 1000
tr.th <- 5
tr.sig2 <- 1
# generate dataset
x <- rnorm(n, tr.th, sqrt(tr.sig2)) #set the number of MC samples
N.sam <- 5000
#i. set hyperparameters for fairly informative priors
th0 <- tr.th
k0 <- 0.01
a <- 1001
b <- tr.sig2/(a-1)
#calculate posterior parameters
m <- (th0 + n*k0*mean(x))/(1+n*k0)
alpha <- a + n/2
beta <- 1/b + sum(x^2)/2 + th0^2/(2*k0) - (th0 + n*k0*mean(x))^2/(2*k0*(1+n*k0))
#MC simulation
sig2.1 <- 1/rgamma(N.sam, alpha, beta)
th.1 <- rnorm(N.sam, m, sqrt(sig2.1/(1/k0 + n)))
par(mar=c(4.5, 4.5, 2.1, 2.1), mfrow=c(1,2))
hist(th.1, col=8, lwd=2, , main="", cex.axis=1.5, cex.lab=1.5)
abline(v=tr.th, lty=2, lwd=3, col=2)
hist(sig2.1, col=8, lwd=2, , main="", cex.axis=1.5, cex.lab=1.5) 
abline(v=tr.sig2, lty=2, lwd=3, col=2)

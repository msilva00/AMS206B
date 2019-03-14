# install.packages("coda")
# install.packages("mvtnorm")
# install.packages("MASS")

#### 1A ####
theta_1 = 1.5 # true value theta1
theta_2 = 2 # true value theta2
mean_z1 = sqrt(theta_2/theta_1)
mean_z2 = sqrt(theta_1/theta_2) + 1/(2*theta_2)

# hyperparams
b = 100
a = mean_z1*b


pdf_z = function(z){
  -(3/2)*log(z) - theta_1*z - theta_2/z
}



#M-H Algorithm
MH_alg1 = function(N){
  MH_samples = rep(NA, N)
  count = 0
  current_z = 1.0
  for(i in 1:N){
    curr_p = pdf_z(current_z) 
    z_new = rgamma(1, a, b)
    p_new = pdf_z(z_new)
    
    accept = exp(p_new + dgamma(current_z,a,b,log = T) - p_new - dgamma(z_new,a,b,log = T))
    if(runif(1) < accept){
      current_z = z_new
      count = count + 1
    }
    MH_samples[i] = current_z
  }
  return(list(MH_samples=MH_samples,count=count))
}
alg1 = MH_alg1(1000)
MH_samples = alg1$MH_samples
count = alg1$count
plot.ts(MH_samples)
mean(MH_samples)
mean(1/MH_samples)
count/length(MH_samples)
mean1_1 = mean(MH_samples)
mean1_2 = mean(1/MH_samples)
acceptance_rate1 = count/length(MH_samples)

#### 1B ####
v = 0.01
pdf_z2 = function(z){
  return(-(1/2)*log(z) - theta_1*z - theta_2/z)
}
MH_RW = function(N){
  N = N
  MH_RW = rep(NA, N)
  a_count = 0
  z_curr = 1.0
  for (i in 1:N) {
    p_curr = pdf_z2(z_curr)
    z_new = exp(log(z_curr) + rnorm(1,0,sqrt(v)))
    p_new = pdf_z2(z_new)
    acceptance = exp(p_new - p_curr)
    if(runif(1) < acceptance){
      z_curr = z_new
      a_count = a_count+1
    }
    MH_RW[i] = z_curr
  }
  return(list(MH_RW=MH_RW, a_count=a_count))
}
alg2 = MH_RW(50000)
MH_RW = alg2$MH_RW
count = alg2$a_count
plot.ts(MH_RW)
(mean2_1 = mean(MH_RW))
(mean2_2=mean(1/MH_RW))
(count2 = count/length(MH_RW))

#### 2A ####
x = read.table("my-data.txt", header = F)[,1]
n = length(x)
sum_x = sum(x)
sum_logx = sum(log(x))
library(coda)
library(mvtnorm)
library(MASS)

nu_condit = function(nu, theta_curr){
  return(nu*(sum_logx + n*log(theta_curr)-1)-n*lgamma(nu)+3*log(nu))
}


sample = NULL
N = 1000
sample$theta = rep(NA,N)
sample$nu = rep(NA,N)
alpha = 1
beta = 1
v = 0.01
theta_curr = 2
nu_curr = 3
set.seed(1)
for(i in 1:N){
  theta_curr = rgamma(1, n*nu_curr + alpha, beta + sum_x)
  nu_new = exp(log(nu_curr) + rnorm(1,0,sqrt(v)))
  pnu_curr = nu_condit(nu_curr, theta_curr)
  pnu_new = nu_condit(nu_new, theta_curr)
  accept = exp(pnu_new - pnu_curr)
  if(runif(1) < accept)
    nu_curr = nu_new
  
  sample$theta[i] = theta_curr
  sample$nu[i] = nu_new
}

effectiveSize(sample$theta)
effectiveSize(sample$nu)


thetas = sample$theta
nus = sample$nu 
mean(nus)
mean(thetas)

plot.ts(thetas, ylab = expression(theta), xlab = "", main = "Trace")
plot.ts(nus, ylab = expression(nu), xlab = "", main = "Trace")

quantile(nus, c(0.025,0.975))
quantile(thetas, c(0.025,0.975))

pcurr = function(nu_curr, theta_curr){
  return((n+nu_curr+2)*log(theta_curr)-n*lgamma(nu_curr) + nu_curr*(sum_logx-1) + 3*log(nu_curr) - theta_curr*(2+sum_x))
}
#### 2B ####
set.seed(1)
V = 0.05*diag(2)
theta_curr = 2
nu_curr = 3
N_test = 500
for(i in 1:N_test){
  nu_new = exp(log(nu_curr) + rnorm(1,0,sqrt(V[1,1])))
  theta_new = exp(log(theta_curr) + rnorm(1,0, V[2,2]))
  p_curr = pcurr(nu_curr, theta_curr)
  p_new = pcurr(nu_new, theta_new)
  accept = exp(p_new - p_curr)
  if(runif(1) < accept){
    nu_curr = nu_new
    theta_curr = theta_new
  }
  sample$theta[i] = theta_curr
  sample$nu[i] = nu_curr
}
V = cov(cbind(sample$nu, sample$theta))


pcurr2 = function(nu,theta){
  return(n*nu*log(theta) - n*lgamma(nu) + nu *(sum_logx -1) + 3*log(nu) + log(theta) - theta * (2 + sum_x))
}
for(i in N_test+1:N){
  new = mvrnorm(1, c(log(nu_curr), log(theta_curr)), V)
  nu_new = exp(new[1])
  theta_new = exp(new[2])
  p_curr = pcurr2(nu_curr, theta_curr)
  p_new = pcurr2(nu_new, theta_new)
  acceptance = exp(p_new - p_curr)
  if(runif(1) < acceptance){
    nu_curr = nu_new
    theta_curr = theta_new
  }
  sample$theta[i] = theta_curr
  sample$nu[i] = nu_curr
}

mean(sample$nu)
mean(sample$theta)
nu2 = sample$nu
theta2 = sample$theta
quantile(nu2, probs = c(0.025,0.975))
quantile(theta2, probs = c(0.025,0.975))
plot.ts(nu2[0:5000])
plot.ts(theta2[5000:10000])
autocorr.plot(nu2, auto.layout = F)


#### 2C ####
h = function(w) {
  a1 = (exp(w[2]) - 1) * sum.log.x + 3 * w[2] - exp(w[2])
  return(-(a1 - n * lgamma(exp(w[2])) + (2 + n * exp(w[2])) * w[1] - exp(w[1]) * (2 + sum.x) ))
}


for(i in N_test+1:N){
  nu_new = exp(log(nu_curr) + rnorm(1,0,sqrt(V[1,1])))
  theta_new = exp(log(theta_curr) + rnorm(1,0,sqrt(V[2,2])))
  p_curr = pcurr(nu_curr = nu_curr, theta_curr = theta_curr)
  p_new =pcurr(nu_curr = nu_new, theta_curr = theta_new)
  
  accept = exp(p_new - p_curr)
  if(runif(1) < accept){
    nu_curr = nu_new
    theta_curr = theta_new
  }
  sample$theta[i] = theta_curr
  sample$nu[i] = nu_curr
}

#### 4 ####
set.seed(2)
# input data
y <-c(4,5,4,1,0,4,3,4,0,6,3,3,
      4,0,2,6,3,3,5,4,5,3,1,4,
      4,1,5,5,3,4,2,5,2,2,3,4,
      2,1,3,2,2,1, 1,1,1,3,0,0,
      1,0,1,1,0,0,3,1,0,3,2,2,0,
      1,1,1,0,1,0,1,0,0,0,2,1,0,
      0,0,1,1,0,2,3,3,1,1,2,1,
      1, 1,1,2,4,2,0,0,0,1,4,0,
      0,0,1,0,0,0,0,0,1,0,0,1,0,1)
n <- length(y)
# specify hyperparameters. In this case using the data.
alpha <- 3
beta <- alpha/mean(y[1:40]) 
gam <- 3
delta <- gam/mean(y[-(1:40)])

# set up MCMC variables
N <- 50000
N.burn <- 5000
sample_save <- NULL
sample_save$m <- rep(NA, N) 
sample_save$theta <- rep(NA, N) 
sample_save$phi <- rep(NA, N)

# initialize chains
theta_curr <- rgamma(1, alpha, beta) 
phi_curr <- rgamma(1, gam, delta) 
m_curr <- 40
# sampling
for(i in 1:N){
  theta_curr <- rgamma(1, sum(y[1:m_curr]) + alpha, m_curr + beta)
  
  phi_curr <- rgamma(1, sum(y[-(1:m_curr)]) + gam, (n-m_curr + delta))
  
  m_new <- sample((1:n), 1, FALSE)
  
  p_curr <- lgamma(sum(y[1:m_curr]) + alpha) - (sum(y[1:m_curr]) + alpha)*log(m_curr + beta) +lgamma(sum(y[-(1:m_curr)]) + gam) - (sum(y[-(1:m_curr)]) + gam)*log(n-m_curr + delta)
  
  p_new <- lgamma(sum(y[1:m_new]) + alpha) - (sum(y[1:m_new]) + alpha)*log(m_new + beta) +lgamma(sum(y[-(1:m_new)]) + gam) - (sum(y[-(1:m_new)]) + gam)*log(n-m_new + delta)
  # calculate acceptance probability and accept/reject acordingly
  accpt.prob <- exp(p_new - p_curr) 
  if(runif(1) < accpt.prob)
  {
    m_curr <- m_new
  }
  # save the current draws
  sample_save$theta[i] <- theta_curr
  sample_save$phi[i] <- phi_curr
  sample_save$m[i] <- m_curr
}

# plots for phi
plot.ts(tail(sample_save$phi,5000), main = "Traceplot", ylab = expression(phi), xlab = "")
abline(h = mean(tail(sample_save$phi,5000)), col = "red")

hist(tail(sample_save$phi,5000), xlab = expression(phi), main = expression("Histogram for " ~ phi))
abline(v = mean(tail(sample_save$phi,5000)), col = "red")

# plots for theta 
plot.ts(tail(sample_save$theta,5000), main = "Traceplot", ylab = expression(theta), xlab = "")
abline(h = mean(tail(sample_save$theta,5000)), col = "red")
hist(tail(sample_save$theta,5000), xlab = expression(theta), main = expression("Histogram for " ~ theta))
abline(v = mean(tail(sample_save$theta,5000)), col = "red")

# plots for m
plot.ts(tail(sample_save$m,5000), main = "Traceplot", ylab = expression(m), xlab = "")
abline(h = mean(tail(sample_save$m,5000)), col = "red")

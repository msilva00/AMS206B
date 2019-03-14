install.packages("coda")
install.packages("mvtnorm")
install.packages("MASS")

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
N = 10000
sample$theta = rep(NA,N)
sample$nu = rep(NA,N)
alpha = 3
beta = 2
v = 0.01
theta_curr = 2
nu_curr = 3

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


thetas = sample$theta
nus = sample$nu 
mean(nus)
mean(thetas)

plot.ts(thetas)
plot.ts(nus)

pcurr = function(nu_curr, theta_curr){
  return((n+nu_curr+2)*log(theta_curr)-n*lgamma(nu_curr) + nu_curr*(sum_logx-1) + 3*log(nu_curr) - theta_curr*(2+sum_x))
}
#### 2B ####
V = 0.05*diag(2)
theta_curr = 2
nu_curr = 3
N_test = 5000
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

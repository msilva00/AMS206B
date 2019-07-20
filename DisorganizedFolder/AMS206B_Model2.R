setwd("~/Desktop")
pdensity_dat = read.csv("pdensity.txt", header = T, sep = "")

names(pdensity_dat) = c("Plot", "Density", "Yield")
head(pdensity_dat)

splitByPlot <- split(pdensity_dat, pdensity_dat$Plot)
data_by_machine = split(Machine[, "measurements"], Machine[,"machine"])

sum_x_counts_by_plot = lapply(splitByPlot, function(x){c(sum(x))})

colSums(splitByPlot[[1]])


# declare hyper-params
a_sigma = 100
b_sigma = 100
c_0 = 100
c_1 = 100
c_2 = 100
d_0 = 100
d_1 = 100
d_2 = 100
v_0 = 100
v_1 = 100
v_2 = 100
v_ell = c(100,100,100)
c_ell = c(100,100,100)
d_ell = c(100,100,100)

N = 3000
# B0i = rep(NA,N)
# B0i[1] = 1
# B1i = rep(NA,N)
# B1i[1] = 1
# B2i = rep(NA,N)
# B2i[1] = 1
# mu0 = rep(NA,N)
# mu0[1] = 1
# mu1 = rep(NA,N)
# mu1[1] = 1
# mu2 = rep(NA,N)
# mu2[1] = 1
# tao0 = rep(NA,N)
# tao0[1] = 1
# tao1 = rep(NA,N)
# tao1[1] = 1
# tao2 = rep(NA,N)
# tao2[1] = 1
# sig2 = rep(NA,N)
# sig2[1] = 1

B0i = rep(1,N)
B1i = rep(1,N)
B2i = rep(1,N)
mu0 = rep(1,N)
mu1 = rep(1,N)
mu2 = rep(1,N)
tao0 = rep(1,N)
tao1 = rep(1,N)
tao2 = rep(1,N)
sig2 = rep(1,N)


Xij = pdensity_dat$Density
Yij = pdensity_dat$Yield

I=10
J = 8
B0i_mean = function(i){
  num = sig2[i-1]* mu0[i-1] + tao0[i-1] * sum(Yij - (B1i[i-1]*Xij + B2i[i-1] * Xij^2))
  den = J*tao0[i-1] + sig2[i-1]
  return(num/den)
}
B0i_sd = function(i){
  num = tao0[i-1]*sig2[i-1]
  den = J*tao0[i-1] + sig2[i-1]
  return(sqrt(num/den))
}
i=2
sum(Xij*(Yij - (B0i[i-1] - B2i[i-1]*Xij^2 ) ))
B1i_mean = function(i){
  num = tao1[i-1] * sum(Xij*(Yij - (B0i[i-1] - B2i[i-1]*Xij^2 ) )) + sig2[i-1]*mu1[i-1]
  den = tao1[i-1]* sum(Xij^2) + sig2[i-1]
  return(num/den)
}
B1i_sd = function(i){
  num = tao1[i-1] * sig2[i-1]
  den = tao1[i-1]* sum(Xij^2) + sig2[i-1]
  sqrt(num/den)
}
i
B2i_mean = function(i){
  num = sig2[i-1]*mu2[i-1] + tao2[i-1] *sum(Xij^2 * (Yij - (B0i[i-1] + B1i[i-1]*Xij) ))
  den = tao2[i-1] * sum(Xij^4) + sig2[i-1]
  return(num/den)
}
B2i_sd = function(i){
  num = tao2[i-1]* sig2[i-1]
  den = tao2[i-1]* sum(Xij^4) + sig2[i-1]
  return(sqrt(num/den))
}
ell = 0

mu_mean = function(i, ell){
  if(ell == 0){
    B_ell = B0i
    tao = tao0
  }
  if(ell == 1){
    B_ell = B1i
    tao = tao1
  }
  if(ell == 2){
    B_ell = B2i
    tao = tao2}
  num = v_ell[ell+1] * sum(B_ell[1])
  den = I*v_ell[ell+1] + tao[1]
  return(num/den)
}
mu_sd = function(i, ell){
  if(ell == 0){
    B_ell = B0i
    tao = tao0
  }
  if(ell == 1){
    B_ell = B1i
    tao = tao1
  }
  if(ell == 2){
    B_ell = B2i
    tao = tao2}
  num = tao[1]*v_ell[ell+1]
  den = I*v_ell[ell+1] + tao[1]
  return(sqrt(num/den))
}
tao_shape = function(i, ell){
  if(ell == 0){
    B_ell = B0i
    mu = mu0
  }
  if(ell == 1){
    B_ell = B1i
    mu = mu1
  }
  if(ell == 2){
    B_ell = B2i
    mu = mu2
  }
  return((I/2) + c_ell[ell+1])
}
i
tao_rate = function(i, ell){
  if(ell == 0){
    B_ell = B0i
    mu = mu0
  }
  if(ell == 1){
    B_ell = B1i
    mu = mu1
  }
  if(ell == 2){
    B_ell = B2i
    mu = mu2
  }
  return(d_ell[ell+1] + (1/2)*sum( B_ell[1] - mu[1] )^2 )
}
  
sig2_rate = function(i){
  return(b_sigma + (1/2) * 
           rnorm(1,0,sig2[1]))
}

#### Gibbs Model II ####
for(i in 2:N){
  # update B0i
  B0i[i] = rnorm(1, B0i_mean(i), sqrt(B0i_sd(i)))
  # update B1i
  B1i[i] = rnorm(1, B1i_mean(i), B1i_sd(i))
  # update B2i
  B2i[i] = rnorm(1, B2i_mean(i), B2i_sd(i))
  # update mu0
  mu0[i] = rnorm(1, mu_mean(i, ell = 0), mu_sd(i, ell = 0))
  # update mu1
  mu1[i] = rnorm(1, mu_mean(i, ell = 1), mu_sd(i, ell = 1))
  # update mu2
  mu2[i] = rnorm(1, mu_mean(i, ell = 2), mu_sd(i, ell = 2))
  # update tao0
  tao0[i] = rinvgamma(1,tao_shape(i, ell = 0), scale = abs(tao_rate(i, ell = 0)))
  # update tao1
  tao1[i] = rinvgamma(1,tao_shape(i, ell = 1), scale=abs(tao_rate(i, ell = 1)))
  # update tao2
  tao2[i] = rinvgamma(1,tao_shape(i, ell = 2), scale=abs(tao_rate(i, ell = 2)))
  # update sig2
  sig2[i] = rinvgamma(1, 50, abs(rnorm(1,10,100)))
}
par(mfrow = c(1,2))
plot.ts(B0i[1500:3000]+3, ylab = "B0i")
plot.ts(B1i[1500:3000]+3, ylab = "B1i")
plot.ts(B2i[1500:3000]+3, ylab = "B2i")
plot.ts(mu0[1500:3000]+3)
plot.ts(mu1[1500:3000])
plot.ts(mu2[1500:3000])
plot.ts(tao0[1500:3000])
plot.ts(tao1[1500:3000])
plot.ts(tao2[1500:3000])
plot.ts(sig2[1500:3000])

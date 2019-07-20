setwd("~/Desktop")
pdensity_dat = read.csv("pdensity.txt", header = T, sep = "")

names(pdensity_dat) = c("Plot", "Density", "Yield")
head(pdensity_dat)


# Y_ij
Dat_byPlot = unstack(pdensity_dat, Yield ~ Plot)
names(Dat_byPlot) = c("Y1j","Y2j","Y3j","Y4j","Y5j","Y6j","Y7j","Y8j","Y9j","Y10j")
Yij_mat = as.matrix(Dat_byPlot)
colnames(Yij_mat) = NULL

# X_ij
Dat_byDensity = unstack(pdensity_dat, Density ~ Plot)
names(Dat_byDensity) = c("X1j","X2j","X3j","X4j","X5j","X6j","X7j","X8j","X9j","X10j")
Xij_mat = as.matrix(Dat_byDensity)
dim(Xij_mat)
colnames(Xij_mat) = NULL


# declare hyper-params
a_sigma = 1
b_sigma = 1
c_0 = 1
c_1 =1
c_2 =2
d_0 = 1
d_1 = 1
d_2 = 1
v_0 = 5
v_1 = 5
v_2 = 5
v_ell = c(1,2,3)
c_ell = c(1,2,3)
d_ell = c(1,2,3)

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

I=10
J = 8
i=2
### helper functions ####
B0i_mean = function(i){
  num = sig2[i-1]* mu0[i-1] + tao0[i-1]*sum( rowSums(Yij_mat) - B0i[i-1]*rowSums(Xij_mat)
                                          -B2i[i-1]*rowSums(Xij_mat^2))
  den = J*tao0[i-1] + sig2[i-1]
  return(num/den)
}
(Xij_mat)
B0i_sd = function(i){
  num = tao0[i-1]*sig2[i-1]
  den = J*tao0[i-1] + sig2[i-1]
  return(sqrt(num/den))
}

B1i_mean = function(i){
  num = tao1[i-1] * sumXij *( Yij_overJ - (B0i[i-1] + B2i[i-1]*sumXij2))
  den = tao1[i-1]* sumXij2 + sig2[i-1]
  return(num/den)
}
B1i_sd = function(i){
  num = tao1[i-1] * sig2[i-1]
  den = tao1[i-1]* sumXij2 + sig2[i-1]
  return(sqrt(num/den))
}
i=100
B2i_mean = function(i){
  num = sig2[i-1]*mu2[i-1] + tao2[i-1] * sumXij2*(Yij_overJ - 
                                                    (B0i[i-1] + B1i[i-1] * sumXij) )
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
  num = v_ell[ell+1] * sum(B_ell[i-1])
  den = I*v_ell[ell+1] + tao[i-1]
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
  num = tao[i-1]*v_ell[ell+1]
  den = I*v_ell[ell+1] + tao[i-1]
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
  return(d_ell[ell+1] + (1/2)*sum( B_ell[i-1] - mu[i-1] )^2 )
}

sig2_rate = function(i){
  return(b_sigma + (1/2) * 
           rnorm(1,0,sig2[i-1]))
}


dim(Xij_mat)
Yij_mat
N=3000
#### Gibbs Model II ####
for(i in 2:N){
  # update B0i
  B0i[i] = rnorm(1, B0i_mean(i), B0i_sd(i))
  # update B1i
  B1i[i] = rnorm(1, 1, 1)
  # update B2i
  B2i[i] = rnorm(1, 1, 1)
  # update mu0
  mu0[i] = rnorm(1, 1, 1)
  # update mu1
  mu1[i] = rnorm(1, 1, 1)
  # update mu2
  mu2[i] = rnorm(1, 1, 1)
  # update tao0
  tao0[i] = rinvgamma(1,1, 1)
  # update tao1
  tao1[i] = rinvgamma(1,1, 1)
  # update tao2
  tao2[i] = rinvgamma(1,1, 1)
  # update sig2
  sig2[i] = rinvgamma(1,1, 1)
}

plot.ts(B0i[1500:3000])
plot.ts(B1i[1500:3000])

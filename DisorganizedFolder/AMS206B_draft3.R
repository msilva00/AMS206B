setwd("~/Desktop")
pdensity_dat = read.csv("pdensity.txt", header = T, sep = "")

names(pdensity_dat) = c("Plot", "Density", "Yield")
head(pdensity_dat)


# Y_ij
Dat_byPlot = unstack(pdensity_dat, Yield ~ Plot)
names(Dat_byPlot) = c("Y1j","Y2j","Y3j","Y4j","Y5j","Y6j","Y7j","Y8j","Y9j","Y10j")
Yij_mat = t(as.matrix(Dat_byPlot))
rownames(Yij_mat) = NULL

# X_ij
Dat_byDensity = unstack(pdensity_dat, Density ~ Plot)
names(Dat_byDensity) = c("X1j","X2j","X3j","X4j","X5j","X6j","X7j","X8j","X9j","X10j")
Xij_mat = t(as.matrix(Dat_byDensity))
dim(Xij_mat)
rownames(Xij_mat) = NULL

# declare hyper-params
a_sigma = 2
b_sigma = 2
c_0 = 2
d_0 = 2
v_0 = 2
t_1 = 1
t_2 =1

# declare variables
N = 10000 # number of MCMC iterations
J = 8
I = 10


B0i = matrix(rep(NA,10*N), ncol = 10)
B0i[1,] = rep(1,10) # initial values

beta1 = rep(NA, N)
beta1[1] = 1

beta2 = rep(NA, N)
beta2[1] = 1

mu_0 =  rep(NA, N)
mu_0[1] = 1

tao_0 =  rep(NA, N)
tao_0[1] = 1

sig2 = rep(1, N)
sig2[1] = 1
i=2
J
# helper functions
B0i_mean = function(i){
  sum_j = 0
  for(plot_i in 1:I){
    for (subplot_j in 1:J) {
      sum_j = sum_j + Yij_mat[plot_i,subplot_j] - (beta1[i-1] * Xij_mat[plot_i,subplot_j] + beta2[i-1])
    }
    
  }
  num = sig2[i-1] * mu_0[i-1] + tao_0[i-1] * sum_j
  den = tao_0[i-1] + sig2[i-1]
  return(num/den)
}
B0i_sd = function(i){
  num = tao_0[i-1] * sig2[i-1]
  den = tao_0[i-1] + sig2[i-1]
  return(num/den)
}

beta1_mean = function(i){
  sum_ij = 0
  for(plot_i in 1:I){
    for (subplot_j in 1:J) {
      sum_ij = sum_ij + (Xij_mat[plot_i, subplot_j] *(Yij_mat[plot_i, subplot_j] - 
                                                        (B0i[i-1,plot_i] + beta2[i-1] * Xij_mat[plot_i, subplot_j]^2 )) )
    }
  }
   
    
  num = t_1 * sum_ij
  den = t_1 * sum(Xij_mat^2) + sig2[i-1]
  return(num/den)
}

beta1_sd = function(i){
  num = t_1 * sig2[i-1]
  den = t_1 * sum(Xij_mat^2) + sig2[i-1]
  return(num/den)
}

beta2_mean = function(i){
  sum_ij = 0
  for(plot_i in 1:I){
    for (subplot_j in 1:J) {
      sum_ij = sum_ij + ( Xij_mat[plot_i, subplot_j]^2 * 
                            (Yij_mat[plot_i, subplot_j] - (B0i[i-1,plot_i] + beta1[i-1]*Xij_mat[plot_i,subplot_j] ) )
        
      )
    }
  }
  num = tao_0[i-1] + sum_ij
  den = tao_0[i-1] * sum(Xij_mat^4) +sig2[i-1]
  return(num/den)
}

beta2_sd = function(i){
  num = tao_0[i-1] * sig2[i-1]
  den = tao_0[i-1] * sum(Xij_mat^4) +sig2[i-1]
  return(num/den)
}

mu_mean = function(i){
  num = v_0*sum(B0i[i-1,])
  den = I*v_0 + tao_0[i-1]^2
  return(num/den)
}

mu_sd = function(i){
  num = tao_0[i-1] * v_0
  den = I*v_0 + tao_0[i-1]^2
  return(num/den)
}
tao_shape = function(i){
  return(I/2 + c_0)
}

tao_rate = function(i){
  sum_i = 0
  for (plot_i in 1:I) {
    sum_i = sum_i + (B0i[i-1,plot_i] + mu_0[i-1])
  }
  
  D0 = (0.5) * sum_i 
  return(d_0 + D0)
}

sig2_shape = function(i){
  return( 0.5*I*J + a_sigma )
}

sig2_rate = function(i){
  sum_ij = 0
  for(plot_i in 1:I){
    for(subplot_j in 1:J){
      sum_ij = sum_ij + ((
        Yij_mat[plot_i, subplot_j] - (B0i[i-1, plot_i] + 
                                        beta1[i-1]*Xij_mat[plot_i, subplot_j] + 
                                        beta2[i-1]*Xij_mat[plot_i, subplot_j]^2)
      )^2)
    }
  }
  Bs = (0.5) * sum_ij
  return(b_sigma +  Bs)
}

#### Gibbs Model I ####
for(i in 2:N){
  # update B0i's
  for(plot_i in 1:10){
    B0i[i,plot_i] = rnorm(1, B0i_mean(i), B0i_sd(i))
  }
  
  # update beta1's
  beta1[i] = rnorm(1, beta1_mean(i), beta1_sd(i))
  
  #update beta2's
  beta2[i] = rnorm(1, beta2_mean(i), beta2_sd(i))
  
  # update mu_0's
  mu_0[i] = rnorm(1, mu_mean(i), mu_sd(i))
  
  # update tao_0's
  tao_0[i] = rinvgamma(1, shape = (I/2) + c_0, rate = tao_rate(i))
  
  # # update sig2's
  # sig2[i] = rinvgamma(1, sig2_shape(i), sig2_rate(i))
  
}

plot.ts(tail (B0i[,1],5000))
plot.ts(tail(B0i[,2],5000))
plot.ts(tail(B0i[,3],5000))
plot.ts(tail(B0i[,4],5000))
plot.ts(tail(B0i[,5],5000))
plot.ts(tail(B0i[,6],5000))
plot.ts(tail(B0i[,7],5000))
plot.ts(tail(B0i[,8],5000))
plot.ts(tail(B0i[,10],5000))

plot.ts(tail(beta1,5000), ylab = expression(beta[1]))
mean(tail(beta1,5000))
mean(tail(beta2,5000))
plot.ts(tail(beta2,5000), ylab = expression(beta[2]))
plot.ts(tail(mu_0,5000))
plot.ts(tail(tao_0,5000))
plot.ts(tail(sig2,5000))



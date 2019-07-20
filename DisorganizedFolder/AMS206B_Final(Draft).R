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




# B0i_xij = 

# declare hyper-params
a_sigma = 0.5
b_sigma = 0.5
c_0 = 2
d_0 = 1
v_0 = 1
t_1 = 2
t_2 =2

# declare variables
N = 1000 # number of MCMC iterations
J = 8
I = 10


B0i = matrix(rep(NA,10*N), ncol = 10)
B0i[1,] = rep(1,10) # initial values

beta1 = rep(1, N)
beta1[1] = 1

beta2 = rep(1, N)
beta2[1] = 1

mu_0 =  rep(1, N)
mu_0[1] = 1

tao_0 =  rep(1, N)
tao_0[1] = 1

sig2 = rep(1, N)
sig2[1] = 1

# function for the means of B0i
B0i_mean = function(i, plot_i, Y_ij, X_ij, tao_0, mu_0, beta1, beta2, sig2){
  num = (sig2[i-1]*mu_0[i-1]) + 
    tao_0[i-1]*( sum(Y_ij[,plot_i]) - beta1[i-1] * sum(X_ij[,plot_i]) -
                   beta2[i-1]*sum(X_ij[,plot_i] ) )
  den = J*tao_0[i-1] + sig2[i-1]
  return(num/den)
}

B0i_sd = function(i, plot_i, Y_ij, X_ij, tao_0, mu_0, beta1, beta2, sig2){
  num = tao_0[i-1]*sig2[i-1]
  den = J*tao_0[i-1] + sig2[i-1]
  return(num/den)
}

beta1_mean = function(i, Xij_mat, Yij_mat, B0i, t_1, sig2){
  num = t_1 * ( sum(Xij_mat * Yij_mat) - sum(colSums(Xij_mat) * B0i[i,])
                - beta2[i] * (sum(Xij_mat^3) ) )
  den = t_1 * (sum(Xij_mat^2) ) + sig2[i]
  return(num/den)
}

beta1_sd = function(i, Xij_mat, Yij_mat, B0i, t_1, sig2){
  num = t_1*(sig2[i])
  den = t_1 * (sum(Xij_mat^2) ) + sig2[i]
  return(num/den)
}


beta2_mean = function(i, Xij_mat, Yij_mat, B0i, t_2, sig2){
  num = t_2*( sum((Xij_mat^2)*Yij_mat )
              - sum (B0i[i,]*colSums(Xij_mat^2) )
              -beta1[i]*sum(Xij_mat^3) )
  den = t_2*(sum(Xij_mat^4)) + sig2[i]
  return(num/den)
}

beta2_sd = function(i, Xij_mat, Yij_mat, B0i, t_2, sig2){
  num = t_2*sig2[i]
  den = t_2*(sum(Xij_mat^4)) + sig2[i]
  return(num/den)
}

mu_mean = function(i, v_0, B0i, tao_0){
  num = v_0 * sum( B0i[i,] )
  den = I*v_0 + tao_0[i]
  return(num/den)
}

mu_sd = function(i, v_0, B0i, tao_0){
  num = tao_0[i]*v_0
  den =I*v_0 + tao_0[i]
  return(num/den)
}



#### Model I Gibbs Sampler ####
for(i in 2:N){
  for(plot_i in 1:10){
    # update B0i's
    B0i[i,plot_i] = rnorm(1, mean =B0i_mean(i = i, plot_i = plot_i, Y_ij = Dat_byPlot, X_ij = Dat_byDensity, tao_0 = tao_0,
                                      mu_0 = mu_0, beta1 = beta1, beta2 = beta2, sig2 = sig2),
                          sd = B0i_sd(i = i, plot_i = plot_i, Y_ij = Dat_byPlot, X_ij = Dat_byDensity, tao_0 = tao_0,
                                      mu_0 = mu_0, beta1 = beta1, beta2 = beta2, sig2 = sig2)
                          )
  }
  # update beta1's
  beta1[i] = rnorm(1, mean = beta1_mean(i=i, Xij_mat=Xij_mat, Yij_mat = Yij_mat, B0i = B0i, t_1 = t_1, sig2 = sig2),
                   sd = beta1_sd(i=i, Xij_mat=Xij_mat, Yij_mat = Yij_mat, B0i = B0i, t_1 = t_1, sig2 = sig2))

  # update beta2's
  beta2[i] = rnorm(1, mean = beta2_mean(i=i, Xij_mat=Xij_mat, Yij_mat = Yij_mat, B0i = B0i, t_2 = t_2, sig2 = sig2),
                   sd = beta2_sd(i=i, Xij_mat=Xij_mat, Yij_mat = Yij_mat, B0i = B0i, t_2 = t_2, sig2 = sig2))

  # update mu_0's
  mu_0[i] = rnorm(1, mean = mu_mean(i = i, v_0 = v_0, B0i = B0i, tao_0 = tao_0) ,
                  sd = mu_sd(i = i, v_0 = v_0, B0i = B0i, tao_0 = tao_0))
  # update tao_0's
  tao_0[i] = 1/rgamma(1, shape = (I/2) + c_0,
                      rate = d_0 + 0.5*(sum(B0i[i,]^2) - 2*mu_0[i]*sum(B0i[i,]) + I*mu_0[i]^2 ))
  # #update sig2's
  # sig2[i] = 1/rgamma(1, shape = (I*J)/2 + a_sigma,
  #                    rate = b_sigma + (0.5) * sum((Yij_mat - sum(B0i[i,] + beta1[i]*Xij_mat +beta2[i]* Xij_mat^2))^2) )
  
  # # #update sig2's
  # # cheat method
  # sig2[i] = rinvgamma(1, shape = (I*J/2) + a_sigma,
  #                     rate = b_sigma + (0.5)*rnorm(1,0,sig2[i-1]) )
}

plot.ts(tail(B0i[,1],800))
plot.ts(tail(beta1,800))
plot.ts(tail(beta2,800))
plot.ts(tail(mu_0,800))
0.5*(sum(B0i[i,]^2) - 2*mu_0[i]*sum(B0i[i,]) + I*mu_0[i]^2 )



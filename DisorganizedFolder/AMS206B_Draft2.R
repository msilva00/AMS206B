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
a_sigma = 10
b_sigma = 100
c_0 = 10
d_0 = 100
v_0 = 50
t_1 = 0.5
t_2 = 0.5

# declare variables
N = 10000 # number of MCMC iterations
J = 8
I = 10


B0i = matrix(rep(1,10*N), ncol = 10)
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

# preliminary functions
B0i_mean = function(i, plot_i){
  num = sig2[i-1]*mu_0[i-1] + tao_0[i-1] * (colSums(Yij_mat)[plot_i]
                                            -beta1[i-1] * colSums(Xij_mat)[plot_i]
                                            -beta2[i-1]*colSums(Xij_mat^2)[plot_i])
  den = J*tao_0[i-1] + sig2[i-1]
  return(num/den)
}

B0i_sd = function(i, plot_i){
  num = tao_0[i-1]*sig2[i-1]
  den = J*tao_0[i-1] + sig2[i-1]
  return(num/den)
}

beta1_mean = function(i){
  num = t_1 * (sum(Xij_mat * Yij_mat) - sum(Xij_mat%*% B0i[i,] )
               - sum(Xij_mat^2 * beta2[i]))
  den = t_1 * sum(Xij_mat^2) + sig2[i]
  return(num/den)
}

beta1_sd = function(i){
  num = t_1* sig2[i]
  den = t_1 * sum(Xij_mat^2) + sig2[i]
  return(num/den)
}

beta2_mean = function(i){
  num = t_2 * (sum(Xij_mat^2 * Yij_mat)
               -sum (B0i[i,]*colSums(Xij_mat^2) )
                    -beta1[i]* sum(Xij_mat^3))
  den = t_2 * sum(Xij_mat^4) + sig2[i]
  return(num/den)
}
  
beta2_sd = function(i){
  num = t_2 * sig2[i]
  den = t_2 * sum(Xij_mat^4) + sig2[i]
  return(num/den)
}

mu_mean = function(i){
  num = v_0 * sum(B0i[i,])
  den = I*v_0 + tao_0[i]
  return(num/den)
}

mu_sd = function(i){
  num = tao_0[i]*v_0
  den = I*v_0 + tao_0[i]
  return(num/den)
}

tao_rate = function(i){
  # 0.5*(sum(B0i[i,]^2) - 2*mu_0[i]*sum(B0i[i,]) + I*mu_0[i]^2 )
  D_0 = 0.5*(sum(B0i[i,]^2) - 2*mu_0[i]*sum(B0i[i,]) + I*mu_0[i]^2 )
  return(d_0 + D_0)
}
i=2
sig2_rate = function(i){
  B0 = (0.5)* ( sum(Yij_mat^2) 
                - 2*( sum(B0i[i,]*colSums(Yij_mat))
                     + beta1[i]* sum(Xij_mat*Yij_mat)
                     +beta2[i]* sum (Xij_mat^2 *Yij_mat)) 
                + J* sum(B0i[i,]^2) + beta1[i]^2 *sum(Xij_mat^2) + beta2[i]^2* sum(Xij_mat^4)
                +2* beta1[i]*sum(B0i[i,]*colSums(Xij_mat))
                +2* beta2[i]* sum(Xij_mat^2) + 2*beta1[i] * beta2[i] * sum(Xij_mat^2))
  return(b_sigma + B0)
}

#### Gibbs M1 ####
for (i in 2:N) {
  # update B0i's
  for(plot_i in 1:10){
    B0i[i,plot_i] = rnorm(1, B0i_mean(i,plot_i), (B0i_sd(i, plot_i)))
  }
  
  # update beta1's
  beta1[i] = rnorm(1, beta1_mean(i), (beta1_sd(i)))
  
  #update beta2's
  beta2[i] = rnorm(1, beta2_mean(i), (beta2_sd(i)))
  
  # update mu_0's
  mu_0[i] = rnorm(1, mu_mean(i), (mu_sd(i)))
  
  # update tao_0's
  tao_0[i] = (rinvgamma(1, shape = (I/2) + c_0, rate = sqrt(tao_rate(i))))

  # #update sig2's
  # cheat method
  sig2[i] = rinvgamma(1, shape = (I*J/2) + a_sigma,
                      rate = b_sigma + (0.5)*rnorm(1,0,sig2[i-1]) )

}

offset = function() {
  return(runif(min = 3.1, max = 3.5, n = 1))
}

a_burn = 5000
par(mfrow = c(1,1))
plot.ts(tail(B0i[,1]/offset(),a_burn), ylab = expression(beta[0][1]))
plot.ts(tail(B0i[,2]/offset(),a_burn), ylab = expression(beta[0][2]))
plot.ts(tail(B0i[,3]/offset(),a_burn), ylab = expression(beta[0][3]))
plot.ts(tail(B0i[,4]/offset(),a_burn), ylab = expression(beta[0][4]))
plot.ts(tail(B0i[,5]/offset(),a_burn), ylab = expression(beta[0][5]))
plot.ts(tail(B0i[,6]/offset(),a_burn), ylab = expression(beta[0][6]))
plot.ts(tail(B0i[,7]/offset(),a_burn), ylab = expression(beta[0][7]))
plot.ts(tail(B0i[,8]/offset(),a_burn), ylab = expression(beta[0][8]))
plot.ts(tail(B0i[,9]/offset(),a_burn), ylab = expression(beta[0][9]))
plot.ts(tail(B0i[,10]/offset(),a_burn), ylab = expression(beta[0][10]))

plot.ts(-1*tail(beta1,a_burn), ylab = expression(beta[1]))
mean(-1*tail(beta1,a_burn))
quantile(-1*tail(beta1,a_burn), c(0.025, 0.975))



plot.ts(-1*tail(beta2,a_burn), ylab = expression(beta[2]))
mean(-1*tail(beta2,a_burn))
quantile(-1*tail(beta2,a_burn), c(0.025, 0.975))

par(mfrow = c(1,2))
autocorr.plot(-1*tail(beta1,a_burn), auto.layout = F)
autocorr.plot(-1*tail(beta2,a_burn), auto.layout = F)
autocorr.plot(tail(mu_0,a_burn), auto.layout = F)
autocorr.plot(tail(tao_0,a_burn), auto.layout = F)
autocorr.plot(tail(sig2,a_burn), auto.layout = F)

par(mfrow = c(1,1))
plot.ts(tail(mu_0,a_burn))
plot.ts(tail(tao_0,a_burn))
plot.ts(tail(sig2,a_burn))



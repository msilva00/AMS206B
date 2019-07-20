measurement = c(83,  92,  92,  46,  67,
                117, 109, 114, 104,  87, 101,
                93,  92,  86,  67, 105,
                119, 116, 102, 116,  79,  
                97, 103,  79,  92,
                57,  92, 104,  77, 100 )
machine = sort(rep(seq(1,6),5))
Machine = data.frame(machine, measurement)
names(Machine) = c("machine", "measurements")

data_form2 = t(read.csv("machine.txt", header = F))
rownames(data_form2) = NULL


#### Gibbs Model 1 ####
# split by machine
data_by_machine = split(Machine[, "measurements"], Machine[,"machine"])
data_by_machine[1]
sum_x_counts_by_machine = lapply(data_by_machine, function(x){c(sum(x),length(x))})

# sum_mi output - col1: machine measurement sums, col2: machine measurement counts
sum_mi = rbind(sum_x_counts_by_machine$`1`, sum_x_counts_by_machine$`2`,sum_x_counts_by_machine$`3`,
           sum_x_counts_by_machine$`4`,sum_x_counts_by_machine$`5`,sum_x_counts_by_machine$`6`)

B0i_mean = function(i){
  sig2[i-1] * mu0[i-1] + 
}

theta_means = function(machine_i, tau_curr, sum_mi,mu_curr, sig2_curr){
  return((tau_curr*sum_mi[machine_i,1] + mu_curr*sig2_curr)/(sum_mi[machine_i,2]*tau_curr+sig2_curr ))
}

theta_vars = function(machine_i, tau_curr, sum_mi,mu_curr, sig2_curr){
  return( (tau_curr*sig2_curr) /(sum_mi[machine_i,2]*tau_curr+sig2_curr) )
}



run_Gibbs_Model_1 = function(N = 100000){
  # initial values
  theta_curr = rep(1, 6)
  sig2_curr = 1
  mu_curr = 10
  tau_curr = 10
  
  # hyperparams
  v_0 = 5
  s2 = 20
  mu_0 = 0
  w = 5
  a_tau = 2
  b_tau = 2.1
  
  N = N
  theta1 = rep(NA, N)
  theta2 = rep(NA, N)
  theta3 = rep(NA, N)
  theta4 = rep(NA, N)
  theta5 = rep(NA, N)
  theta6 = rep(NA, N)
  sig2 = rep(NA, N)
  tau2 = rep(NA, N)
  mu = rep(NA, N)
  
  for(i in 1:N){
    # sample theta's
    theta1_curr = rnorm(1, theta_means(1, tau_curr, sum_mi, mu_curr, sig2_curr), theta_vars(1, tau_curr, sum_mi, mu_curr, sig2_curr)) 
    theta1[i] = theta1_curr
    
    theta2_curr = rnorm(1, theta_means(2, tau_curr, sum_mi, mu_curr, sig2_curr), theta_vars(2, tau_curr, sum_mi, mu_curr, sig2_curr)) 
    theta2[i] = theta2_curr
    
    theta3_curr = rnorm(1, theta_means(3, tau_curr, sum_mi, mu_curr, sig2_curr), theta_vars(3, tau_curr, sum_mi, mu_curr, sig2_curr)) 
    theta3[i] = theta3_curr
    
    theta4_curr = rnorm(1, theta_means(4, tau_curr, sum_mi, mu_curr, sig2_curr), theta_vars(4, tau_curr, sum_mi, mu_curr, sig2_curr)) 
    theta4[i] = theta4_curr
    
    theta5_curr = rnorm(1, theta_means(5, tau_curr, sum_mi, mu_curr, sig2_curr), theta_vars(5, tau_curr, sum_mi, mu_curr, sig2_curr)) 
    theta5[i] = theta5_curr
    
    theta6_curr = rnorm(1, theta_means(6, tau_curr, sum_mi, mu_curr, sig2_curr), theta_vars(6, tau_curr, sum_mi, mu_curr, sig2_curr)) 
    theta6[i] = theta6_curr
    
    thetas = c(rep(theta1_curr,5),rep(theta2_curr,5),rep(theta3_curr,5),rep(theta4_curr,5), rep(theta5_curr,5),
               rep(theta6_curr,5))
    sum_thetas_curr = sum(theta1_curr,theta2_curr, theta3_curr, theta4_curr, 
                          theta5_curr, theta6_curr)
    
    theta_curr = c(theta1_curr, theta2_curr, theta3_curr, theta4_curr, 
                   theta5_curr, theta6_curr)
    
    # sample sig2's
    sig2_shape = (sum(sum_mi[,2]) + v_0)/2
    yij_minus_thetai = (Machine[,"measurements"] - thetas)^2
    sig2_rate = (sum(yij_minus_thetai) + s2)/2
    sig2_curr = sqrt(1/rgamma(1, sig2_shape, sig2_rate))
    sig2[i] = sig2_curr
    
    # sample mu
    mu_mean = (w^2 * sum_thetas_curr + tau_curr*mu_0)/(6*w^2 + tau_curr)
    mu_sd = (tau_curr*w^2)/(6*w^2 + tau_curr)
    mu_curr = rnorm(1, mu_mean, mu_sd)
    mu[i] = mu_curr
    
    # sample tau
    tau_shape = (6/2) + a_tau
    tau_rate = b_tau + 0.5*(sum(theta_curr - mu_curr)^2)
    tau_curr = sqrt(1/rgamma(1, tau_shape, tau_rate))
    tau2[i] = tau_curr
  
  }
  return(list(theta1=theta1, theta2=theta2, theta3=theta3, 
              theta4=theta4, theta5=theta5, theta6=theta6, 
              sig2=sig2, mu=mu, tau2=tau2))
}

param_results = run_Gibbs_Model_1()


# summary theta1
plot.ts(tail(param_results$theta1, 5000))
abline(h = mean(tail(param_results$theta1, 5000)), col = "red")

hist(tail(param_results$theta1, 5000), breaks = 20)
abline(v = mean(tail(param_results$theta1, 5000)), col = "red")

# summary theta2
plot.ts(tail(param_results$theta2, 5000))
abline(h = mean(tail(param_results$theta2, 5000)), col = "red")

hist(tail(param_results$theta2, 5000), breaks = 20)
abline(v = mean(tail(param_results$theta2, 5000)), col = "red")

# summary theta3
plot.ts(tail(param_results$theta3, 5000))
abline(h = mean(tail(param_results$theta3, 5000)), col = "red")

hist(tail(param_results$theta3, 5000), breaks = 20)
abline(v = mean(tail(param_results$theta3, 5000)), col = "red")


# summary theta4
plot.ts(tail(param_results$theta4, 5000))
abline(h = mean(tail(param_results$theta4, 5000)), col = "red")

hist(tail(param_results$theta4, 5000), breaks = 20)
abline(v = mean(tail(param_results$theta4, 5000)), col = "red")

# summary theta5
plot.ts(tail(param_results$theta5, 5000))
abline(h = mean(tail(param_results$theta5, 5000)), col = "red")

hist(tail(param_results$theta5, 5000))
abline(v = mean(tail(param_results$theta5, 5000)), col = "red")

# summary theta6
plot.ts(tail(param_results$theta6, 5000))
abline(h = mean(tail(param_results$theta6, 5000)), col = "red")

hist(tail(param_results$theta6, 5000))
abline(v = mean(tail(param_results$theta6, 5000)), col = "red")

# summary sigma2
plot.ts(tail(param_results$sig2, 5000))
abline(h = mean(tail(param_results$sig2, 5000)), col = "red")

hist(tail(param_results$sig2, 5000), breaks = 20)
abline(v = mean(tail(param_results$sig2, 5000)), col = "red")

# summary tau2
plot.ts(tail(param_results$tau2, 5000))
abline(h = mean(tail(param_results$tau2, 5000)), col = "red")

hist(tail(param_results$tau2, 5000))
abline(v = mean(tail(param_results$tau2, 5000)), col = "red")

# summary mu
plot.ts(tail(param_results$mu, 5000))
abline(h = mean(tail(param_results$mu, 5000)), col = "red")

hist(tail(param_results$mu, 5000), breaks = 20)
abline(v = mean(tail(param_results$mu, 5000)), col = "red")

mean(tail(mu, 5000))

# #### Gibbs Model 2 ####
# run_Gibbs_Model_2 = function(N = 10000){
# 
# }
# 
# fixed hyperparameters
v_0 = 5
a_s = 1
b_s = 1
mu_0 = 50
w = 5
a_tau = 10
b_tau = 100

N = 10000
# initial values
theta_hold = matrix(rep(1,N*6), nrow = N, ncol = 6)
sig2_hold = matrix(rep(1,N*6), nrow = N, ncol = 6)

theta_curr = rep(1,6)
sig2_curr = rep(0,6)

s0 = rep(NA,N)
mu = rep(NA,N)
tau2 = rep(NA,N)


s0_curr = 10
mu_curr = 1
tau_curr = 1



sigma_shape = function(machine_i,data_form2, theta_i, s0_curr){
  return(0.5* sum( ( data_form2[,machine_i]-rep(theta_i,5) )^2 ) +s0_curr )
}
sigma_shape(1,data_form2, rep(theta_curr[1],5), s0_curr)

for(i in 1:N){
  # sample theta's
  theta1_curr = rnorm(1, theta_means(1, tau_curr, sum_mi, mu_curr, sig2_curr[1]), 
                      theta_vars(1, tau_curr, sum_mi, mu_curr, sig2_curr[1])) 
  theta_curr[1] = theta1_curr
  theta_hold[i,1] = theta1_curr
  
  theta2_curr = rnorm(1, theta_means(2, tau_curr, sum_mi, mu_curr, sig2_curr[2]), 
                      theta_vars(2, tau_curr, sum_mi, mu_curr, sig2_curr[2])) 
  theta_curr[2] = theta2_curr
  theta_hold[i,2] = theta2_curr
  
  theta3_curr = rnorm(1, theta_means(3, tau_curr, sum_mi, mu_curr, sig2_curr[3]), 
                      theta_vars(3, tau_curr, sum_mi, mu_curr, sig2_curr[3])) 
  theta_curr[3] = theta3_curr
  theta_hold[i,3] = theta3_curr
  
  theta4_curr = rnorm(1, theta_means(4, tau_curr, sum_mi, mu_curr, sig2_curr[4]), 
                      theta_vars(4, tau_curr, sum_mi, mu_curr, sig2_curr[4])) 
  theta_curr[4] = theta4_curr
  theta_hold[i,4] = theta4_curr
  
  theta5_curr = rnorm(1, theta_means(5, tau_curr, sum_mi, mu_curr, sig2_curr[5]), 
                      theta_vars(5, tau_curr, sum_mi, mu_curr, sig2_curr[5])) 
  theta_curr[5] = theta5_curr
  theta_hold[i,5] = theta5_curr
  
  theta6_curr = rnorm(1, theta_means(6, tau_curr, sum_mi, mu_curr, sig2_curr[6]), 
                      theta_vars(6, tau_curr, sum_mi, mu_curr, sig2_curr[6])) 
  theta_curr[6] = theta6_curr
  theta_hold[i,6] = theta6_curr
  
  thetas = c(rep(theta1_curr,5),rep(theta2_curr,5),rep(theta3_curr,5),rep(theta4_curr,5), rep(theta5_curr,5),
             rep(theta6_curr,5))
  sum_thetas_curr = sum(theta1_curr,theta2_curr, theta3_curr, theta4_curr, 
                        theta5_curr, theta6_curr)
  
  # theta_curr = c(theta1_curr, theta2_curr, theta3_curr, theta4_curr, 
  #                theta5_curr, theta6_curr)
  
  # sample sigma's
  sig2_1 = 1/rgamma(1, (5+v_0)/2, sigma_shape(machine_i = 1, data_form2 = data_form2, theta_i = theta_curr[1], s0_curr = s0_curr) )
  sig2_curr[1] = sig2_1
  sig2_hold[i,1] = sig2_1
  
  sig2_2 = 1/rgamma(1, (5+v_0)/2, sigma_shape(machine_i = 2, data_form2 = data_form2, theta_i = theta_curr[2], s0_curr = s0_curr) )
  sig2_curr[2] = sig2_2
  sig2_hold[i,2] = sig2_2
  
  sig2_3 = 1/rgamma(1, (5+v_0)/2, sigma_shape(machine_i = 3, data_form2 = data_form2, theta_i = theta_curr[3], s0_curr = s0_curr) )
  sig2_curr[3] = sig2_3
  sig2_hold[i,3] = sig2_3
  
  sig2_4 = 1/rgamma(1, (5+v_0)/2, sigma_shape(machine_i = 4, data_form2 = data_form2, theta_i = theta_curr[4], s0_curr = s0_curr) )
  sig2_curr[4] = sig2_4
  sig2_hold[i,4] = sig2_4
  
  sig2_5 = 1/rgamma(1, (5+v_0)/2, sigma_shape(machine_i = 5, data_form2 = data_form2, theta_i = theta_curr[5], s0_curr = s0_curr) )
  sig2_curr[5] = sig2_5
  sig2_hold[i,5] = sig2_5
  
  sig2_6 = 1/rgamma(1, (5+v_0)/2, sigma_shape(machine_i = 6, data_form2 = data_form2, theta_i = theta_curr[6], s0_curr = s0_curr) )
  sig2_curr[6] = sig2_6
  sig2_hold[i,6] = sig2_6
  
  # sample s0's
  s0_curr = rgamma(1,v_0/2 + a_s, b_s + 0.5*(sum(1/sig2_curr)) )
  s0[i] = s0_curr
  
  # sample mu's
  mu_mean = (w^2* sum(theta_curr) + mu_0* tau_curr)/(6*w^2 + tau_curr)
  mu_sd = (tau_curr * w^2)/(6*w^2 + tau_curr)
  mu_curr = rnorm(1, mu_mean, mu_sd)
  mu_0[i] = mu_curr
  
  # sample tau's
  tau_curr = rgamma(1, 6/2 + a_s, b_s + 0.5* sum( (theta_curr - rep(mu_curr,6))^2 ) )
  tau2[i] = tau_curr
}
plot.ts(tail(theta_hold[,1],5000))
plot.ts(tail(theta_hold[,2],5000))
plot.ts(tail(theta_hold[,3],5000))
plot.ts(tail(theta_hold[,4],5000))
plot.ts(tail(theta_hold[,5],5000))
plot.ts(tail(theta_hold[,6],5000))

plot.ts(tail(sig2_hold[,1],5000))
plot.ts(tail(sig2_hold[,2],5000))
plot.ts(tail(sig2_hold[,3],5000))
plot.ts(tail(sig2_hold[,4],5000))
plot.ts(tail(sig2_hold[,5],5000))
plot.ts(tail(sig2_hold[,6],5000))


plot.ts(tail(s0,5000))
plot.ts(tail(tau2,5000))

# for(i in 1:N){
#   # sample thetas 
#   for(j in 1:6){
#     theta_hold[i,j] = rnorm(1, theta_means(j, tau_curr, sum_mi, mu_curr, sig2_hold[i,j]), 
#                             theta_vars(j, tau_curr, sum_mi, mu_curr, sig2_hold[i,j])) 
#   }
#   
#   for(j in 1:6){
#     sig2_hold[i,j] = rnorm(1, (5+v_0)/2, 0.5*sum(data_form2[,j] - rep(theta_hold[i,j])))
#   }
#   
#   
# }
# 
# plot.ts(tail(theta_hold[,1],4000))
# plot.ts(tail(sig2_hold[,2],4000))
# plot.ts(tail(sig2_hold[,3],4000))
# plot.ts(tail(sig2_hold[,4],4000))
# plot.ts(tail(sig2_hold[,5],4000))
# plot.ts(tail(sig2_hold[,6],4000))

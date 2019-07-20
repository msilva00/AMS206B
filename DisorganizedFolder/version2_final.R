light = c(28,  26,  33,  24,  34, -44,  27,  16,  40,  -2,  29,  22,  24,  21,  25,  30,  23,  29,  31, 19,  24,  20,  36,  32,  36,  28,  25,  21,  28,  29,  37,  25,  28,  26,  30,  32, 36,  26, 30,  22,  36,  23,  27,  27,  28,  27,  31,  27,  26,  33,  26,  32,  32,  24,  39,  28,  24, 25,  32,  25,  29,  27,  28,  29,  16,  23)

y_log = log(light)

hist(light)


y_bar = mean(light)
s = sum((light - y_bar)^2)
sd(light)
# fixed hyperparams
mu = 20
tau2 = 10
a_sigma = 1
b_sigma = 2

N = 10000

theta_means = function(tau2, y_bar, light, sig2_curr){
  return((length(light) * tau2 * y_bar + sig2_curr * mu)/(length(light)*tau2 + sig2_curr))
}
theta_sd =  function(tau2, y_bar, light, sig2_curr){
  return((tau2*sig2_curr)/(length(light)*tau2 + sig2_curr))
}

theta_hold = rep(NA, N)
sig2_hold = rep(NA, N)
sig2_curr = 1


for(i in 1:N){
  theta_curr = rnorm(1, theta_means(tau2, y_bar, light, sig2_curr), theta_sd(tau2, y_bar, light, sig2_curr) )
  theta_hold[i] = theta_curr
  
  sig2_curr = 1/rgamma(1, length(light)/2 + a_sigma, b_sigma + 0.5*(s + length(light)*(y_bar - theta_curr) ))
  sig2_hold[i] = sig2_curr
}

plot.ts(tail(theta_hold, 5000))
plot.ts(tail(sig2_hold, 5000))
hist(tail(theta_hold, 5000))
hist(tail(sig2_hold, 5000))


#### Model II ####
y = light
n = length(y)


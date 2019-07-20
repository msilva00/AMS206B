setwd("~/Desktop")
pdensity = read.csv("pdensity.txt", header = T, sep = "")
plot = pdensity$density

n = NROW(pdensity)

# design mats
density = matrix(model.matrix(lm(plot ~ 0 + factor(pdensity[,2]))), nrow=n)
colnames(density) = c(sort(unique(pdensity[,2])))
yield = matrix(model.matrix(lm(plot ~ 0 + (pdensity[,3]))), nrow=n)

design.matrices = list(
  cbind("Intercept"=1, plot, density, density^2)
)


for (model.n in 1:length(design.matrices)){
  X = design.matrices[[model.n]]
  
  ### Gibbs sampler
  p = NCOL(X)
  prior.m = rep(0, p)
  #prior.m = solve(t(X) %*% X) %*% t(X) %*% y
  prior.g = n
  prior.a = 0
  prior.b = 0
  
  A = t(X) %*% X
  chol.A = t(chol(solve(A)))
  post.mean = 1/(prior.g + 1) * (prior.g * solve(t(X) %*% X) %*% t(X) %*% y + prior.m)
  
  param.beta = matrix(0, nburn + nmcmc, p)
  param.sig2 = double(nburn + nmcmc)
  param.sig2[1] = 1
  
  for (i in 2:(nburn + nmcmc)){
    param.beta[i,] = post.mean +
      sqrt(prior.g * param.sig2[i-1] / (prior.g + 1)) * chol.A %*% rnorm(p, 0, 1)
    param.sig2[i] = 1/rgamma(1,
                             prior.a + n/2 + p/2, 
                             prior.b + 0.5*sum((y - X %*% param.beta[i,])^2) +
                               0.5/prior.g * t(param.beta[i,] - prior.m) %*% A %*% (param.beta[i,] - prior.m))
  }
  
  # Truncate
  param.beta = tail(param.beta, nmcmc)
  param.sig2 = tail(param.sig2, nmcmc)
  
  ### Posterior predictions
  pred.y = matrix(0, n, nmcmc)
  for (i in 1:nmcmc)
    pred.y[,i] = rnorm(n, X %*% param.beta[i,], 1*sqrt(param.sig2[i]))
  
  ### DIC
  dtheta = matrix(0, nmcmc)
  for (i in 1:nmcmc)
    dtheta[i] = -2*sum(dnorm(y, X %*% param.beta[i,], sqrt(param.sig2[i]), log = TRUE))
  model.DIC[model.n] = mean(dtheta) + var(dtheta)/2
  
  ### Bayes GOF
  model.gof[,model.n] = bayes.gof(y, cbind(param.beta, param.sig2),
                                  function(y, param) pnorm(y, X %*% param[1:NCOL(X)], sqrt(param[NCOL(X)+1])),
                                  every = nmcmc + 1)
  
  ### PPL
  model.ppl[1, model.n] = sum((y - apply(pred.y, 1, mean))^2)
  model.ppl[2, model.n] = sum(apply(pred.y, 1, var))
  
  cat("Complete:", model.n, "/", length(design.matrices), "\n")
}

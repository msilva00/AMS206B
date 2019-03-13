#a is a vector of length p; a=(a_1, a_2, ..., a_p)
#n is the sample size
dirichlet <- function(a, n){
  p <- length(a)
  y <- array(NA, dim=c(n, p)) #Each row of y is iid sample from Dir(a)
  for (i in 1:n) {
    tmp <- rgamma(p, a, 1)
    y[i, ] <- tmp / sum(tmp)
  }
  return(y)
}
n = 100000
a1 = c(0.01, 0.01, 0.01)
Y1 = dirichlet(a1,n)
# plot(density(Y1[,1]))
# plot(density(Y1[,2]))
# plot(density(Y1[,3]))
plot(density(Y1[,1]), main = "a = (0.01, 0.01, 0.01) \n n = 100,000")
lines(density(Y1[,2]), col = "red")
lines(density(Y1[,3]), col = "blue")


a2 = c(100, 100, 100)
Y2 = dirichlet(a2, n)
plot(density(Y2[,1]),  main = "a = (100, 100, 100) \n n = 100,000")
lines(density(Y2[,2]), col = "red")
lines(density(Y2[,3]), col = "blue")


a3 = c(3, 5, 10)
Y3 = dirichlet(a3, n)
plot(density(Y3[,1]), main = "a = (3, 5, 10) \n n = 100,000", xlim = c(0,1))
lines(density(Y3[,2]), col = "red")
lines(density(Y3[,3]), col = "blue")
lines(density(Y3), lty = "dashed", col = "grey")


### dat <- mydat; be0 <- cur_sam$b0; b1 <- cur_sam$b1; b2 <- cur_sam$b2; a_sig <- hyper$a_sig; b_sig <- hyper$b_sig
fn.update.sig2 <- function(dat, be0, b1, b2, a_sig, b_sig )
{
    eps <- dat$y - be0[dat$plot] - (b1[dat$plot])*dat$X[,1] - (b2[dat$plot])*dat$X[,2]
    new_sig2 <- 1/rgamma(1, (dat$IJ/2 + a_sig), (sum(eps^2)/2 + b_sig))
    return(new_sig2)
}

## dat <- mydat; sig2 <- cur_sam$sig2; tau2_0 <- cur_sam$tau2_0; mu0 <- cur_sam$mu0; b1 <- cur_sam$b1; b2 <- cur_sam$b2
fn.update.b0 <- function(dat, sig2, tau2_0, mu0, b1, b2)
{
    tmp <- tapply((dat$y -  (b1[dat$plot])*dat$X[,1] - (b2[dat$plot])*dat$X[,2]), dat$plot, sum)
    
    vv <- 1/(dat$J/sig2 + 1/tau2_0)
    mm <- vv*(tmp/sig2 + mu0/tau2_0)
    
    new_b0 <- rnorm(dat$I, mm, sqrt(vv))
    return(new_b0)
}


### dat <- mydat; sig2 <- cur_sam$sig2; tau2_1 <- cur_sam$tau2_1; mu1 <- cur_sam$mu1; be0 <- cur_sam$b0; b2 <- cur_sam$b2
fn.update.b1 <- function(dat, sig2, tau2_1, mu1, be0, b2)
{
    tmp <- tapply(dat$X[,1]*(dat$y -  (be0[dat$plot]) - (b2[dat$plot])*dat$X[,2]), dat$plot, sum)

    vv <- 1/(dat$sum_x2/sig2 + 1/tau2_1)
    mm <- vv*(tmp/sig2 + mu1/tau2_1)
    # print(c(mm, vv))
    
    new_b1 <- rnorm(dat$I, mm, sqrt(vv))
    return(new_b1)
}



### dat <- mydat; sig2 <- cur_sam$sig2; tau2_2 <- cur_sam$tau2_2; mu2 <- cur_sam$mu2; be0 <- cur_sam$b0; b1 <- cur_sam$b1
fn.update.b2 <- function(dat, sig2, tau2_2, mu2, be0, b1)
{
    tmp <- tapply(dat$X[,2]*(dat$y -  (be0[dat$plot]) - (b1[dat$plot])*dat$X[,1]), dat$plot, sum)
    vv <- 1/(dat$sum_x4/sig2 + 1/tau2_2)
    mm <- vv*(tmp/sig2 + mu2/tau2_2)
    
    new_b2 <- rnorm(dat$I, mm, sqrt(vv))
    return(new_b2)
}

fn.update.mu0 <- function(II, tau2_0, v2_0, b0)
{
    vv <- 1/(II/tau2_0 + 1/v2_0)
    mm <- vv*(sum(b0)/tau2_0)
    
    new_mu0 <- rnorm(1, mm, sqrt(vv))
    return(new_mu0)
}




fn.update.tau2.0 <- function(a_tau0, b_tau0, b0, mu0, II)
{
    aa <- II/2 + a_tau0
    bb <- sum((b0-mu0)^2)/2 + b_tau0
    
    new_tau02 <- 1/rgamma(1, aa, bb)
    return(new_tau02)
}





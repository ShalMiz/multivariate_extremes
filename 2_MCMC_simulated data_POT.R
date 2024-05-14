library(expint)
library(mvtnorm)
library(rmutil) #levy
set.seed(432)

# --------------------------- Generating data  ---------------------------

d <- 4
N <- 500
sigma <- matrix(c(1,0.3,0.4,0.2,
                  0.3,1,0.5,0.3,
                  0.4,0.5,1,0.6,
                  0.2,0.3,0.6,1), ncol=d)
W <- matrix(c(1,0,0,0,
              1,1,0,0,
              0,1,1,0,
              1,0,0,1), nrow=d)
phi <- c(0.52,0.56,0.6,0.65)

g <- function(Z){
  x <- 1/(1-pnorm(Z))
  return(x)
}

Z <- array(NA,dim = c(N,d))
X <- array(NA,dim = c(N,d))
Y <- array(NA,dim = c(N,d))
S <- array(NA,dim = c(N,d))


for (i in 1:N) {
  Z[i,] <- rmvnorm(1, mean=c(0,0,0,0), sigma=sigma)
  S[i,] <- rlevy(d, m=0, s=1)
  R <- t(W)%*%S[i,]
  for (j in 1:d) {
    X[i,j] <- (R[j]^phi[j])*g(Z[i,j])
  }
}

tau <- 1

for (i in 1:N) {
  Y[i,] <- X[i,] + rnorm(4,0,tau)
}

L <- apply(Y, 2, function(x) quantile(x,0.9))
#L <- quantile(Y, 0.9)


# MCMC Utility functions -------------------------------------------------------

log.post.X <- function(X_j, phi, S_j, sigma){
  R = t(W) %*% S_j
  z = (R^phi)/X_j
  if(all(z > 0 & z < 1)){
    g_inv_vec = qnorm(1-z)
    log_d_ginv_vec = log(z/X_j) - dnorm( g_inv_vec, 0, 1)
  }else{
    return(-Inf)
  }
  lprob = mvtnorm::dmvnorm(as.vector(g_inv_vec), c(0,0,0,0), sigma, log = TRUE) + sum(log_d_ginv_vec)
  return(lprob)
}

log.post.Y <- function(X_j, L, t){
  lprob <- array(NA, dim = d)
  for (i in 1:d) {
    if(X_j[i] <= L[i])
      lprob[i] <- pnorm(L[j], mean = 0, sd = t^2, log.p = TRUE)
    else
      lprob[i] <- dnorm(X_j[i], mean = 0, sd = t^2, log = TRUE)
  }
  return(sum(lprob))
}


g.inv <- function(x,R,phi){
  out1 <- 1-(R^phi)/x
  suppressWarnings(out2 <- sapply(c(out1),qnorm))
  return(out2)
}

log.post.s <- function(x,phi,s,sigma){
  R <- t(W)%*%s
  Z_tmp <- g.inv(x,R,phi)
  out2 <- sapply(c(Z_tmp),dnorm)
  if(any(is.na(out2))){
    return(-Inf)
  }
  out3 <- ((R^phi)/x^2)/out2
  if(any(s<0)){
    out4 = 0
  }else{
    out4 <- sapply(c(s),dlevy)
  }
  out1 <- dmvnorm(Z_tmp, mean=c(0,0,0,0), sigma=sigma, log = TRUE)
  out <- out1 + sum(log(out3) + sum(log(out4)))
  return(out)
}

log.post.phi <- function(x,phi,s, sigma){
  if(any(phi<0.45 | phi >0.7)){
    return(-Inf)
  }
  R <- t(W)%*%s
  Z_tmp <- g.inv(x,R,phi)
  out2 <- sapply(c(Z_tmp),dnorm)
  if(any(is.na(out2))){
    return(-Inf)
  }
  out3 <- ((R^phi)/x^2)/out2
  # if(is.null(dim(Z_tmp))){
  #   out1 <- 0
  #   out <- out1 + sum(log(out3))
  #   return(out)
  # }
  out1 <- dmvnorm(Z_tmp, mean=c(0,0,0,0), sigma=sigma, log = TRUE)
  out <- out1 + sum(log(out3))
  return(out)
}

log.post.phi.total <- function(X,phi,S,sigma){
  n <- nrow(X)
  out <- array(0,n)
  for (i in 1:n) {
    out[i] <- log.post.phi(X[i,],phi,S[i,],sigma)
  }
  return(sum(out))
}

log.post.tau <- function(Y, X, tau, L){
  lprob <- array(NA, dim = dim(X))
  for (i in 1:nrow(X)) {
    for (j in 1:d) {
      if(X[i,j] <= L[j])
        lprob[i,j] <- pnorm(L[j], mean = 0, sd = tau^2, log.p = TRUE)
      else
        lprob[i,j] <- dnorm(X[i,j], mean = 0, sd = tau^2, log = TRUE)
    }
  }
  return(sum(lprob))
}


log.post.sigma <- function(X, phi, S, sigma){
  lprob <- log.post.x.total(X, phi, S, sigma)
  return(lprob)
}


log.post.x.total <- function(X, phi, S, sigma){
  lprob <- array(NA, dim = nrow(X))
  for (i in 1:nrow(X)) {
    lprob[i] <- log.post.X( X[j,], phi, S[i,], sigma)
  }
  return(sum(lprob))
}

cov.s.fun <- function(dat){
  cov.s.mat <- apply(dat, 2, cov)
  cov.s.list <- list()
  for( iter in 1:ncol(cov.s.mat)){
    cov.s.list[[iter]] <- matrix(cov.s.mat[,iter],4,4)
  }
  return(cov.s.list)
}





# Adaptive MCMC -----------------------------------------------------------
cov.s <- list()
cov.X <- list()
for( iter in 1:N){
  cov.s[[iter]] <- diag(4)
  cov.X[[iter]] <- diag(4)
}

cov.phi <- diag(4)
phi.initial <- c(0.52,0.56,0.6,0.65)
S.initial <- S
M <- 600
k <- 100

X.mc  <- array(NA,dim = c(M,N,4))
s.mc  <- array(NA,dim = c(M,N,4))
phi.mc <- array(NA,dim = c(M,4))
tau.mc <- array(NA,dim = M)
sigma.mc <- array(NA, dim = c(4,4,M))

acc.X <- array(0, dim = c(M, N))
acc.s <- array(0, dim = c(M, N))
acc.phi  <- array(0,M)
acc.tau  <- array(0,M)
acc.sigma  <- array(0,M)

r.X  <- array(NA,dim = c(M,N))
r.s  <- array(NA,dim = c(M,N))
r.phi  <- array(NA,M)
r.tau  <- array(NA,M)
r.sigma  <- array(NA,M)

d <- 4
c.0 <- 1
c.1 <- 0.8
r.opt <- 0.234
sig.phi <- 0.00001
sig.s <- array(2.4^2/d,N)
sig.X <- array(2.4^2/d,N)
t <- 1

for (i in 1:M) {
  print(paste0("print i is: ",i))
  if(i %% 100 == 0) cat("Done with", i, "iterations\n")

  # update X
  for (j in 1:N){
    X.t <- as.vector(X[j,] + rmvnorm(1,c(0,0,0,0),sig.X[j]*cov.X[[j]]) )
    r <- exp(log.post.X( X.t, phi, S[j,], sigma) - log.post.X( X[j,], phi, S[j,], sigma) + log.post.Y(X.t, L, tau) - log.post.Y(X[j,], L, tau))
    r.X[i,j] <- r
    if (is.na(r)){
      X.mc[i,j,] <- X[j,]
    }
    if (!is.na(r) & runif(1) < r){
      X[j,] <- X.t
      acc.X[i,j] <- 1
    }
    X.mc[i,j,] <- X[j,]
  }

  # update S
  for (j in 1:N){
    s.t <- as.vector(S[j,] + rmvnorm(1,c(0,0,0,0),sig.s[j]*cov.s[[j]]))
    r <- exp(log.post.s(X[j,], phi, s.t, sigma) - log.post.s(X[j,], phi, S[j,], sigma))
    r.s[i,j] <- r
    if (is.na(r)){
      s.mc[i,j,] <- S[j,]
    }
    if (!is.na(r) & runif(1) < r){
      S[j,] <- s.t
      acc.s[i,j] <- 1
    }
    s.mc[i,j,] <- S[j,]
  }

  # update phi
  phi.t <- as.vector( phi + rmvnorm(1,c(0,0,0,0),sig.phi*cov.phi))
  print(paste0("Phi is: ",phi.t))
  r.phi[i] <- exp(log.post.phi.total(X, phi.t, S, sigma)/N - log.post.phi.total(X, phi, S, sigma)/N)
  print(paste0("r.phi is: ",r.phi[i]))
  if (is.na(r.phi[i])){
    phi.mc[i,] <- phi
  }
  if (!is.na(r.phi[i]) & runif(1) < r.phi[i]){
    phi <- phi.t
    acc.phi[i] <- 1
  }
  phi.mc[i,] <- phi

  # update tau
  tau.t <- tau + runif(1, 0, 1)
  r <- exp(log.post.tau(Y,X, tau.t, L) - log.post.tau(Y,X, tau, L))
  print(paste0("r is: ",r))
  r.tau[i] <- r
  if (is.na(r)){
    tau.mc[i] <- tau
  }
  if( runif(1) < r)
    tau <- tau.t
  tau.mc[i] <- tau

  # update sigma
  eigen_values <- rep(-1,4)
  while (any(eigen_values < 0)) {
    sigma.t <- matrix(0, nrow = 4, ncol = 4)
    for (i in c(2,3,4,7,8,12)) {
      sigma.t[i] <- sigma[i] + runif(1,0,1)
    }
    sigma.t <- sigma.t + t(sigma.t)
    diag(sigma.t) <- 1
    eigen_values <- eigen(sigma.t)$values
  }
  #r <- exp(log.post.sigma(sigma.t) - log.post.sigma(sigma) + log.post.x.total(X, phi, S, sigma.t) - log.post.x.total(X, phi, S, sigma))
  r <- exp(log.post.sigma(X, phi, S, sigma.t) - log.post.sigma(X, phi, S, sigma))
  r.sigma[i] <- r
  if (is.na(r)){
    sigma.mc[,,i] <- sigma
  }
  if( runif(1) < r)
    sigma <- sigma.t
  sigma.mc[,,i] <- sigma


  #adaptive part
  if(i %% k == 0){
    gamma.1 <- 1/(t+3)^c.1
    gamma.2 <- c.0*gamma.1

    r.X.est <- apply(acc.X[(i-k+1):i,],2, function(x) sum(x,na.rm = TRUE)/k)
    cov.X.est <- cov.s.fun(X.mc[(i-k+1):i,,])
    cov.X.list <- list()
    for (iter in 1:length(cov.X.est)) {
      cov.X.list[[iter]] <- cov.X[[iter]] + gamma.1*(cov.X.est[[iter]] - cov.X[[iter]])
    }
    cov.X <- cov.X.list
    log.sig.sq.new.X <- log(sig.X^2) + gamma.2*(r.X.est-r.opt)
    sig.X <- sqrt(exp(log.sig.sq.new.X))

    r.phi.est <- sum(acc.phi[(i-k+1):i], na.rm = TRUE)/k
    cov.phi.est <- cov(phi.mc[(i-k+1):i, ])
    cov.phi <- cov.phi + gamma.1*(cov.phi.est - cov.phi)
    log.sig.sq.new.phi <- log(sig.phi^2) + gamma.2*(r.phi.est-r.opt)
    sig.phi <- sqrt(exp(log.sig.sq.new.phi))

    r.s.est <- apply(acc.s[(i-k+1):i,],2, function(x) sum(x,na.rm = TRUE)/k)
    cov.s.est <- cov.s.fun(s.mc[(i-k+1):i,,])
    cov.s.list <- list()
    for (iter in 1:length(cov.s.est)) {
      cov.s.list[[iter]] <- cov.s[[iter]] + gamma.1*(cov.s.est[[iter]] - cov.s[[iter]])
    }
    cov.s <- cov.s.list
    log.sig.sq.new.s <- log(sig.s^2) + gamma.2*(r.s.est-r.opt)
    sig.s <- sqrt(exp(log.sig.sq.new.s))

    t <- t + 1
  }
}


model_name <- 'phi_mcmc'

pdf(paste0("Posterior of phi model_",model_name,".pdf"),height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=FALSE)
plot(density(phi.mc[,1]), main = expression(phi[1]))
plot(density(phi.mc[,2]), main = expression(phi[2]))
plot(density(phi.mc[,3]), main = expression(phi[3]))
plot(density(phi.mc[,4]), main = expression(phi[4]))
dev.off()

pdf(paste0("phi with delta model_",model_name,".pdf"),height=6,width=7.5)
  par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=FALSE)
  plot(1:M, phi.mc[,1], type = "l", ylab = expression(phi[1]))
  # lines(fit@sim$samples[[2]]$`phi[1]`,col = "blue")
  # lines(fit@sim$samples[[3]]$`phi[1]`,col = "green")
  # lines(fit@sim$samples[[4]]$`phi[1]`,col = "magenta4")

  plot(1:M, phi.mc[,2], type = "l", ylab = expression(phi[2]))
  # lines(fit@sim$samples[[2]]$`phi[2]`,col = "blue")
  # lines(fit@sim$samples[[3]]$`phi[2]`,col = "green")
  # lines(fit@sim$samples[[4]]$`phi[2]`,col = "magenta4")

  plot(1:M, phi.mc[,3], type = "l", ylab = expression(phi[3]))
  # lines(fit@sim$samples[[2]]$`phi[3]`,col = "blue")
  # lines(fit@sim$samples[[3]]$`phi[3]`,col = "green")
  # lines(fit@sim$samples[[4]]$`phi[3]`,col = "magenta4")

  plot(1:M, phi.mc[,4], type = "l", ylab = expression(phi[4]))
  # lines(fit@sim$samples[[2]]$`phi[4]`,col = "blue")
  # lines(fit@sim$samples[[3]]$`phi[4]`,col = "green")
  # lines(fit@sim$samples[[4]]$`phi[4]`,col = "magenta4")
  #mtext(paste0("phi with delta ",delta), side = 3,line = - 2,outer = TRUE)
dev.off()

pdf(paste0("Posterior of S model_",model_name,".pdf"),height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=FALSE)
plot(density(s.mc[,2,1]), main = "s[2,1]")
plot(density(s.mc[,2,2]), main = "s[2,2]")
plot(density(s.mc[,2,3]), main = "s[2,3]")
plot(density(s.mc[,2,4]), main = "s[2,4]")
dev.off()

pdf(paste0("Posterior of sigma model_",model_name,".pdf"),height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=FALSE)
plot(density(sigma.mc[1,2,]), main = "sigma[1,2]")
plot(density(sigma.mc[1,3,]), main = "sigma[1,3]")
plot(density(sigma.mc[1,4,]), main = "sigma[1,4]")
plot(density(sigma.mc[2,3,]), main = "sigma[2,3]")
dev.off()

pdf(paste0("Posterior of X[2,] model_",model_name,".pdf"),height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=FALSE)
plot(density(X.mc[2,1,]), main = "X[2,1]")
plot(density(X.mc[2,2,]), main = "X[2,2]")
plot(density(X.mc[2,3,]), main = "X[2,3]")
plot(density(X.mc[2,4,]), main = "X[2,4]")
dev.off()

pdf(paste0("Posterior of tau model_",model_name,".pdf"),height=6,width=7.5)
plot(density(tau.mc), main = expression(tau))
abline(v = tau, col = "red")
dev.off()


pdf(paste0("s with delta model_",model_name,".pdf"),height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=FALSE)
plot(s.mc[,2,1], type = "l",  ylab = "s[2,1]")
# lines(fit@sim$samples[[2]]$`s[2,1]`,col = "blue")
# lines(fit@sim$samples[[3]]$`s[2,1]`,col = "green")
# lines(fit@sim$samples[[4]]$`s[2,1]`,col = "magenta4")
abline(h=S.initial[2,1], col = "red")

plot(s.mc[,2,2], type = "l", ylab = "s[2,2]")
# lines(fit@sim$samples[[2]]$`s[2,2]`,col = "blue")
# lines(fit@sim$samples[[3]]$`s[2,2]`,col = "green")
# lines(fit@sim$samples[[4]]$`s[2,2]`,col = "magenta4")
abline(h=S.initial[2,2], col = "red")

plot(s.mc[,2,3], type = "l", ylab = "s[2,3]")
# lines(fit@sim$samples[[2]]$`s[2,3]`,col = "blue")
# lines(fit@sim$samples[[3]]$`s[2,3]`,col = "green")
# lines(fit@sim$samples[[4]]$`s[2,3]`,col = "magenta4")
abline(h=S.initial[2,3], col = "red")

plot(s.mc[,2,4], type = "l", ylab = "s[2,4]")
# lines(fit@sim$samples[[2]]$`s[2,4]`,col = "blue")
# lines(fit@sim$samples[[3]]$`s[2,4]`,col = "green")
# lines(fit@sim$samples[[4]]$`s[2,4]`,col = "magenta4")
abline(h=S.initial[2,4], col = "red")
dev.off()


pdf(paste0("sigma with delta model_",model_name,".pdf"),height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=FALSE)
plot(sigma.mc[1,2,], type = "l",  ylab = "sigma[1,2]")
# lines(fit@sim$samples[[2]]$`sigma[1,2]`,col = "blue")
# lines(fit@sim$samples[[3]]$`sigma[1,2]`,col = "green")
# lines(fit@sim$samples[[4]]$`sigma[1,2]`,col = "magenta4")
abline(h=sigma[1,2], col = "red")

plot(sigma.mc[1,2,], type = "l", ylab = "sigma[1,3]")
# lines(fit@sim$samples[[2]]$`sigma[1,3]`,col = "blue")
# lines(fit@sim$samples[[3]]$`sigma[1,3]`,col = "green")
# lines(fit@sim$samples[[4]]$`sigma[1,3]`,col = "magenta4")
abline(h=sigma[1,3], col = "red")

plot(sigma.mc[1,2,], type = "l", ylab = "sigma[1,4]")
# lines(fit@sim$samples[[2]]$`sigma[1,4]`,col = "blue")
# lines(fit@sim$samples[[3]]$`sigma[1,4]`,col = "green")
# lines(fit@sim$samples[[4]]$`sigma[1,4]`,col = "magenta4")
abline(h=sigma[1,4], col = "red")

plot(sigma.mc[1,2,], type = "l", ylab = "sigma[2,3]")
# lines(fit@sim$samples[[2]]$`sigma[2,3]`,col = "blue")
# lines(fit@sim$samples[[3]]$`sigma[2,3]`,col = "green")
# lines(fit@sim$samples[[4]]$`sigma[2,3]`,col = "magenta4")
abline(h=sigma[2,3], col = "red")
# mtext(paste0("sigma with delta ", delta), side = 3,line = - 2,outer = TRUE)
dev.off()

pdf(paste0("X with delta model_",model_name,".pdf"),height=6,width=7.5)
par(mfrow = c(2,2),mar=c(3,4,5,3),xpd=FALSE)
plot(X.mc[2,1,], type = "l",  ylab = "X[2,1]")
# lines(fit@sim$samples[[2]]$`X[2,1]`,col = "blue")
# lines(fit@sim$samples[[3]]$`X[2,1]`,col = "green")
# lines(fit@sim$samples[[4]]$`X[2,1]`,col = "magenta4")
abline(h=X[2,1], col = "red")

plot(X.mc[2,2,], type = "l", ylab = "X[2,2]")
# lines(fit@sim$samples[[2]]$`X[2,2]`,col = "blue")
# lines(fit@sim$samples[[3]]$`X[2,2]`,col = "green")
# lines(fit@sim$samples[[4]]$`X[2,2]`,col = "magenta4")
abline(h=X[2,2], col = "red")

plot(X.mc[2,3,], type = "l", ylab = "X[2,3]")
# lines(fit@sim$samples[[2]]$`X[2,3]`,col = "blue")
# lines(fit@sim$samples[[3]]$`X[2,3]`,col = "green")
# lines(fit@sim$samples[[4]]$`X[2,3]`,col = "magenta4")
abline(h=X[2,3], col = "red")

plot(X.mc[2,4,], type = "l", ylab = "X[2,4]")
# lines(fit@sim$samples[[2]]$`X[2,4]`,col = "blue")
# lines(fit@sim$samples[[3]]$`X[2,4]`,col = "green")
# lines(fit@sim$samples[[4]]$`X[2,4]`,col = "magenta4")
abline(h=X[2,4], col = "red")
# mtext(paste0("X[,2,] with delta ", delta), side = 3,line = - 2,outer = TRUE)
dev.off()

pdf(paste0("tau  model_",model_name,".pdf"),height=6,width=7.5)
plot(tau.mc, type = "l",  ylab = "tau")
# lines(fit@sim$samples[[2]]$`tau`,col = "blue")
# lines(fit@sim$samples[[3]]$`tau`,col = "green")
# lines(fit@sim$samples[[4]]$`tau`,col = "magenta4")
abline(h=tau, col = "red")
# mtext(paste0("tau with delta ", delta), side = 3,line = - 2,outer = TRUE)
dev.off()


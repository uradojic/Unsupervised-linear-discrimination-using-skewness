# R code for the paper
# "Unsupervised linear discrimination using skewness"
# by US IN SOME ORDER

# Notations correspond to the paper notation

# Required R packages
# - mvtnorm
# - ICtest
# - MASS
# - parallel
# - ggplot2

# general two component mixture model
# input:
# n: number of random samples desired
# alpha1: mixing proportion for first group
# h: p-vector difference from the origin for the locations, cannot be zero
# Sigma: full rank p times p matrix
#
# output:
# n x p matrix with random samples from the distribution
# attribute G gives the cluster membership
library(mvtnorm)
binGMM2 <- function(n, alpha1, h, Sigma = diag(length(h)))
{
  n1 <- floor(n*alpha1)
  n2 <- n-n1
  alpha2 <- 1-alpha1
  mu1 <- -alpha2*h
  mu2 <- alpha1*h
  x1 <- rmvnorm(n1, mean=mu1, sigma=Sigma)
  x2 <- rmvnorm(n2, mean=mu2, sigma=Sigma)
  X <- rbind(x1,x2)
  attr(X, which = "G") <- factor(rep(LETTERS[1:2],c(n1,n2) ))
  X
}

TEST <- binGMM2(100, 0.4, 1:3) 
attr(x = TEST, which = "G")


# Function to compute C2 (cov with different divisor)
# input:
# X: n x p data matrix
# output:
# p x p matrix 
C2 <- function(X)
{
  n <- nrow(X)
  X <- sweep(X, 2, colMeans(X), "-")
  crossprod(X)/n
}

# Function to compute c3 
# input:
# X: n x p data matrix
# output:
# p vector
c3 <- function(X)
{
  X <- sweep(X, 2, colMeans(X), "-")
  W <- rowSums(X^2)
  colMeans(X*W)
}

# Function for BETA 
# input:
# alpha: probability of group 1
# output:
# numeric value
BETA <- function(alpha)
{
  alpha*(1-alpha)
}

# Function for GAMMA
# input:
# alpha: probability of group 1
# output:
# numeric value
GAMMA <- function(alpha)
{
  alpha2 <- 1 - alpha
  alphas <- c(alpha,alpha2)
  max(alphas)- min(alphas)
}


# Function to compute C for affine equivariant moment estimator
# input:
# alpha: probability of group 1
# value: tau (function of h and Sigma)
# p: dimension of data
# output:
# numeric value
C_R <- function(alpha, tau, p)
{
  beta <- BETA(alpha)
  p1 <- 2*(p+1)*(1+beta*tau)^4
  p2 <- beta^2*(1-4*beta)*tau^3
  p3 <- (1+beta*tau)
  p4 <- beta*tau^2 + 6*beta*tau + 2
  
  p1/p2 + p3*p4/p2
}

# Function to compute C for TOBI
# input:
# alpha: probability of group 1
# value: tau (function of h and Sigma)
# output:
# numeric value
C_L <- function(alpha, tau)
{
  beta <- BETA(alpha)
  p1 <- (1+beta*tau)
  p2 <- beta*tau^2 + 6*beta*tau + 2
  p3 <- beta^2*(1-4*beta)*tau^3
  p1 * p2/p3
}

# Function to compute C for LDA
# input:
# alpha: probability of group 1
# value: tau (function of h and Sigma)
# output:
# numeric value
C_LDA <- function(alpha, tau)
{
  beta <- BETA(alpha)
  p1 <- (1+beta*tau)
  p2 <- beta*tau
  p1/p2
}

# Moment based estimator theta_m

# input:
# X: n x p data matrix
# alpha1: mixing proportion for first group
# output:
# p-vector with unit norm
ThetaM <- function(X, alpha1)
{
  n <- nrow(X)
  p <- ncol(X)
  gamma <- GAMMA(alpha1)
  beta <- BETA(alpha1)
  center <- colMeans(X)
  X <- sweep(X, 2, center, "-")
  C2 <- crossprod(X)/n
  W <- rowSums(X^2)
  c3 <- colMeans(X*W)
  n.c3.2 <- sum(c3^2)
  part1 <- (C2- beta^(1/3)*gamma^(-2/3)*n.c3.2^(-2/3)*tcrossprod(c3))
  part2 <- beta^(-1/3)*gamma^(-1/3)*n.c3.2^(-1/3) * c3
  THETA <- solve(part1, part2)
  THETA
  THETA / sqrt(sum(THETA^2))
}


# Affine equivariant moment based estimator theta_R
# input:
# X: n x p data matrix
# output:
# p-vector with unit norm
ThetaR <- function(X)
{
  n <- nrow(X)
  center <- colMeans(X)
  X <- sweep(X, 2, center, "-")
  C2 <- crossprod(X)/n
  EVD.C2 <- eigen(C2, symmetric = TRUE)
  C2.inv.sqrt <- EVD.C2$vectors %*% tcrossprod(diag((1/EVD.C2$values)^0.5), 
                                               EVD.C2$vectors)
  Xst <- tcrossprod(X, C2.inv.sqrt)
  W <- rowSums(Xst^2)
  c3 <- colMeans(Xst*W)
  THETA <- c3 %*% C2.inv.sqrt  
  as.vector(THETA / sqrt(sum(THETA^2)))
}


# Tobi estimator, theta_L
# input:
# X: n x p data matrix
# output:
# p-vector with unit norm
ThetaL <- function(X)
{
  n <- nrow(X)
  p <- ncol(X)
  center <- colMeans(X)
  X <- sweep(X, 2, center, "-")
  C2 <- crossprod(X)/n
  EVD.C2 <- eigen(C2, symmetric = TRUE)
  C2.inv.sqrt <- EVD.C2$vectors %*% tcrossprod(diag((1/EVD.C2$values)^0.5), 
                                               EVD.C2$vectors)
  Xst <- tcrossprod(X, C2.inv.sqrt)
  
  tobi <- matrix(0, p, p)
  for(i in 1:p){
    t_k <- crossprod(Xst,sweep(Xst, 1, Xst[, i], "*"))/n
    tobi <- tobi + crossprod(t_k)
  }
  tobi.evd <- eigen(tobi, symmetric = TRUE)
  theta <- C2.inv.sqrt %*% tobi.evd$vectors[,1] 
  as.numeric(theta / sqrt(sum(theta^2)))
}

# 3-Jade estimator theta_J
# input:
# X: n x p data matrix
# init: initial starting value  
# eps: tolerance limit
# iter: maxiter
# output:
# p-vector with unit norm

ThetaJ <- function(X, init =rep(1, ncol(X)), eps = 1e-06, maxiter = 100)
{
  n <- nrow(X)
  p <- ncol(X)
  center <- colMeans(X)
  X <- sweep(X, 2, center, "-")
  C2 <- crossprod(X)/n
  EVD.C2 <- eigen(C2, symmetric = TRUE)
  C2.inv.sqrt <- EVD.C2$vectors %*% tcrossprod(diag((1/EVD.C2$values)^0.5), 
                                               EVD.C2$vectors)
  Xst <- tcrossprod(X, C2.inv.sqrt)
  
  u_0 <- init/sqrt(sum(init^2))
  diff <- Inf
  iter <- 1
  while(diff > eps)
  {
    if (iter == maxiter) 
      stop("maxiter reached without convergence")
    u_new <- rep(0,p)
    for(i in 1:p){
      t_k <- crossprod(Xst,sweep(Xst, 1, Xst[, i], "*"))/n
      u_i <- t(u_0) %*% t_k %*% (u_0) %*% t(u_0) %*% t_k 
      u_new <- u_new + u_i
    }
    u_new <- u_new/sqrt(sum(u_new^2))
    diff <- sqrt(sum((u_new-u_0)^2))
    u_0 <- as.numeric(u_new)
    iter=iter+1
  }
  theta <-  u_0 %*% C2.inv.sqrt 
  as.numeric(theta / sqrt(sum(theta^2)))
}

# skewnness based PP-estimator theta_P
# wrapper for function NGPP from package ICtest
# input:
# X: n x p data matrix
# output:
# p-vector with unit norm
library(ICtest)
ThetaP <- function(X)
{
  OP <- NGPP(X, k=1, alpha=1.0, method="defl", maxiter=1000)$W[1,]
  as.numeric(OP / sqrt(sum(OP^2))) 
}


# LDA estimator 
# wrapper for function lda from package MASS
# input:
# X: n x p data matrix
# G: n-vector with specifying the two group labels
# output:
# p-vector with unit norm
library(MASS)
ThetaLDA <- function(X, G)
{
  LDA <- lda(G ~ X)
  u_LDA<-LDA$scaling
  as.numeric(u_LDA/ sqrt(sum(u_LDA^2)))
}

# Some examples
set.seed(1)
h <- c(2,2,2)
X <- binGMM2(2000, 0.15, h)

OM <- ThetaM(X, 0.15)
OR <- ThetaR(X)
OL <- ThetaL(X)
OJ <- ThetaJ(X, init=OL)
OP <- ThetaP(X)
OLDA <- ThetaLDA(X, G=attr(X, "G"))
OM
OR
OL
OJ
OLDA
true <- h / sqrt(sum(h^2))

# performance measure
# function to compute angle between to unit vectors
MSI <- function(a,b) abs(sum(a*b))  / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) )
MSI1 <- MSI(true,OM)
MSI2 <- MSI(true,OR)
MSI3 <- MSI(true,OL)
MSI4 <- MSI(true,OJ)
MSI5 <- MSI(true,OP)
MSI6 <- MSI(true,OLDA)
MSI1
MSI2
MSI3
MSI4
MSI5
MSI6

###############################################################################
# Functions for Simulation 1:
###############################################################################

# Function to find an orthogonal p-vector having unit length
# input:
# m: p vector
# output:
# unit p-vector orthogonal to m
find_t <- function(m) {
  p <- length(m)
  u <- rep(0,p)
  idx <- which(m != 0)[1]
  if (idx < p) {
    u[idx + 1] <- 1
  } else {
    u[1] <- 1
  }
  
  proj_m_u <- sum(u * m) / sum(m * m)
  u_orth <- u - proj_m_u * m
  norm_u_orth <- sqrt(sum(u_orth^2))
  t <- u_orth / norm_u_orth
  t
}

# Function to for a single iteration to get C
# input:
# n: integer, sample size
# a: mixing proportion in [0,1]
# m: p-vector, value for h
# t: p-vector, unit length vector orthogonal to m.
# output:
# vector for C, order
# LDA, theta_R, theta_L, theta_P and theta_J
SingleIterationC <- function(n,a,m,t)
{
  X <- binGMM(n=n, alpha1=a, h=m)
  OR <- ThetaR(X)
  OL <- ThetaL(X)
  OJ <- tryCatch(ThetaJ(X), error = function(e) NA)
  OP <- tryCatch(NGPP(X, k=1, alpha=1.0, method="defl", maxiter=1000)$W[1,], error = function(e) NA)
  OP <- tryCatch(OP/sqrt(sum(OP^2)), error = function(e) NA)
  G <- attr(X, "G")
  LDA <- lda(G ~ X)
  u_LDA<-LDA$scaling
  OLDA <-  as.numeric(u_LDA/ sqrt(sum(u_LDA^2)))
  CLDA <- (t%*%OLDA)
  CL <- (t%*%OL)
  CR <- (t%*%OR)
  CJ <- tryCatch((t%*%OJ), error = function(e) NA)
  CP <- tryCatch((t%*%OP), error = function(e) NA)
  c(CLDA,CR,CL,CP,CJ)
}

# checks
ms <- rep(2,3)
ms %*% ms
ts <- find_t(ms)  
ts %*% ms

# Obtaining estimates based on 10000 repetitions
# example for a=0.3. For other a values change the value accordingly
# is slow
mrep <- 10000
set.seed(1)
n_500_03 <- replicate(mrep,SingleIterationC(n=500, a=0.3, m=ms, t=ts))
n_1000_03 <- replicate(mrep,SingleIterationC(n=1000, a=0.3, m=ms, t=ts))
n_2000_03 <- replicate(mrep,SingleIterationC(n=2000, a=0.3, m=ms, t=ts))
n_4000_03 <- replicate(mrep,SingleIterationC(n=4000, a=0.3, m=ms, t=ts))
n_8000_03 <- replicate(mrep,SingleIterationC(n=8000, a=0.3, m=ms, t=ts))
n_16000_03 <- replicate(mrep,SingleIterationC(n=16000, a=0.3, m=ms, t=ts))
n_32000_03 <- replicate(mrep,SingleIterationC(n=32000, a=0.3, m=ms, t=ts))
n_64000_03 <- replicate(mrep,SingleIterationC(n=64000, a=0.3, m=ms, t=ts))

n_500_03_b <- t(n_500_03)
n_1000_03_b <- t(n_1000_03)
n_2000_03_b <- t(n_2000_03)
n_4000_03_b <- t(n_4000_03)
n_8000_03_b <- t(n_8000_03)
n_16000_03_b <- t(n_16000_03) 
n_32000_03_b <- t(n_32000_03)
n_64000_03_b <- t(n_64000_03)

obs_500_03 <- colSums(is.na(n_500_03_b))
obs_1000_03 <- colSums(is.na(n_1000_03_b))
obs_2000_03 <- colSums(is.na(n_2000_03_b))
obs_4000_03 <- colSums(is.na(n_4000_03_b))
obs_8000_03 <- colSums(is.na(n_8000_03_b))
obs_16000_03 <- colSums(is.na(n_16000_03_b))
obs_32000_03 <- colSums(is.na(n_32000_03_b))
obs_64000_03 <- colSums(is.na(n_64000_03_b))

obs_500_03 
obs_1000_03 
obs_2000_03
obs_4000_03 
obs_8000_03 
obs_16000_03
obs_32000_03 
obs_64000_03 

VARS_500_03 <- apply(n_500_03_b, 2, var, na.rm = TRUE)
VARS_1000_03 <- apply(n_1000_03_b, 2, var, na.rm = TRUE)
VARS_2000_03 <- apply(n_2000_03_b, 2, var, na.rm = TRUE)
VARS_4000_03 <- apply(n_4000_03_b, 2, var, na.rm = TRUE)
VARS_8000_03 <- apply(n_8000_03_b, 2, var, na.rm = TRUE)
VARS_16000_03 <- apply(n_16000_03_b, 2, var, na.rm = TRUE)
VARS_32000_03 <- apply(n_32000_03_b, 2, var, na.rm = TRUE)
VARS_64000_03 <- apply(n_64000_03_b, 2, var, na.rm = TRUE)

RES_500_03 <- 500*VARS_500_03
RES_1000_03 <- 1000*VARS_1000_03
RES_2000_03 <- 2000*VARS_2000_03
RES_4000_03 <- 4000*VARS_4000_03
RES_8000_03 <- 8000*VARS_8000_03
RES_16000_03 <- 16000*VARS_16000_03
RES_32000_03 <- 32000*VARS_32000_03
RES_64000_03 <- 64000*VARS_64000_03

RES_03 <- rbind(
  RES_500_03,
  RES_1000_03,
  RES_2000_03,
  RES_4000_03,
  RES_8000_03,
  RES_16000_03,
  RES_32000_03,
  RES_64000_03
)
colnames(RES_03) <- c("CLDA","CR","CL","CP","CJ")
# final result
RES_03
# true values
C_R(alpha=0.3, tau =sum(ms^2), p=3)
C_L(alpha=0.3, tau =sum(ms^2))
C_LDA(alpha=0.3, tau =sum(ms^2))



###############################################################################
# Functions for Simulation 2:
###############################################################################

# Function to generate a random p times p matrix A
# from N(0,1) which makes sure that AA' can be inverted

rSigma <- function(p, kappa_limit = 100000) {
  while (TRUE) {
    A <- matrix(rnorm(p^2), nrow = p)
    Sigma <- A %*% t(A)
    cond_number <- kappa(Sigma)
    if (cond_number < kappa_limit) {
      Sigma.inv <- solve(Sigma)
      break 
    }
  }
  return(list(A = A, Sigma = Sigma, Sigma.inv = Sigma.inv))
}
# input:
# p: positive integer
# kappa_limit: limit until when it is considered that AA' is invertible 
# output:
# list with elemets:
# - A: p x p matrix
# - Sigma: p x p matrix (AA')
# - Sigma.inv: p x p matrix (inverse of Sigma)
rSigma(3)


# Function for simulations for finite sample performance evaluations
# input:
# m: positive integer, number of repetitions
# n: positive integer, number of sample sizes
# alphas: vector of alpha values (mixing probability, need all between 0 and 1)
# taus: vector of tau values (need all be positive)
# p: integer, giving data dimension
# seed: setting the seed for the simulation
# output:
# data frame with columns:
# - N=n, sample size
# - p=p, dimension of the data 
# - Sigma, character value, "I", or "S", if data was mixed with A (S) or not (I)
# - Alpha=alpha, mixing proportion
# - Tau=tau, 
# - METHOD=method, factor with levels "LDA","ThetaM","ThetaR","ThetaL","ThetaJ","ThetaP"
# - MSI = MSI value - will return NA when method did not converge
simu3momentsSI <- function(m, n, alphas=c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),
                           taus=1:20, p=3, seed = 1)
{
  n_alpha <- length(alphas)
  n_tau <- length(taus)
  n_comp <- m*n_alpha*n_tau 
  
  resM <- numeric(n_comp)
  resR <- numeric(n_comp)
  resL <- numeric(n_comp)
  resP <- numeric(n_comp)
  resJ <- numeric(n_comp)
  resLDA <- numeric(n_comp)
  alpha <- numeric(n_comp)
  tau <- numeric(n_comp)
  
  resM_I <- numeric(n_comp)
  resR_I <- numeric(n_comp)
  resL_I <- numeric(n_comp)
  resP_I <- numeric(n_comp)
  resJ_I <- numeric(n_comp)
  resLDA_I <- numeric(n_comp)
  
  
  set.seed(seed)
  iter <- 0
  a <- rep(1,p)
  
  for (aa in alphas)
  {
    for (tt in taus)
    {
      h.tt <- sqrt(tt)*a/sqrt(p)
      
      for (mm in 1:m)
      {
        iter <- iter + 1
        RAND <- rSigma(p=p)
        X <- binGMM2(n, aa, h.tt, Sigma=diag(p))
        X <- sweep(X,2,colMeans(X),"-")
        AX <- X %*% (RAND$A)
        
        ThetaTrue.x <- diag(p) %*% h.tt  
        ThetaTrue.x <- ThetaTrue.x / sqrt(sum(ThetaTrue.x^2))
        
        
        ThetaTrue.ax <- solve(RAND$A)%*%c(ThetaTrue.x)
        ThetaTrue.ax <- ThetaTrue.ax / sqrt(sum(ThetaTrue.ax^2))
        
        
        resM[iter] <- tryCatch(MSI(ThetaTrue.ax,ThetaM(AX,aa)), error = function(e) NA)
        resR[iter] <- tryCatch(MSI(ThetaTrue.ax,ThetaR(AX)), error = function(e) NA)
        TL <- ThetaL(AX)
        resL[iter] <- tryCatch(MSI(ThetaTrue.ax,TL), error = function(e) NA)
        resP[iter] <- tryCatch(MSI(ThetaTrue.ax,ThetaP(AX)), error = function(e) NA)
        resJ[iter] <- tryCatch(MSI(ThetaTrue.ax,ThetaJ(AX,init=TL)), error = function(e) NA)
        resLDA[iter] <- tryCatch(MSI(ThetaTrue.ax,ThetaLDA(AX,attr(X, "G"))), error = function(e) NA)
        
        resM_I[iter] <- tryCatch(MSI(ThetaTrue.x,ThetaM(X,aa)), error = function(e) NA)
        resR_I[iter] <- tryCatch(MSI(ThetaTrue.x,ThetaR(X)), error = function(e) NA)
        TL <- ThetaL(X)
        resL_I[iter] <- tryCatch(MSI(ThetaTrue.x,TL), error = function(e) NA)
        resP_I[iter] <- tryCatch(MSI(ThetaTrue.x,ThetaP(X)), error = function(e) NA)
        resJ_I[iter] <- tryCatch(MSI(ThetaTrue.x,ThetaJ(X,init=TL)), error = function(e) NA)
        resLDA_I[iter] <- tryCatch(MSI(ThetaTrue.x,ThetaLDA(X,attr(X, "G"))), error = function(e) NA)
        
        alpha[iter] <- aa
        tau[iter] <- tt
      }
    }
  }
  method= factor(rep(c("LDA","ThetaM","ThetaR","ThetaL","ThetaJ","ThetaP"), each=n_comp))
  RES_S <- data.frame(N=n, p=p,Sigma="S", Alpha=alpha, Tau=tau, METHOD=method,
                      MSI = c(resLDA,resM,resR,resL,resJ,resP))
  RES_I <- data.frame(N=n, p=p,Sigma="I", Alpha=alpha, Tau=tau, METHOD=method,
                      MSI = c(resLDA_I,resM_I,resR_I,resL_I,resJ_I,resP_I))
  rbind(RES_I,RES_S)
}


TEST <- simu3momentsSI(3, 2000, alphas=c(0.05,0.1,0.15),
                       taus=4:7, p=3, seed = 1)
# to reproduce simulation in paper chose approbriate settings. Is slow!

# Function for simulations for finite sample performance evaluations with parallel processing 
# input:
# m: positive integer, number of repetitions
# n: positive integer, number of sample sizes
# alphas: vector of alpha values (mixing probability, need all between 0 and 1)
# taus: vector of tau values (need all be positive)
# p: integer, giving data dimension
# seed: setting the seed for the simulation
# ncores: number of cores used in for parallel computations
# output:
# data frame with columns:
# - N=n, sample size
# - p=p, dimension of the data 
# - Sigma, character value, "I", or "S", if data was mixed with A (S) or not (I)
# - Alpha=alpha, mixing proportion
# - Tau=tau, 
# - METHOD=method, factor with levels "LDA","ThetaM","ThetaR","ThetaL","ThetaJ","ThetaP"
# - MSI = MSI value - will return NA when method did not converge
library(parallel)                                   
simu3momentsSIpar <- function(m, n, alphas = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45),
                           taus = 1:20, p = 3, seed = 1, ncores = 2) {
  
  n_alpha <- length(alphas)
  n_tau <- length(taus)
  n_comp <- m * n_alpha * n_tau
  
  resM <- numeric(n_comp)
  resR <- numeric(n_comp)
  resL <- numeric(n_comp)
  resP <- numeric(n_comp)
  resJ <- numeric(n_comp)
  resLDA <- numeric(n_comp)
  alpha <- numeric(n_comp)
  tau <- numeric(n_comp)
  
  resM_I <- numeric(n_comp)
  resR_I <- numeric(n_comp)
  resL_I <- numeric(n_comp)
  resP_I <- numeric(n_comp)
  resJ_I <- numeric(n_comp)
  resLDA_I <- numeric(n_comp)
  
  iter <- 0
  a <- rep(1, p)
  
  # Helper function for a single iteration
  run_single_iteration <- function(mm, aa, tt, p, n, a) {
    h.tt <- sqrt(tt) * a / sqrt(p)
    RAND <- rSigma(p = p)
    X <- binGMM2(n, aa, h.tt, Sigma = diag(p))
    X <- sweep(X, 2, colMeans(X), "-")
    AX <- X %*% (RAND$A)
    
    ThetaTrue.x <- diag(p) %*% h.tt
    ThetaTrue.x <- ThetaTrue.x / sqrt(sum(ThetaTrue.x^2))
    
    ThetaTrue.ax <- solve(RAND$A) %*% c(ThetaTrue.x)
    ThetaTrue.ax <- ThetaTrue.ax / sqrt(sum(ThetaTrue.ax^2))
    
    results <- list(
      resM = NA, resR = NA, resL = NA, resP = NA, resJ = NA, resLDA = NA,
      resM_I = NA, resR_I = NA, resL_I = NA, resP_I = NA, resJ_I = NA, resLDA_I = NA,
      alpha = aa, tau = tt
    )
    
    # Calculate the statistics
    results$resM <- tryCatch(MSI(ThetaTrue.ax, ThetaM(AX, aa)), error = function(e) NA)
    results$resR <- tryCatch(MSI(ThetaTrue.ax, ThetaR(AX)), error = function(e) NA)
    TL <- ThetaL(AX)
    results$resL <- tryCatch(MSI(ThetaTrue.ax, TL), error = function(e) NA)
    results$resP <- tryCatch(MSI(ThetaTrue.ax, ThetaP(AX)), error = function(e) NA)
    results$resJ <- tryCatch(MSI(ThetaTrue.ax, ThetaJ(AX, init = TL)), error = function(e) NA)
    results$resLDA <- tryCatch(MSI(ThetaTrue.ax, ThetaLDA(AX, attr(X, "G"))), error = function(e) NA)
    
    # Results for identity matrix case
    results$resM_I <- tryCatch(MSI(ThetaTrue.x, ThetaM(X, aa)), error = function(e) NA)
    results$resR_I <- tryCatch(MSI(ThetaTrue.x, ThetaR(X)), error = function(e) NA)
    TL <- ThetaL(X)
    results$resL_I <- tryCatch(MSI(ThetaTrue.x, TL), error = function(e) NA)
    results$resP_I <- tryCatch(MSI(ThetaTrue.x, ThetaP(X)), error = function(e) NA)
    results$resJ_I <- tryCatch(MSI(ThetaTrue.x, ThetaJ(X, init = TL)), error = function(e) NA)
    results$resLDA_I <- tryCatch(MSI(ThetaTrue.x, ThetaLDA(X, attr(X, "G"))), error = function(e) NA)
    
    return(results)
  }
  
  # Create a cluster for parallel processing
  cl <- makeCluster(ncores)
  
  # Load the required packages on each worker
  clusterEvalQ(cl, {
    library(MASS)
    library(ICtest)
    library(mvtnorm)
  })
  
  # Set reproducible random seeds for each worker
  clusterSetRNGStream(cl, seed)
  
  # Export the necessary variables and functions to each worker in the cluster
  clusterExport(cl, c("rSigma", "binGMM2", "MSI", "ThetaM", "ThetaR", "ThetaL", 
                      "ThetaP", "ThetaJ", "ThetaLDA", "sweep", "solve", "diag",
                      "BETA","GAMMA","C2","c3"))
  
  # Iterate over alpha and tau
  for (aa in alphas) {
    for (tt in taus) {
      # Create tasks for parallel processing (loop over m)
      iter_results <- parLapply(cl, 1:m, function(mm) {
        run_single_iteration(mm, aa, tt, p, n, a)
      })
      
      # Store the results from parallel computation
      for (i in 1:length(iter_results)) {
        iter <- iter + 1
        
        resM[iter] <- iter_results[[i]]$resM
        resR[iter] <- iter_results[[i]]$resR
        resL[iter] <- iter_results[[i]]$resL
        resP[iter] <- iter_results[[i]]$resP
        resJ[iter] <- iter_results[[i]]$resJ
        resLDA[iter] <- iter_results[[i]]$resLDA
        
        resM_I[iter] <- iter_results[[i]]$resM_I
        resR_I[iter] <- iter_results[[i]]$resR_I
        resL_I[iter] <- iter_results[[i]]$resL_I
        resP_I[iter] <- iter_results[[i]]$resP_I
        resJ_I[iter] <- iter_results[[i]]$resJ_I
        resLDA_I[iter] <- iter_results[[i]]$resLDA_I
        
        alpha[iter] <- iter_results[[i]]$alpha
        tau[iter] <- iter_results[[i]]$tau
      }
    }
  }
  
  # Stop the cluster after processing
  stopCluster(cl)
  
  method <- factor(rep(c("LDA", "ThetaM", "ThetaR", "ThetaL", "ThetaJ", "ThetaP"), each = n_comp))
  RES_S <- data.frame(N = n, p = p, Sigma = "S", Alpha = alpha, Tau = tau, METHOD = method,
                      MSI = c(resLDA, resM, resR, resL, resJ, resP))
  RES_I <- data.frame(N = n, p = p, Sigma = "I", Alpha = alpha, Tau = tau, METHOD = method,
                      MSI = c(resLDA_I, resM_I, resR_I, resL_I, resJ_I, resP_I))
  
  return(rbind(RES_I, RES_S))
}

Test <- simu3momentsSIpar(12,1000, alphas=c(0.05,0.1,0.15),
                          taus=4:7, p=3)
# to reproduce simulation in paper chose approbriate settings. Is still slow!

                                 
#################################################################################################
# Code for the figure for a random guess MSI
#################################################################################################

# Function to simulate random guess MSI vlaues as a function of dimension p 
# If one estimates the optimal direction as a random guess in the unit sphere,
# then the distribution of the corresponding MSI is the same as for | e_1 ' U  | where U is uniformly random in the unit sphere.
# The distribution of (e_1' U)^2 is Beta(1/2, (p - 1)/2), so this lets us estimate the distribution of the random guess MSIs numerically.
# Estimates the (0.05, 0.50, 0.95)-quantiles of the random MSI based on "reps" repetitions
# input:
# -p: dimension of the problem
# -reps: number of repetitions for the simulation
# output:
# vecctor of length 3 containing the (0.05, 0.50, 0.95)-quantiles                                 
estim_random_guess_q <- function(p, reps){
  x <- sqrt(rbeta(reps, 1/2, (p - 1)/2))
  quantile(x, probs = c(0.05, 0.50, 0.95))
}

# To obtain the figure from the paper                                  
# Setting the dimensions p
p_set <- 2:100                        
# Setting the number of reps per dimension
reps <- 200000
# Initializing the result-object
res <- NULL
# Setting the seed for reproducibility                                 
set.seed(1234)
# Loop over dimensions and collect results
for(p in p_set){
  res <- rbind(res, c(p, estim_random_guess_q(p, reps)))
  print(p)
}
colnames(res) <- c("p", "q1", "q2", "q3")

# The plot
library(ggplot2)
ggplot(res, aes(x = p, y = q2)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), fill = "grey70") +
  geom_line() +
  labs(x = "Dimension p", y = "Distribution of MSI under random guess") +
  theme_bw()
                                 
                               

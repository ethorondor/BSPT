rm(list=ls())
options(scipen=999,digits = 4)
library(strucchange)
library(rlist)
library(tidyverse)
library(parallel)
library(timeSeries)
library(zoo)


nSim <- 100
getParamSpace <- function(param){
  mu1.space <- seq(param$mu1$prior[1],param$mu1$prior[2],by=param$mu1$step)
  return(mu1.space)
}
getEst.grid <- function(theObject){
  Y <- theObject$vec
  rho0 <- theObject$rho0
  p21 <- theObject$p21
  #param.space <- getParamSpace(theObject@param)
  pp <- theObject$param.space$ps.positive
  pn <- theObject$param.space$ps.negative
  for(i in 1:nrow(pp)){
    r <- getSumLoglk(pp[i])
    pp$sumLogL[i] <- r[1]
    pp$pi[i] <- r[2]
  }
  for(i in 1:nrow(pn)){
    r <- getSumLoglk(pn[i])
    pn$sumLogL[i] <- r[1]
    pn$pi[i] <- r[2]
  }
  if(sum(pp$SumLogL)>sum(pn$sumLogL)){
    ps <- pp
  }else{
    ps <- pn
  }
  ps$w <- exp(ps$sumLogL)/sum(exp(ps$sumLogL))
  rho.mean = sum(ps$rho*ps$w)
  mu1.mean = sum(ps$mu1*ps$w)
  pi <- getSumLoglk(rho.mean,mu1.mean,Y,rho0,p21)
  return(c(rho.mean, mu1.mean, pi[[2]]))
}
getLikelihood = function(mu1){
  #initialize result can improve performance
  l <- cbind(dnorm(Y,mean=0,sd=1),dnorm(Y,mean=mu1,sd=1))
  pi <- matrix(data=0,nrow=length(Y),ncol=2)
  f <- vector(mode = 'numeric',length=length(Y))
  # state    1     2
  # p = [  1-rho,  p21
  #      rho,1-p21]
  p <- matrix(c(1-rho,p21 ,rho,1-p21),nrow=2,ncol=2,byrow = TRUE)
  f[1] <- c(1,1)%*%(p%*%rho0*l[1,])
  pi[1,] <- (p%*%rho0*l[1,])/f[1]
  for(tmp.c in 2:length(Y)){
    f[tmp.c] <- c(1,1)%*%(p%*%pi[tmp.c-1,]*l[tmp.c,])
    pi[tmp.c,] <- (p%*%pi[tmp.c-1,]*l[tmp.c,])/f[tmp.c]
  }
  rm(l)
  return(list(f,pi))
}

getSumLoglk <- function(mu1){
  ll <- getLikelihood(mu1)
  return( c(sum(log(ll[[1]][2:length(Y)])),ll[[2]][length(Y),2] ))
}
result <- vector(mode = "numeric",length = nSim)
param.positive = list(mu1=list(prior = c(0.5,1.50), step = 0.01))
param.negative = list(mu1=list(prior = c(-1.50,-0.50), step = 0.01))
ps <- list(ps.positive=getParamSpace(param.positive),ps.negative = getParamSpace(param.negative))
df1 <- data.frame(y=rnorm(400))
df1[200:400,"y"] <- df1[200:400,"y"]+0.8
sr =  data.frame(rho=vector(mode = "numeric",length = length(df1$y)),
                       mu1=vector(mode = "numeric",length = length(df1$y)),
                       pi=vector(mode = "numeric",length = length(df1$y)),
                       piStar=vector(mode = "numeric",length = length(df1$y)))
bs <- list("Bsquid",
                  vec = df1$y[1:300],
                  nSim = 5000,
                  rho0 = c(0.5,0.5),
                  rho = 0.01,
                  p21 = 0.05,
                  param = list(mu1=list(prior = c(0.5,1.50), step = 0.1)),
                  param.space = ps,
                  cnFctr = 0.01,
                  mcCore = 10
)
Y <<- bs$vec
rho0 <<- bs$rho0
p21 <<- bs$p21
rho <<- bs$rho
getSumLoglk(0.5)
sp <- getParamSpace(bs$param)
mu <- vapply(sp,getSumLoglk,numeric(2))
df <- data.frame(p = sp,sumlikelihood=mu[1,],pi = mu[2,])
est <- getEst.grid(bs)
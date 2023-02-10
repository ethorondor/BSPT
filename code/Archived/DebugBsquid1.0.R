rm(list=ls())
options(scipen=999,digits = 4)
library(strucchange)
library(rlist)
library(tidyverse)
library(parallel)
library(timeSeries)
library(zoo)

#setwd("/mnt/MyDoc/Dropbox/Research/MonitoringStructureBreaks/code")

setClass(
  "Bsquid",
  slots = c(
    vec = "numeric",
    nSim = "numeric",
    rho0 = "numeric",
    pi21 = "numeric",
    param = "list" ,
    cnFctr = "numeric"
  )
)

setGeneric(name = "Bsquid.main",
           def = function(theObject)
           {
             standardGeneric("Bsquid.main")
           })
setMethod("Bsquid.main",
          signature = "Bsquid",
          definition = function(theObject)
          {
            Y <<- theObject@vec
            rho0 <<- theObject@rho0
            rho <<- theObject@rho
            p21 <<- theObject@p21
            param.space <- theObject@param.space
            est <- getEst.grid(param.space)
            #est <- getEst(theObject)
            piStar <- getPistar(est[1],theObject@cnFctr)
            return(c(est[1],est[2],piStar))
          }
)

############################################################################
############################################################################
############################################################################
df1 <- data.frame(y=rnorm(300))
df1[150:300,"y"] <- df1[150:300,"y"]+0.8
param = list(mu1.positive=list(prior=c(0.50,1.50),step=0.01),
             mu1.negative=list(prior=c(-1.50,-0.50),step=0.01))
getLikelihood = function(mu1=numeric(), others=list() ){
  Y = others$Y
  rho = others$rho
  rho0 = others$rho0
  p21 = others$p21
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
  return(c(sum(log(f)), pi[length(Y),2]))
}

  others <- list(
      Y = df1$y[1:200],
      rho = 0.01,
      rho0 = c(0.5,0.5),
      p21 = 0.05
  )

  pm <- c(seq(param$mu1.positive$prior[1],param$mu1.positive$prior[2],by=param$mu1.positive$step),
          seq(param$mu1.negative$prior[1],param$mu1.negative$prior[2],by=param$mu1.negative$step))
  pp <- t(sapply(pm,getLikelihood, others=others, simplify = TRUE,USE.NAMES = FALSE))
  p.df <- data.frame(mu1=pm,loglk=pp[,1],pi=pp[,2])
  p.df$w <- exp(p.df$loglk)/sum(exp(p.df$loglk))
  mu1 = sum(p.df$mu1*p.df$w)
  pi = sum(p.df$pi*p.df$w)
############################################################################
############################################################################
############################################################################

# the optimization function will take the parameter space as an input.
# the value of pi_prime will be calculated for each theta
# then aggregate the pi_prime across the parameter space with the given weight
# there are two important spaces in this function: the data space and the pi space
getPistar=function(param.space,c){
  ##specify data space by specify the upper and lower bound
  rho = 0.01
  lowBound <- -3
  upperBound <- 3
  inc = 0.01
  smplSp <- seq(from=lowBound, to=upperBound, by=inc)
  pi <- seq(from=0,to=1,by=0.01)
  pi <- pi[which(pi>0)]
getTranMatrix <- function(mu1=numeric()){
    piPrime <- matrix(data=0,nrow=length(pi),ncol=length(smplSp))
    Z1 <- dnorm(smplSp,mean=0,sd=1)
    Z2 <- dnorm(smplSp,mean=mu1,sd=1)
    P1 <- Z1/sum(Z1)
    P2 <- Z2/sum(Z2)
    lr <- Z2/Z1
    #piPrime is a matrix with dimension: size of pi by size of data space
    #the matrix is defined as given pi what is pi_prime at each data point
    T.tmp <- matrix(0,nrow=length(pi),ncol = length(pi))
    for(j in 1:length(pi)){
      like <- round(100*(lr*(pi[j]+rho*(1-pi[j])))/(lr*(pi[j]+rho*(1-pi[j]))+(1-rho)*(1-pi[j])))
      like[which(like=='inf')] <- 1
      piPrime.like <- data.frame(l = like, p1=P1, p2=P2)
      pp <- (piPrime.like%>%group_by(l)%>%summarize(pp1=sum(p1),pp2=sum(p2)))
      pp$P <- (1-pi[j])*(1-rho)*pp$pp1+(pi[j]+(1-pi[j])*rho)*pp$pp2
      pp <- pp[which(pp$l>0),]
      pp$P <- pp$P/sum(pp$P)
      for(i in 1:nrow(pp)){
        T.tmp[j,pp$l[i]] = pp$P[i]
      }
    }
    return(T.tmp)
}
T <- matrix(0,nrow=length(pi),ncol = length(pi))
for(i in 1:nrow(param.space)){
  T <- getTranMatrix(mu1=param.space$mu1[i])*param.space$w[i]+T
}
  h <- c*T%*%pi+1+(c-1)*pi
  g <- as.matrix(1-pi)
  r <- cbind(g,h)
  Q <- pmin(g,h)
  ep = 1
  cnt = 0
  while(ep > 0.001 & cnt <500){
    Q1 <- Q
    cnt = cnt+1
    h <- T%*%Q+c*pi
    Q <- pmin(g,h)
    ep <- max(abs(Q1-Q),0,na.rm=TRUE)
    diff <- (g-h)^2
    if(cnt < 499){
      piStar <- pi[which.min(diff)]
    }else{
      piStar <- NA
    }
  }
  return(piStar)
}
piStar <- getPistar(param.space <- p.df[,c(1,4)], c = 0.01)


#######################################################################
#######################################################################
#######################################################################
df1 <- data.frame(y=rnorm(300))
df1[150:300,"y"] <- df1[150:300,"y"]+0.8
result = data.frame(mu1=numeric(300),pi=numeric(300))
for(i in 100:300){
  others <- list(
    Y = df1$y[1:i],
    rho = 0.03,
    rho0 = c(0.5,0.5),
    p21 = 0.00
  )
  
  pm <- c(ps$mu1.pos,ps$mu1.neg)
  pp <- t(sapply(pm,getLikelihood, others=others, simplify = TRUE,USE.NAMES = FALSE))
  p.df <- data.frame(mu1=pm,loglk=pp[,1],pi=pp[,2])
  p.df$w <- exp(p.df$loglk)/sum(exp(p.df$loglk))
  result$mu1[i] = sum(p.df$mu1*p.df$w)
  result$pi[i] = sum(p.df$pi*p.df$w)
}
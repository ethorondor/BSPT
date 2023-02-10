rm(list=ls())
#options(scipen=999,digits = 4)
library(strucchange)
library(rlist)
library(tidyverse)
library(parallel)
library(timeSeries)
library(zoo)

##############################################################################
##############################################################################
##############################################################################

Bsquid.main = function(theObject)
{
  others <- list(
    Y = theObject$vec,
    rho0 = theObject$rho0,              
    rho = theObject$rho,
    p21 = theObject$p21
  )
  pm <- c(seq(theObject$param$mu1.positive$prior[1],theObject$param$mu1.positive$prior[2],by=theObject$param$mu1.positive$step),
          seq(theObject$param$mu1.negative$prior[1],theObject$param$mu1.negative$prior[2],by=theObject$param$mu1.negative$step))
  pp <- t(sapply(pm,getLikelihood, others=others, simplify = TRUE,USE.NAMES = FALSE))
  p.df <- data.frame(mu1=pm,loglk=pp[,1],pi=pp[,2])
  p.df$w <- exp(p.df$loglk-mean(p.df$loglk))/(sum(exp(p.df$loglk-mean(p.df$loglk))))
  #mu1 = sum(p.df$mu1*p.df$w)
  #Pi = sum(p.df$pi*p.df$w)
  #piStar <- getPistar(param.space <- p.df[,c(1,4)], c=theObject@cnFctr, rho=theObject@rho)
  pistar <- getPistar(param.space <- p.df[,c(1,4)], c=theObject$cnFctr, rho=theObject$rho)
  return(c(sum(p.df$mu1*p.df$w), sum(p.df$pi*p.df$w), pistar)
  )
}


############################################################################
############################################################################
############################################################################
getLikelihood = function(mu1, others ){
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


getTranVector <- function(pi,T.lst){
  rho <- T.lst$rho
  t.tmp <- T.lst$t.tmp 
  df <- T.lst$df
  like <- round(100*(df$lr*(pi+rho*(1-pi)))/(df$lr*(pi+rho*(1-pi))+(1-rho)*(1-pi)))
  like[which(like=='inf')] <- 1
  df$l <- like
  pp <- (df%>%group_by(l)%>%summarize(pp1=sum(P1),pp2=sum(P2)))
  pp$P <- (1-pi)*(1-rho)*pp$pp1+(pi+(1-pi)*rho)*pp$pp2
  pp <- pp[which(pp$l>0),]
  pp$P <- pp$P/sum(pp$P)
  t.tmp[pp$l] = pp$P
  return(t.tmp)
}

getTranMatrix <- function(arg1,arg2,arg3,arg4){
  mu1 = arg1
  rho = arg2
  smplSp = arg3
  pi.v = arg4
  T.lst <- list(df=data.frame(
    Z1 = dnorm(smplSp,mean=0,sd=1),
    Z2 = dnorm(smplSp,mean=mu1,sd=1),
    P1 = dnorm(smplSp,mean=0,sd=1)/sum(dnorm(smplSp,mean=0,sd=1)),
    P2 = dnorm(smplSp,mean=mu1,sd=1)/sum(dnorm(smplSp,mean=mu1,sd=1)),
    lr = dnorm(smplSp,mean=mu1,sd=1)/dnorm(smplSp,mean=0,sd=1)  ),
    t.tmp = vector(mode='numeric',length=length(pi.v)),
    rho = rho
  )
  return(sapply(pi.v,getTranVector,T.lst)) 
  #piPrime is a matrix with dimension: size of pi by size of data space
  #the matrix is defined as given pi what is pi_prime at each data point
}

############################################################################
############################################################################
############################################################################

# the optimization function will take the parameter space as an input.
# the value of pi_prime will be calculated for each theta
# then aggregate the pi_prime across the parameter space with the given weight
# there are two important spaces in this function: the data space and the pi space
getPistar=function(param.space,c,rho){
  ##specify data space by specify the upper and lower bound
  #param.space <- p.df[,c(1,4)] 
  #c = 0.01
  rho <- rho
  lowBound <- -3
  upperBound <- 3
  inc = 0.01
  smplSp <- seq(from=lowBound, to=upperBound, by=inc)
  P <- seq(from=0,to=1,by=0.01)
  pi.v <- P[which(P>0)]
  # get weighted average of transition matrix across parameter space, weighted with posterior parameter distribution
  mx.pm <- lapply(param.space$mu1, getTranMatrix,arg2=rho,arg3=smplSp,arg4=pi.v)
  T <- Reduce("+",Map("*", mx.pm, param.space$w))
  
  h <- c*T%*%pi.v+1+(c-1)*pi.v
  g <- as.matrix(1-pi.v)
  r <- cbind(g,h)
  Q <- pmin(g,h)
  ep = 1
  cnt = 0
  while(ep > 0.001 & cnt <500){
    Q1 <- Q
    cnt = cnt+1
    h <- T%*%Q+c*pi.v
    Q <- pmin(g,h)
    ep <- max(abs(Q1-Q),0,na.rm=TRUE)
    diff <- (g-h)^2
    if(cnt < 499){
      piStar <- pi.v[which.min(diff)]
    }else{
      piStar <- NA
    }
  }
  return(piStar)
}
##############################################################################
####################### Single Run ###########################################
##############################################################################
nSim <- 1
simLength <<- 600
simStart <<- 160
brkTm <<-150
df1 <- data.frame(y=rnorm(simLength,mean=0,sd=1))
df1[brkTm:simLength,"y"] <- df1[brkTm:simLength,"y"]+0.8
bs <- list(vec = df1$y,
           rho0 = c(0.5,0.5),
           rho = 0.01,
           pi21 = 0,
           param = list(mu1.positive=list(prior=c(0.50,1.50),step=0.01),
                        mu1.negative=list(prior=c(-1.50,-0.50),step=0.01)),
           cnFctr = 0.01
)
rs <- Bsquid.main(bs)
############################################################################
############################################################################
############################################################################

#nSim <- 1
#simLength <<- 600
#simStart <<- 160
#brkTm <<-150
#result <- vector(mode = "numeric",length = nSim)

#  s <- data.frame(mu1=numeric(),pi=numeric(),piStar=numeric(),count=numeric())
#  df1 <- data.frame(y=rnorm(simLength,mean=0,sd=1))
#  df1[brkTm:simLength,"y"] <- df1[brkTm:simLength,"y"]+0.8
#  bs <- list()
#  for(i in simStart+1:length(df1$y)){
#    bs <- list(vec = df1$y[1:i],
#              rho0 = c(0.5,0.5),
#              rho = 0.01,
#              pi21 = 0,
#              param = list(mu1.positive=list(prior=c(0.50,1.50),step=0.01),
#                           mu1.negative=list(prior=c(-1.50,-0.50),step=0.01)),
#              cnFctr = 0.005
#    )
#    s[i-simStart,] <- c(Bsquid.main(bs),i)
#    if(s[i-simStart,2]>s[i-simStart,3]){
#      break
#    }
#  }
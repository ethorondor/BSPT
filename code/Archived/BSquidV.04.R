library(tidyverse)
library(parallel)

#library(zoo)
setClass(
  "Bsquid",
  slots = c(
    vec = "numeric",
    nSim = "numeric",
    rho0 = "numeric",
    p21 = "numeric",
    param = "list" ,
    param.space = "list",
    cnFctr = "numeric",
    mcCore = "numeric"
  )
)
############################################################################
############################################################################
############################################################################
# df1 <- data.frame(y=rnorm(300))
# df1[150:300,"y"] <- df1[150:300,"y"]+0.8
# simu <- new("Bsquid", vec = df1$y,
#             nSim = 10000,
#             rho0 = c(0.5,0.5),
#             p21 = 0.05,
#             param = list(rho=list(prior = c(0.00,1.00), step = 0.01),
#                          mu1=list(prior = c(0.5,1.50), step = 0.10)
#             ),
#             cnFctr = 0.01,
#             mcCore = 14
# )

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
            p21 <<- theObject@p21
            param.space <- theObject@param.space
            est <- getEst.grid(param.space)
            #est <- getEst(theObject)
            piStar <- getPistar(est[1],est[2],theObject@cnFctr)
            return(c(est[1],est[2],est[3],piStar))
          }
)

getEst.grid <- function(param.space){
  pp <- param.space$ps.positive
  pn <- param.space$ps.negative
  for(i in 1:nrow(pp)){
    r <- getSumLoglk(pp$rho[i],pp$mu1[i])
    pp$sumLogL[i] <- r[1]
    pp$pi[i] <- r[2]
  }
  for(i in 1:nrow(pn)){
    r <- getSumLoglk(pn$rho[i],pn$mu1[i])
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

getParamSpace <- function(param){
  rho.space <- list(seq(param$rho$prior[1],param$rho$prior[2],by=param$rho$step))
  mu1.space <- list(seq(param$mu1$prior[1],param$mu1$prior[2],by=param$mu1$step))
  s <- length(rho.space[[1]])*length(mu1.space[[1]])
  r = data.frame(rho=vector(mode = 'numeric',length=s),
                 mu1=vector(mode = 'numeric',length = s),
                 sumLogL=vector(mode='numeric',length=s),
                 pi = vector(mode='numeric',length = s))
  c=1
  for(j in rho.space[[1]]){
    for(k in mu1.space[[1]]){
      r$rho[c] = j
      r$mu1[c] = k
      c=c+1
    }
  }
  rm(rho.space,mu1.space,s)
  return(r)
}

getLikelihood = function(rho, mu1){
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

getSumLoglk <- function(rho,mu1){
  ll <- getLikelihood(rho,mu1)
  return( c(sum(log(ll[[1]][2:length(Y)])),ll[[2]][length(Y),2] ))
}

##there are two important spaces in this function: the data space and the pi space
getPistar=function(rho,mu1,c){
  #  #specify data space by specify the upper and lower bound
  lowBound <- -3
  upperBound <- 3
  inc = 0.01
  smplSp <- seq(from=lowBound, to=upperBound, by=inc)
  pi <- seq(from=0,to=1,by=0.01)
  pi<-pi[which(pi>0)]
  #  #initialize vector and matrix
  piPrime <- matrix(data=0,nrow=length(pi),ncol=length(smplSp))
  Z1 <- dnorm(smplSp,mean=0,sd=1)
  Z2 <- dnorm(smplSp,mean=mu1,sd=1)
  
  P1 <- Z1/sum(Z1)
  P2 <- Z2/sum(Z2)
  
  lr <- Z2/Z1
  #piPrime is a matrix with dimension: size of pi by size of data space
  #the matrix is defined as given pi what is pi_prime at each data point
  for(j in 1:length(pi)){
    like <- (lr*(pi[j]+rho*(1-pi[j])))/(lr*(pi[j]+rho*(1-pi[j]))+(1-rho)*(1-pi[j]))
    like[which(like=='inf')]<- 1
    piPrime[j,]=t(round(100*like))
  }
  #initialize matrix
  Pp1 <- matrix(0,nrow=length(smplSp),ncol=length(pi))
  Pp2 <- Pp1
  T1 <- matrix(0,nrow=length(pi),ncol = length(pi))
  T2 <- T1
  N_T <- T1
  for(j in 1:length(pi)){
    Pp1 <- matrix(0,nrow=length(smplSp),ncol=length(pi))
    Pp2 <- matrix(0,nrow=length(smplSp),ncol=length(pi))
    for(i in 1:length(smplSp)){
      #Pp1 is data space by pi space, it is defined as give pi, the probability of getting the pi_prime given state 1
      Pp1[i,piPrime[j,i]]=P1[i]
      #Pp1 is data space by pi space, it is defined as give pi, the probability of getting the pi_prime given state 2
      Pp2[i,piPrime[j,i]]=P2[i]
    }
    T1[j,]<-colSums(Pp1,na.rm = FALSE,dims = 1)
    T2[j,]<-colSums(Pp2,na.rm = FALSE,dims = 1)
  }
  
  for(j in 1:length(pi)){
    N_T[j,] <- pi[j]*T2[j,]+(1-(pi[j]))*(1-rho)*T1[j,]+(1-pi[j])*rho*T2[j]
  }
  
  T_tmp <- rep(1,length(pi))
  T_tmp[rowSums(N_T)>1] <- rowSums(N_T)[rowSums(N_T)>1]
  T <- N_T/T_tmp
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

#rs <- Bsquid.main(simu)



getEst = function(theObject){
  pm.rho <- vector(mode = "numeric",length = theObject@nSim)
  pm.mu1 <- vector(mode="numeric",length = theObject@nSim)
  param.rho <- runif(theObject@nSim,min=theObject@param[[1]]$prior[1],max=theObject@param[[1]]$prior[2])
  param.mu1 <-  runif(theObject@nSim,min=theObject@param[[2]]$prior[1],max=theObject@param[[2]]$prior[2])
  for(i in 1:theObject@nSim){
    likelihood <- getLikelihood(param.rho[i],param.mu1[i],theObject@vec,theObject@rho0,theObject@p21)
    sumLogL <- sum(log(likelihood[1][[1]][2:length(theObject@vec)]))
    if(sumLogL == '-Inf' | sumLogL == 'NaN'){sumLogL = -99999}
    if(i==1){
      pSwch = 1
    }else{
      pSwch = exp(sumLogL - sumLogL_0)
    }
    if(runif(1)<pSwch){
      tem.rho = param.rho[i]
      tem.mu1 = param.mu1[i]
      sumLogL_0 = sumLogL
    }
    pm.rho[i] <- tem.rho
    pm.mu1[i] <- tem.mu1
  }
  rho.mean <- mean(pm.rho)
  mu1.mean <- mean(pm.mu1)
  estL <- getLikelihood(rho.mean,mu1.mean,theObject@vec,theObject@rho0,theObject@p21)
  pi <- estL[[2]][length(theObject@vec),2]
  return(list(rho.mean=rho.mean,mu1.mean=mu1.mean,pi=pi))
}


############################################################################
############################################################################
############################################################################
#getResult <- function(result){
#  brk <- min(result$opt$cnt[result$opt$pi>result$opt$pi.star],na.rm = TRUE)
#  return(brk)
#}

#simu.result <- mclapply(simu,BSquid.main,mc.cores = 10L)
#simuResult <- lapply(simu.result,getResult)
#df <- data.frame(matrix(unlist(simuResult), nrow=length(simuResult), byrow=T))
#colnames(df) <- "StrBrk"
#simu.result <- list(strbrk = df, 
#                    specification="t=1, cnFctr=0.01,nSim=5000,rho0=0.5")
#save(simu.result,file="simu.result_1_01.rdata")

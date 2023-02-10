library(tidyverse)
library(parallel)

#library(zoo)
setClass(
  "BSquid",
  slots = c(
    ts = "list",
    nSim = "numeric",
    rho0 = "numeric",
    pi12 = "numeric",
    model = "list" ,
    cnFctr = "numeric",
    mcCore = "numeric"
  )
)

setGeneric(name = "BSquid.main",
           def = function(theObject)
           {
             standardGeneric("BSquid.main")
           })

setMethod("BSquid.main",
          signature = "BSquid",
          definition = function(theObject)
          {
            Y = theObject@ts$y
            pi12 = theObject@pi12
            nsim = theObject@nSim
            rho0  = theObject@rho0
            cnFctr = theObject@cnFctr
            model = theObject@model
            OT.output <- matrix(data=NA,nrow = length(Y),ncol=7)
            #create a list of input vector
            lst.Y <- list()
            for(temCount in 1:(length(Y)-100)){
              lst.Y[[temCount]] <- Y[1:100+temCount]
            }
            lst.est <- function(Y){
              paramEst(model,Y,pi12,nsim,rho0)
            }
            output <- mclapply(lst.Y,lst.est, mc.cores = getOption("mc.cores",theObject@mcCore))
            OT <- function(estModel){
              BSquid.OT(estModel,c=cnFctr)
            }
            OT.opt <- mclapply(output,OT, mc.cores = getOption("mc.cores",theObject@mcCore))
            opt <- data.frame(cnt=seq(1:(length(Y)-100)),
                              p11=1:(length(Y)-100),
                              mu1=1:(length(Y)-100),
                              pi = 1:(length(Y)-100),
                              pi.star=1:(length(Y)-100),
                              y=1:(length(Y)-100))
            for(i in 1:(length(Y)-100)){
              opt$cnt[i] <- i+100
              opt$p11[i] <- output[[i]]$param$p11$value
              opt$mu1[i]  <- output[[i]]$param$mu1$value
              opt$pi[i] <- output[[i]]$output[3]
              opt$pi.star[i] <- OT.opt[[i]]
              opt$y[i] <- Y[i+100]
            }
            OPT <- list(opt=opt,y=Y)
            return(OPT)
          }
)


ranGenParam <- function(paramLst){
  paramLst$value <- round(runif(1,min=paramLst$prior[1],max=paramLst$prior[2]),2)
  return(paramLst)
}
paramEst = function(model,Y,pi12,nsim,rho0){
  pm <- matrix(data=0,nrow=nsim,ncol=length(model$param))
  temParam <- matrix(data=0,nrow=1,ncol=length(model$param))
  ranParam <- cbind(runif(nsim,min=model$param[[1]]$prior[1],max=model$param[[1]]$prior[2]),
                    runif(nsim,min=model$param[[2]]$prior[1],max=model$param[[2]]$prior[2])
                    )
  for(i in 1:nsim){
    for (j in 1:length(model$param)) {
      model$param[[j]]$value = ranParam[i,j]
    }
    likelihood <- mcLikelihood(model,Y,rho0,pi12)
    sumLogL <- sum(log(likelihood[1][[1]][2:length(Y)]))
            if(sumLogL == '-Inf' | sumLogL == 'NaN'){sumLogL = -99999}
            if(i==1){
                        pSwch = 1
            }else{
                        LLR = exp(sumLogL - sumLogL_0)
                        pSwch <- min(1,LLR)
            }
            if(runif(1)<pSwch){
                for (j in 1:length(model$param)) {
                    temParam[1,j] = model$param[[j]]$value
                }
              sumLogL_0 = sumLogL
            }
     pm[i,] <- temParam
  }
  
  pm <- as.data.frame(pm)
  colnames(pm) <- model$paramNames
  pmSpaceBar <- colMeans(pm)
  for (j in 1:length(model$param)) {
    model$param[[j]]$value = pmSpaceBar[j]
    model$output[j] = pmSpaceBar[j]
  }
  estL <- mcLikelihood(model,Y,rho0,pi12)

  model$output[length(model$output)] <- estL[[2]][length(Y),2]
  return(model)
}

mcLikelihood = function(model, Y, rho0, p12){
  N <- matrix(data=0,nrow = length(Y),ncol = 2)
  pi <- matrix(data=0,nrow=length(Y),ncol=2)
  f <- vector(mode = 'numeric',length=length(Y))
  p <- matrix(c(model$param[[1]]$value,p12 ,1-model$param[[1]]$value,1-p12),nrow=2,ncol=2,byrow = TRUE)
  for(tmp.c in 1:length(Y)){
    if(tmp.c == 1){
      N[tmp.c,] <- c(returnErrorDnorm(model,Y[tmp.c],1),returnErrorDnorm(model,Y[tmp.c],2))
      f[tmp.c] <- c(1,1)%*%(p%*%rho0*N[tmp.c,])
      pi[tmp.c,] <- (p%*%rho0*N[tmp.c,])/f[tmp.c]
    }else{
      N[tmp.c,] <- c(returnErrorDnorm(model,Y[tmp.c],1),returnErrorDnorm(model,Y[tmp.c],2))
      f[tmp.c] <- c(1,1)%*%(p%*%pi[tmp.c-1,]*N[tmp.c,])
      pi[tmp.c,] <- (p%*%pi[tmp.c-1,]*N[tmp.c,])/f[tmp.c]
    }
  }
  return(list(f,pi))
}

returnErrorDnorm = function(model,y,state){
  if(state==1){
    fError <- dnorm(y,mean=0,sd=1)
  }else{
    fError <- dnorm(y,mean=model$param[[2]]$value,sd=1)
  }
  return(fError)
}
##there are two important spaces in this function: the data space and the pi space
BSquid.OT=function(model,c){
    rho <- 1-model$param[[1]]$value
#  #specify data space by specify the upper and lower bound
    #sim1 <- eval(parse(text=model$model.null))
    #sim2 <- eval(parse(text=model$model.alt))
    #sim <- c(sim1,sim2)
    lowBound <- -3
    upperBound <- 3
    inc = 0.01
    Y <- seq(from=lowBound, to=upperBound, by=inc)
    pi <- seq(from=0,to=1,by=0.01)
    pi<-pi[which(pi>0)]
#  #initialize vector and matrix
    Z1 <- vector(mode='numeric',length=length(Y))
    Z2<-Z1
    piPrime <- matrix(data=0,nrow=length(pi),ncol=length(Y))

    for(n in 1: length(Y)){
      Z1[n]<-returnErrorDnorm(model,Y[n],1)
      Z2[n]<-returnErrorDnorm(model,Y[n],2)
    }

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
  Pp1 <- matrix(0,nrow=length(Y),ncol=length(pi))
  Pp2 <- Pp1
  T1 <- matrix(0,nrow=length(pi),ncol = length(pi))
  T2 <- T1
  N_T <- T1
  for(j in 1:length(pi)){
    Pp1 <- matrix(0,nrow=length(Y),ncol=length(pi))
    Pp2 <- matrix(0,nrow=length(Y),ncol=length(pi))
    for(i in 1:length(Y)){
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
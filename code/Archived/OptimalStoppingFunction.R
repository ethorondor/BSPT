rm(list = ls())
cnfctor <- list(0.01, 0.05, 0.07, 0.1)
pi11 <- list(seq(0,1,by=0.01))
mu1 <- list(seq(0.5,1.5,by=0.1))


getOTInput <- function(cnfctor,pi11,mu1){
  s <- length(cnfctor)*length(pi11[[1]])*length(mu1[[1]])
  r = data.frame(cnfctor=1:s,pi11=1:s,mu1=1:s)
  c=1
  for(i in cnfctor){
    for(j in pi11[[1]]){
      for(k in mu1[[1]]){
            r$cnfctor[c] <- i
            r$pi11[c] <- j
            r$mu1[c] <- k
            c=c+1
      }
    }
  }  
  return(r)
}

getPistar=function(est,c){
  rho <- 1-est$pi11.mean
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
  Z2 <- dnorm(smplSp,mean=est$mu1.mean,sd=1)
  
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
OTInput <- getOTInput(cnfctor,pi11,mu1)
est <- list(pi11.mean=OTInput$pi11[1],mu1.mean=OTInput$mu1[1])
c <- OTInput$cnfctor[1]
est <- list(pi11.mean=0.97,mu1.mean=0.8)
c <- 0.01
pistar <- getPistar(est,c)
system.time(getPistar(est,c))
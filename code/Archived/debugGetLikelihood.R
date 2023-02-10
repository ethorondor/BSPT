getLikelihood = function(rho, mu1 , Y, rho0, p21){
  #initialize result can improve performance
  l <- cbind(dnorm(Y,mean=0,sd=1),dnorm(Y,mean=mu1,sd=1))
  pi <- matrix(data=0,nrow=length(Y),ncol=2)
  f <- vector(mode = 'numeric',length=length(Y))
  # state    1     2
  # p = [  p11,  p21
  #      1-p11,1-p21]
  p <- matrix(c(1-rho,p21 ,rho,1-p21),nrow=2,ncol=2,byrow = TRUE)
  f[1] <- c(1,1)%*%(p%*%rho0*l[1,])
  pi[1,] <- (p%*%rho0*l[1,])/f[1]
  for(tmp.c in 2:length(Y)){
    f[tmp.c] <- c(1,1)%*%(p%*%pi[tmp.c-1,]*l[tmp.c,])
    pi[tmp.c,] <- (p%*%pi[tmp.c-1,]*l[tmp.c,])/f[tmp.c]
  }
  return(list(f,pi))
}

getSumLoglk <- function(rho,mu1){
  ll <- getLikelihood(rho,mu1,Y,rho0,p21)
  return( c(sum(log(ll[[1]][2:length(Y)])),ll[[2]][length(Y),2] ))
}

df1 <- data.frame(y=rnorm(300,mean=0,sd=1))
df1[100:300,"y"] <- df1[100:300,"y"]+0.8

lk <- getLikelihood(0.01,0.8,df1$y,c(0.5,0.5),0.05)
sum(log(lk[[1]][1:300]))
View(lk[[2]])


rm(list=ls())
options(scipen=999,digits = 4)
library(strucchange)
library(rlist)
library(tidyverse)
library(parallel)
library(timeSeries)
library(zoo)

setwd("/mnt/MyDoc/Dropbox/Research/MntStrBrk/code")

source("Bsquid.1.0.R")
nSim <- 1000
simLength <<- 600
simStart <<- 100
brkTm <<-110
result <- vector(mode = "numeric",length = nSim)

monStrBrk <- function(c){
  s <- data.frame(mu1=numeric(),pi=numeric(),piStar=numeric(),count=numeric())
  df1 <- data.frame(y=rnorm(simLength,mean=0,sd=1))
  df1[brkTm:simLength,"y"] <- df1[brkTm:simLength,"y"]+0.8
  bs <- list()
  for(i in simStart+1:length(df1$y)){
    bs <- new("Bsquid",
              vec = df1$y[1:i],
              rho0 = c(0.5,0.5),
              rho = 0.01,
              p21 = 0,
              param = list(mu1.positive=list(prior=c(0.50,1.50),step=0.01),
                           mu1.negative=list(prior=c(-1.50,-0.50),step=0.01)),
              cnFctr = 0.01
              )
    s[i-simStart,] <- c(Bsquid.main(bs),i)
    if(s[i-simStart,2]>s[i-simStart,3]){
      break
    }
  }
  if(any(s$pi>s$piStar)){
    return(min(s$count[s$pi>s$piStar],na.rm = TRUE))
  }else{
    return('NA')
  }
}
res = list()
for(i in 1:nSim){
  res[[i]] = monStrBrk(i)
}
#res <- mclapply(seq(1:nSim),monStrBrk,mc.cores = 12)
rs <- vector(mode = "numeric",length=length(res))
for(i in 1:length(res)){
  rs[i] = res[[i]]
}  

rs <- list(rs, des ="110/600 c=0.01,p21=0.00")

save(rs,file="../output/simOP110.01.01.Rdata")
#singal run
#df1 <- data.frame(y=rnorm(400))
#df1[200:400,"y"] <- df1[200:400,"y"]+0.8
#bs <- new("Bsquid",
#          vec = df1$y,
#          rho0 = c(0.5,0.5),
#          rho = 0.01,
#          pi21 = 0,
#          param = list(mu1.positive=list(prior=c(0.50,1.50),step=0.01),
#                       mu1.negative=list(prior=c(-1.50,-0.50),step=0.01)),
#          cnFctr = 0.01
#)

#sr <- Bsquid.main(bs)
#if(sr[2]>sr[3]){
#  break
#}
#system.time(Bsquid.main(bs))  


#system.time(mclapply(s,Bsquid.main,mc.cores = 10))

rm(list=ls())
options(scipen = 999,digits = 4)
require(tidyverse)
require(parallel)
require(timeSeries)
require(zoo)
require(strucchange)
# simulation with OLS-CUSUM
n = 10000
opt <- vector(mode = 'numeric',length=n)
for(i in 1:n){
    df <- data.frame(y=rnorm(1000))
    df[150:300,"y"] <- df[150:300,"y"]+0.8
    me1 <- mefp(y~1,data=df[1:100,,drop=FALSE],type="OLS-CUSUM",h=1,alpha = 0.07)
    me2 <- monitor(me1,data=df)
    opt[i] = me2$breakpoint
}
na = opt[is.na(opt)]
nna = na.omit(opt)
fa = nna[nna<150]
delay = nna[nna>=150]
mean(delay,na.rm=TRUE)
df <- load('/mnt/MyDoc/Dropbox/Research/MntStrBrk/data/GermanM1.rda')

rm(list=ls())
options(scipen=999,digits = 4)
library(strucchange)
library(rlist)
library(tidyverse)
library(parallel)
library(timeSeries)
library(zoo)

setwd("/mnt/MyDoc/Dropbox/Research/MonitoringStructureBreaks/code")

tst <- list(
  nsim = 100,
  n = 100,
  T = 3,
  t = 1
)

#####################################################################################
################################ Bsquid #############################################
#####################################################################################
source("BSquidV.03.R")
# create a list of Bsquid objects
#simuResult <- data.frame(c1=0,c2=0,c3=0)
smpl <- tst$T*tst$n
brk <- tst$t*tst$n+1
simu <- list()
for(j in 1:10){
  df1 <- data.frame(y=rnorm(smpl))
  df1[brk:smpl,"y"] <- df1[brk:smpl,"y"]+0.8
  simu[[j]] <- new("BSquid",
                        ts = df1,
                        nSim = 1000,
                        rho0 = c(0.5,0.5),
                        pi12 = 0.05,
                        model=list( paramNames = c("p11","mu"),
                                    param = list(p11=list(prior = c(0.00,1.00), value = 0.50),
                                                 mu1=list(prior = c(0.5,1.50), value = 0.00)),
                                    output = c(0,1,2)
                        ),
                        cnFctr = 0.01,
                        mcCore = 14
)
}

getResult <- function(result){
  brk <- min(result$opt$cnt[result$opt$pi>result$opt$pi.star],na.rm = TRUE)
  return(brk)
}

simu.result <- mclapply(simu,BSquid.main,mc.cores = 10L)
simuResult <- lapply(simu.result,getResult)
df <- data.frame(matrix(unlist(simuResult), nrow=length(simuResult), byrow=T))
colnames(df) <- "StrBrk"
simu.result <- list(strbrk = df, 
                    specification="t=1, cnFctr=0.01,nSim=5000,rho0=0.5")
save(simu.result,file="simu.result_1_01.rdata")

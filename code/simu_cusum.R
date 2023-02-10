library(strucchange)
library(rlist)
library(parallel)
nsim <- 10000
simParam = list(c(110,0.8,0),
                c(150,0.8,0),
                c(200,0.8,0),
                c(300,0.8,0),
                c(110,0.3,0.5),
                c(150,0.3,0.5),
                c(200,0.3,0.5),
                c(300,0.3,0.5),
                c(110,0.1,0.7),
                c(150,0.1,0.7),
                c(200,0.1,0.7),
                c(300,0.1,0.7),
                c(110,0.0,0.8),
                c(150,0.0,0.8),
                c(200,0.0,0.8),
                c(300,0.0,0.8))

mb <- function(tau,a,g){
  df1 <- data.frame(y=rnorm(tau))
  for(i in tau:900){
    df1[i+1,'y'] = a + g*df1[i,'y'] + rnorm(1)
  }
  me1 <- mefp(y~1,data=df1[1:100,,drop=FALSE],type="OLS-CUSUM",alpha=0.05,h=0.05,dynamic=FALSE,rescale = FALSE)
  me2 <- monitor(me1, data=df1)
  return(me2$breakpoint)
}

crtCUSUM <- function(pm){
  tau = pm[1]
  a = pm[2]
  g = pm[3]
  bt = data.frame(bt = 0)
  for (n in 1:nsim){
    bt[n,"bt"] <- mb(tau,a,g)
  }
  df <- bt[!is.na(bt)]
  typeii <- (nsim-length(df))/nsim
  typei <- length(df[df<tau])/nsim
  ed <- mean(df[df>tau])-tau
  ed.sd <- sd(df[df>tau])
  
  op = list("tau"=tau,
            "a" = a,
            "g" = g,
            "expected delay" = ed,
            "sd" = ed.sd,
            "typeI" = typei,
            "typeII" = typeii)
  return(op)
}

l <- mclapply(simParam,crtCUSUM,mc.cores=10L)
list.save(l,"/home/elrond/Dropbox/Research/MntStrBrk/output/simuCUSUM.rdata")

load("/home/elrond/Dropbox/Research/MntStrBrk/output/simuCUSUM.rdata")

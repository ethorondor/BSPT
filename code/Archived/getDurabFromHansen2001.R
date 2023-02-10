#########################################################################
##  BREAKS.R
##  
##  Bruce E. Hansen
##  Department of Economics
##  Social Science Building
##  University of Wisconsin
##  Madison, WI 53706-1393
##  behansen@wisc.edu
##  http://www.ssc.wisc.edu/~bhansen/
##  
##  Replicates empirical work from
##  "The New Econometrics of Structural Change:
##  Dating Breaks in U.S. Labor Productivity" 
##  by Bruce E. Hansen
##  
#########################################################################
rm(list=ls())
setwd('/mnt/MyDoc/Dropbox/Research/MntStrBrk/reference/hansen2001')
seas <- 12		# seasonal frequency of data	
k <- 1		# AR order 		
trim <- .03		# trim % for estimation 	
trimsup <- .05	# trim % for sup test		
trimexp <- .05	# trim % for exp test		
bootw <- 0		# # of fixed-regressor bootstrap replications (set to 0 if not desired) 
conf <- .9		# Confidence Level 		

tit <- "Labor Productivity - "   # string: title #
indx <- 1 	      # variable index - 1 to 11 #
dat1 <- read.table("ip.dat")
dat2 <- read.table("workers.dat")
dat3 <- read.table("hours.dat")
dat <- log(dat1)-log(dat2)-log(dat3)
dat <- dat[,indx]
dat <- as.matrix(dat[(1-is.nan(dat))>0])
dat <- as.matrix((dat[2:nrow(dat)]-dat[1:(nrow(dat)-1)])*seas*100)
year_1 <- 1947+1/12
year_2 <- year_1-1/12+nrow(dat)/seas

sample_1 <- year_1   # first year to use in sample : should be "year_1" or larger #
sample_2 <- year_2

tr1 <- seas*(sample_1-year_1)
tr2 <- seas*(year_2-sample_2)
dy <- as.matrix(dat[(1+tr1):(nrow(dat)-tr2)])
#dy is the annulized growth rate in percentage form, the first month is 1947 Apr
write.csv(dy,file='hansen_durab.csv')
n <- nrow(dy)
year1 <- as.matrix(seq(sample_1+1/seas,sample_1+n/seas,1/seas)) # "date" vector, should correspond to data in "dy" #

tit2 <- rbind(
"Durables",
"SIC24",
"SIC25",
"SIC32",
"SIC33",
"SIC34",
"SIC35",
"SIC36",
"SIC37",
"SIC38",
"SIC39")

#********************************************#

#for (main in 1:1){
title <- paste(tit,tit2[indx],sep="")

x11()
mtit <- title
ytit <- "Percentage Change - Annual Rate"
plot(year1,dy,type="l",ann=0)
title(main=mtit,ylab=ytit) 

source("pv_sup.R")
source("pv_exp.R")
source("st_chang.R")

# Constant Parameter Model #
yy <- as.matrix(dy[(k+1):n])
year1 <- as.matrix(year1[(k+1):n])
t <- n-k
tr <- seq(1,t,1)
ot <- matrix(1,t,1) 
xar <- x_ar(dy,k)
x <- cbind(xar,ot)
rho <- qr.solve(x,yy)
e0 <- yy-x%*%rho
ee0 <- t(e0)%*%e0
sig <- ee0/(t-k-1)
xx <- solve(t(x)%*%x) 
xe <- x*(e0%*%matrix(1,1,ncol(x)))
v <- xx%*%(t(xe)%*%xe)%*%xx
se <- as.matrix(sqrt(diag(v)))

names <- as.matrix(paste("DY(-",seq(1,k,1),")",sep=""))
names <- rbind(names,"C")
names <- format(names,digits=4)
Trho <- format(rho,digits=4)
Tse <- format(se,digits=4)
cat ("\n")
cat (title,"\n")
cat ("Sample for years ", sample_1, " ", sample_2, "\n")
cat ("Number of AR Lags: ", k,"\n")
cat ("\n")
cat ("Constant Parameter AR Model","\n")
cat ("Variable  ","Estimate  ","St Error  ","\n")
for (j in 1:nrow(rho)) cat (names[j],"  ",Trho[j],"  ",Tse[j],"\n")
cat ("sigma     ",sqrt(sig),"\n")
cat ("\n")
cat ("\n")
cat ("************************************","\n")
cat ("\n")

#***************************************#
# Structural Change Model #

out <- splitest(yy,x,0,trim)
kest <- out$kest
ss <- out$ss
d <- as.matrix((tr<=kest))
y1 <- as.matrix(yy[d]) 
x1 <- as.matrix(x[d%*%matrix(1,1,ncol(x))>0])
x1 <- matrix(x1,nrow(x1)/ncol(x),ncol(x))
t1 <- nrow(y1)
y2 <- as.matrix(yy[(1-d)>0]) 
x2 <- as.matrix(x[(1-d)%*%matrix(1,1,ncol(x))>0])
x2 <- matrix(x2,nrow(x2)/ncol(x),ncol(x))
t2 <- nrow(y2)

# Regime 1 Estimates #
rho1 <- qr.solve(x1,y1)
e1 <- y1-x1%*%rho1
ee1 <- t(e1)%*%e1
sig1 <- ee1/(t1-k-1)
q1 <- t(x1)%*%x1 
xx <- solve(q1) 
q1 <- q1/t1
xe1 <- x1*(e1%*%matrix(1,1,ncol(x1)))
v1 <- xx%*%(t(xe1)%*%xe1)%*%xx
se1 <- as.matrix(sqrt(diag(v1)))
sumrho <- 1-sum(rho1[1:k])
mu1 <- rho1[k+1]/sumrho
h1 <-  rbind((matrix(1,k,1)*mu1),1)/sumrho
vmu1 <- t(h1)%*%v1%*%h1

# Regime 2 Estimates #
rho2 <- qr.solve(x2,y2)
e2 <- y2-x2%*%rho2
ee2 <- t(e2)%*%e2
sig2 <- ee2/(t2-k-1)
q2 <- t(x2)%*%x2 
xx <- solve(q2) 
q2 <- q2/t2
xe2 <- x2*(e2%*%matrix(1,1,ncol(x2)))
v2 <- xx%*%(t(xe2)%*%xe2)%*%xx
se2 <- as.matrix(sqrt(diag(v2)))
mu2 <- rho2[k+1]/sumrho
h2 <-  rbind((matrix(1,k,1)*mu2),1)/sumrho
vmu2 <- t(h2)%*%v2%*%h2

tm <- (mu2-mu1)/sqrt(vmu1+vmu2)
e <- rbind(e1,e2)

# Confidence Interval for Breakdate #
delta <- rho2-rho1
xsi <- (t(delta)%*%q2%*%delta)/(t(delta)%*%q1%*%delta)
phi <- xsi*sig2/sig1
out <- bai_qnt(xsi,phi,1-conf)
c1 <- out$c1
c2 <- out$c2 
ll <- (t(delta)%*%q1%*%delta)/sig1
k1 <- kest-1-ceiling(c2/ll)
k2 <- kest+1-floor(c1/ll)
k1 <- max(rbind(k1,1))
k2 <- min(rbind(k2,t))

# Graph SSE #
ss <- ss/t
ytit <- "Residual Variance"
xtit <- "Breakdate"
xtics <- cbind(1950,2000) 
mtit <- rbind("Least Squares Breakdate Estimation:",
        "Residual Variance as a Function of Breakdate")
x11()
plot(year1,ss,type="l",xlim=xtics,ann=0)
title(main=mtit,ylab=ytit,xlab=xtit) 

Trho1 <- format(rho1,digits=4)
Tse1 <- format(se1,digits=4)
Trho2 <- format(rho2,digits=4)
Tse2 <- format(se2,digits=4)

cat ("\n")
cat ("Complete Structural Change AR Model", "\n")
cat ("\n")
cat ("Estimated Breakdate ", year1[kest], "\n")
cat ("Confidence Interval ", year1[k1]," ",year1[k2],"\n")
cat ("Confidence % ", conf,"\n")
cat ("Trimming Percentage  ", trim, "\n")
cat ("t test for equality of means ", tm,"\n")
cat ("\n")
cat ("Regime 1", "\n")
cat ("Mean     ", mu1," ",sqrt(vmu1),"\n")
cat ("Variable  ","Estimate  ","St Error  ","\n")
for (j in 1:nrow(rho1)) cat (names[j],"  ",Trho1[j],"  ",Tse1[j],"\n")
cat ("sigma     ",sqrt(sig1),"\n")
cat ("\n")
cat ("Regime 2", "\n")
cat ("Mean     ", mu2," ",sqrt(vmu2),"\n")
cat ("Variable  ","Estimate  ","St Error  ","\n")
for (j in 1:nrow(rho2)) cat (names[j],"  ",Trho2[j],"  ",Tse2[j],"\n")
cat ("sigma     ",sqrt(sig2),"\n")
cat ("\n")
cat ("************************************","\n")
cat ("\n")
cat ("\n")

#***************************************#
# Testing for Structural Change #
# Compute W Statistics #
out <- supw(yy,x,trimsup,trimexp)
f <- out$fsup
supf <- out$supf
expf <- out$expf
kest <- out$kest

# Graph F Statistics #
cv <- cbind((rbind((matrix(0,k+1,1)+qa_crit(1,trimsup)),qa_crit(k+1,trimsup))),
            (rbind((matrix(0,k+1,1)+chi_crit(1)),chi_crit(k+1))))
names2 <- rbind(names,"Joint")
xtit <- "Breakdate"
ytit <- "Chow Test"
for (j in 1:k){
  x11()
  mtit <- paste("Testing for Structural Change",title," ",names2[j],sep="")
  z1 <- matrix(1,t,1)*cv[j,1]
  z2 <- matrix(1,t,1)*cv[j,2]
  ytics <- range(rbind(as.matrix(f[,j]),z1,z2)) 
  plot(year1,f[,j],type="l",ylim=ytics,lty=1,col=1,ann=0)
  lines(year1,z1,lty=2,col=2)
  lines(year1,z2,lty=3,col=3)
  title(main=mtit,ylab=ytit,xlab=xtit)  
  legend("topright",c("Chow Test Sequence","Andrews Critical Value",
                      "Critical Value"),lty=c(1,2,3),col=c(1,2,3)) 
}
x11()
mtit <- paste("Testing for Structural Change",title," ","Intercept",sep="")
z1 <- matrix(1,t,1)*cv[k+1,1]
z2 <- matrix(1,t,1)*cv[k+1,2]
ytics <- range(rbind(as.matrix(f[,k+1]),z1,z2)) 
plot(year1,f[,k+1],type="l",ylim=ytics,lty=1,col=1,ann=0)
lines(year1,z1,lty=2,col=2)
lines(year1,z2,lty=3,col=3)
title(main=mtit,ylab=ytit,xlab=xtit)  
legend("topright",c("Chow Test Sequence","Andrews Critical Value",
                    "Critical Value"),lty=c(1,2,3),col=c(1,2,3)) 
x11()
mtit <- rbind("Testing for Structural Change of Unknown Timing:",
              "Chow Test Sequence as a Function of Breakdate")
ytics <- cbind(0,24)
xtics <- cbind(1950,2000)
z1 <- matrix(1,t,1)*cv[k+2,1]
z2 <- matrix(1,t,1)*cv[k+2,2]
plot(year1,f[,k+2],type="l",ylim=ytics,xlim=xtics,lty=1,col=1,ann=0)
lines(year1,z1,lty=2,col=2)
lines(year1,z2,lty=3,col=3)
title(main=mtit,ylab=ytit,xlab=xtit)  
legend("topright",c("Chow Test Sequence","Andrews Critical Value",
                    "Critical Value"),lty=c(1,2,3),col=c(1,2,3)) 

# Andrews' P-Values #
pvs1 <- matrix(0,k+2,1) 
pve1 <- matrix(0,k+2,1)
for (j in 1:(k+1)){ 
  pvs1[j,] <- pv_sup(supf[j],1,trimsup)
  pve1[j,] <- pv_exp(expf[j],1,trimexp)
}
pvs1[k+2,] <- pv_sup(supf[k+2],k+1,trimsup)
pve1[k+2,] <- pv_exp(expf[k+2],k+1,trimexp)

cat ("\n")
cat ("Test Breakdates","\n")
cat ("\n")
cat ("Variable","BreakEst","\n")
cat ("\n")
for (j in 1:nrow(names2)) cat (names2[j],"  ",year1[kest[j]],"\n")
cat ("\n")
cat ("\n")
cat ("Trimming Percentage for SupW ", trimsup,"\n")
cat ("Trimming Percentage for ExpW ", trimexp,"\n")
cat ("\n")
cat ("\n")

statfs <- cbind(supf,pvs1)
statfe <- cbind(expf,pve1)
vnames <- cbind("Variable","SupW","AsyP")
vnamee <- cbind("Variable","ExpW","AsyP")

if (bootw >0){
  cat ("Number of Fixed Regressor Bootstrap Replications ", bootw,"\n")
  # Fixed Regressor Bootstrap #
  boots <- matrix(0,bootw,k+2)
  boote <- matrix(0,bootw,k+2) 
  for (b in 1:bootw){
    yb <- e*as.matrix(rnorm(t))
    out <- supw(yb,x,trimsup,trimexp)
    boots[b,] <- t(out$supf)
    boote[b,] <- t(out$expf)
  }
  pvs2 <- colMeans(boots>(matrix(1,bootw,1)%*%t(supf)))
  pve2 <- colMeans(boote>(matrix(1,bootw,1)%*%t(expf)))
  statfs <- cbind(statfs,pvs2)
  statfe <- cbind(statfe,pve2)
  vnames <- cbind(vnames,"FixRegP")
  vnamee <- cbind(vnamee,"FixRegP")
}

# Print Test Results #
Tvnames <- format(vnames,digits=4)
Tvnamee <- format(vnamee,digits=4)
Tstatfs <- cbind(format(names2,digits=4),format(statfs,digits=4))
Tstatfe <- cbind(format(names2,digits=4),format(statfe,digits=4))
cat ("SupW Tests on Regression Parameters","\n") 
cat (Tvnames,"\n")
for (j in 1:nrow(statfs)) cat(Tstatfs[j,],"\n")
cat ("\n")
cat ("ExpW Tests on Regression Parameters","\n")
cat (Tvnamee,"\n")
for (j in 1:nrow(statfe)) cat(Tstatfe[j,],"\n")


#***********************************#
# Variance Equation  - LS Residuals #
#***********************************#

e2 <- e0*e0
out <- splitest(e2,ot,0,trim)
kest <- out$kest 
ss <- out$ss
d <- as.matrix((tr<=kest))
xx <- cbind(d,(1-d))
sigs <- qr.solve(xx,e2)
u <- e2-xx%*%sigs
u1 <- as.matrix(u[1:kest])
u2 <- as.matrix(u[kest+1:t])
t1 <- kest 
t2 <- t-kest
sig1 <- (t(u1)%*%u1)/(t1-k-1)
sig2 <- (t(u2)%*%u2)/(t2-k-1)
out <- bai_qnt(1,1,1-conf)
c1 <- out$c1
c2 <- out$c2
ll <- sig1/((sigs[1]-sigs[2])^2)
k1 <- kest-1-ceiling(c2*ll)
k2 <- kest+1-floor(c1*ll)
k1 <- max(rbind(k1,1))
k2 <- min(rbind(k2,t))

cat ("\n")
cat ("********************************************","\n")
cat ("\n")
cat ("\n")
cat ("Break in Error Variance ", "\n")
cat ("Estimated Breakdate ", year1[kest],"\n")
cat ("Confidence Interval ", year1[k1]," ",year1[k2],"\n")
cat ("Confidence % ", conf, "\n")
cat ("Regime 1 SD ", sqrt(sigs[1]),"\n")
cat ("Regime 2 SD ", sqrt(sigs[2]),"\n")
cat ("\n")
cat ("\n")
cat ("********************************************","\n")
cat ("\n")
cat ("\n")

# Graph SSE #
ss <- ss/t
ytit <- "Residual Variance"
xtit <- "Breakdate"
mtit <- rbind(title,"Estimation of Breakdate in Variance")
x11()
plot(year1,ss,type="l",ann=0)
title(main=mtit,ylab=ytit,xlab=xtit)   

# Test Constancy of Variance #
out <- supw(e2,ot,trimsup,trimexp)
f <- out$fsup
supf <- out$supf
expf <- out$expf
kest <- out$kest

supf <- supf[1]
expf <- expf[1]
kest <- kest[1]
tests <- rbind(supf,expf)
pvvar <- rbind(pv_sup(supf,1,trimsup),pv_exp(expf,1,trimexp))
cvvar <- qa_crit(1,trimsup)

cat ("\n")
cat ("Tests for Constancy of Error Variance ","\n")
cat ("\n")
cat ("\n")
cat ("     ","Stat   ","AsyP","\n")
cat ("SupF ",tests[1],pvvar[1],"\n")
cat ("ExpF ",tests[2],pvvar[2],"\n")
cat ("Test Breakdate ", year1[kest],"\n")
cat ("\n")
cat ("\n")
cat ("********************************************","\n")
cat ("\n")
cat ("\n")

# Graph F Statistics #
x11()
mtit <- rbind(title,"Wald Sequence for Variance")
z0 <- matrix(1,t,1)*cvvar
ytics <- range(rbind(as.matrix(f[,1]),z0)) 
plot(year1,f[,1],type="l",lty=1,col=1,ylim=ytics,ann=0)
lines(year1,z0,lty=2,col=2)
title(main=mtit)  
legend("topright",c("Wald Sequence","Asymptotic Critical Value"),lty=c(1,2,3),col=c(1,2,3)) 

#**********************************************#
}



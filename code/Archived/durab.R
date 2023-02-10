rm(list=ls())
library(strucchange)
# US labor productivity description
# US labor productivity in the manufacturing/durables sector.


#'durab' is a multivariate monthly time series from 1947(3) to 2001(4) with variables
#y growth rate of the Industrial Production Index to average weekly labor hours in the manufacturing/durables sector,

#lag 1 of the series 'y',
#The data set is available from Bruce Hansen's homepage <URL:
#http://www.ssc.wisc.edu/~bhansen/>. For more information see Hansen (2001).


#     Hansen B. (2001), The New Econometrics of Structural Change:
#     Dating Breaks in U.S. Labor Productivity, _Journal of Economic
#     Perspectives_, *15*, 117-128.

#     Zeileis A., Leisch F., Kleiber C., Hornik K. (2002), Monitoring
#     Structural Change in Dynamic Econometric Models, Report 64, SFB
#     "Adaptive Information Systems and Modelling in Economics and
#     Management Science", Vienna University of Economics, <URL:
#     http://www.wu-wien.ac.at/am/reports.htm#78>.
data(durab)
     ## use AR(1) model as in Hansen (2001) and Zeileis et al. (2002)
durab.model <- y ~ lag

## historical tests
## OLS-based CUSUM process
ols <- efp(durab.model, data = durab, type = "OLS-CUSUM")
plot(ols)
## F statistics
fs <- Fstats(durab.model, data = durab, from = 0.1)
plot(fs)
## F statistics based on heteroskadisticy-consistent covariance matrix
#fsHC <- Fstats(durab.model, data = durab, from = 0.1, cov.type = "HC")
#plot(fsHC)

## monitoring
Durab <- window(durab, start=1964, end = c(1979, 12))
ols.efp <- efp(durab.model, type = "OLS-CUSUM", data = Durab)
newborder <- function(k) 1.5778*k/192
ols.mefp <- mefp(ols.efp, period=2)
ols.mefp2 <- mefp(ols.efp, border=newborder)
Durab <- window(durab, start=1964)
ols.mon <- monitor(ols.mefp)
ols.mon2 <- monitor(ols.mefp2)
plot(ols.mon)
lines(boundary(ols.mon2), col = 2)

v <- ols.mefp2$efpprocess
plot(v)
     ## dating
     bp <- breakpoints(durab.model, data = durab)
     summary(bp)
     plot(summary(bp))

     plot(ols)
     lines(breakpoints(bp, breaks = 1), col = 3)
     lines(breakpoints(bp, breaks = 2), col = 4)
     plot(fs)
### Name: durab
### Title: US Labor Productivity
### Aliases: durab
### Keywords: datasets
require(zoo)
### ** Examples
load(file='/mnt/MyDoc/Dropbox/Research/MntStrBrk/data/SP2001.rda')
write.zoo(SP2001, '/mnt/MyDoc/Dropbox/Research/MntStrBrk/data/SP2001.csv', quote = FALSE, sep = ",")
require(strucchange)
data(durab)
durab = as.data.frame(durab)
write.csv(durab,file='durab.csv')
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
fsHC <- Fstats(durab.model, data = durab, from = 0.1, cov.type = "HC")
plot(fsHC)

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

## dating
bp <- breakpoints(durab.model, data = durab)
summary(bp)
plot(summary(bp))

plot(ols)
lines(breakpoints(bp, breaks = 1), col = 3)
lines(breakpoints(bp, breaks = 2), col = 4)
plot(fs)
lines(breakpoints(bp, breaks = 1), col = 3)
lines(breakpoints(bp, breaks = 2), col = 4)




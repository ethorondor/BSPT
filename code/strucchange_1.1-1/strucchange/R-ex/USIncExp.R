### Name: USIncExp
### Title: Income and Expenditures in the US
### Aliases: USIncExp
### Keywords: datasets

### ** Examples

## These example are presented in the vignette distributed with this
## package, the code was generated by Stangle("strucchange-intro.Rnw")

###################################################
### chunk number 1: data
###################################################
library(strucchange)
data(USIncExp)
plot(USIncExp, plot.type = "single", col = 1:2, ylab = "billion US$")
legend(1960, max(USIncExp), c("income", "expenditures"),
       lty = c(1,1), col = 1:2, bty = "n")

###################################################
### chunk number 2: subset
###################################################
library(strucchange)
data(USIncExp)
library(ts)
USIncExp2 <- window(USIncExp, start = c(1985,12))

###################################################
### chunk number 3: ecm-setup
###################################################
coint.res <- residuals(lm(expenditure ~ income, data = USIncExp2))
coint.res <- lag(ts(coint.res, start = c(1985,12), freq = 12), k = -1)
USIncExp2 <- cbind(USIncExp2, diff(USIncExp2), coint.res)
USIncExp2 <- window(USIncExp2, start = c(1986,1), end = c(2001,2))
colnames(USIncExp2) <- c("income", "expenditure", "diff.income",
                         "diff.expenditure", "coint.res")
ecm.model <- diff.expenditure ~ coint.res + diff.income

###################################################
### chunk number 4: ts-used
###################################################
plot(USIncExp2[,3:5], main = "")

###################################################
### chunk number 5: efp
###################################################
ocus <- efp(ecm.model, type="OLS-CUSUM", data=USIncExp2)
me <- efp(ecm.model, type="ME", data=USIncExp2, h=0.2)

###################################################
### chunk number 6: efp-boundary
###################################################
bound.ocus <- boundary(ocus, alpha=0.05)

###################################################
### chunk number 7: OLS-CUSUM
###################################################
plot(ocus)

###################################################
### chunk number 8: efp-boundary2
###################################################
plot(ocus, boundary = FALSE)
lines(bound.ocus, col = 4)
lines(-bound.ocus, col = 4)

###################################################
### chunk number 9: ME-null
###################################################
plot(me, functional = NULL)

###################################################
### chunk number 10: efp-sctest
###################################################
sctest(ocus)

###################################################
### chunk number 11: efp-sctest2
###################################################
sctest(ecm.model, type="OLS-CUSUM", data=USIncExp2)

###################################################
### chunk number 12: Fstats
###################################################
fs <- Fstats(ecm.model, from = c(1990, 1), to = c(1999,6), data = USIncExp2)

###################################################
### chunk number 13: Fstats-plot
###################################################
plot(fs)

###################################################
### chunk number 14: pval-plot
###################################################
plot(fs, pval=TRUE)

###################################################
### chunk number 15: aveF-plot
###################################################
plot(fs, aveF=TRUE)

###################################################
### chunk number 16: Fstats-sctest
###################################################
sctest(fs, type="expF")

###################################################
### chunk number 17: Fstats-sctest2
###################################################
sctest(ecm.model, type = "expF", from = 49, to = 162, data = USIncExp2)

###################################################
### chunk number 18: mefp
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1), end = c(1989,12))
me.mefp <- mefp(ecm.model, type = "ME", data = USIncExp3, alpha = 0.05)

###################################################
### chunk number 19: monitor1
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1), end = c(1990,12))
me.mefp <- monitor(me.mefp)

###################################################
### chunk number 20: monitor2
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1))
me.mefp <- monitor(me.mefp)
me.mefp

###################################################
### chunk number 21: monitor-plot
###################################################
plot(me.mefp)

###################################################
### chunk number 22: mefp2
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1), end = c(1989,12))
me.efp <- efp(ecm.model, type = "ME", data = USIncExp3, h = 0.5)
me.mefp <- mefp(me.efp, alpha=0.05)

###################################################
### chunk number 23: monitor3
###################################################
USIncExp3 <- window(USIncExp2, start = c(1986, 1))
me.mefp <- monitor(me.mefp)

###################################################
### chunk number 24: monitor-plot2
###################################################
plot(me.mefp)





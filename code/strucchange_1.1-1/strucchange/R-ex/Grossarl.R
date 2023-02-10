### Name: Grossarl
### Title: Marriages and Births in Grossarl
### Aliases: Grossarl
### Keywords: datasets

### ** Examples

data(Grossarl)

## illegitimate births
######################
## lm + MOSUM
plot(Grossarl$fraction)
fm.min <- lm(fraction ~ politics, data = Grossarl)
fm.max <- lm(fraction ~ politics + morals + nuptiality + lag.marriages,
             data = Grossarl)
fm.final <- step(fm.max)
lines(ts(fitted(fm.min), start = 1700), col = 3)
lines(ts(fitted(fm.final), start = 1700), col = 4)
mos.min <- efp(fraction ~ politics, data = Grossarl, type = "OLS-MOSUM")
mos.final <- efp(fraction ~ politics + morals + nuptiality, data = Grossarl,
                 type = "OLS-MOSUM")
plot(mos.min)
lines(mos.final, lty = 2)

## dating
bp <- breakpoints(fraction ~ 1, data = Grossarl, h = 0.1)
summary(bp)
## RSS, BIC, AIC
plot(bp)
plot(0:8, AIC(bp), type = "b")

## probably use 5 (or maybe 6) breakpoints and compare with
## coding of the factors as used by us
##
## politics                   1803      1816 1850
## morals      1736 1753 1771 1803
## nuptiality                 1803 1810 1816      1883
##
## m = 5            1753 1785           1821 1856 1878
## m = 6       1734 1754 1785           1821 1856 1878
##              6    2    5              1    4    3

fm.bp <- lm(fraction ~ breakfactor(breakpoints(bp, breaks = 6)),
            data = Grossarl)

plot(Grossarl$fraction)
lines(fitted(fm.final), col = 3)
lines(fitted(fm.bp), col = 4)


## marriages
############
## lm + MOSUM
plot(Grossarl$marriages)
fm.min <- lm(marriages ~ politics, data = Grossarl)
fm.final <- lm(marriages ~ politics + morals + nuptiality, data = Grossarl)
lines(ts(fitted(fm.min), start = 1700), col = 3)
lines(ts(fitted(fm.final), start = 1700), col = 4)
mos.min <- efp(marriages ~ politics, data = Grossarl, type = "OLS-MOSUM")
mos.final <- efp(marriages ~ politics + morals + nuptiality, data = Grossarl,
                 type = "OLS-MOSUM")
plot(mos.min)
lines(mos.final, lty = 2)

## dating
bp <- breakpoints(marriages ~ 1, data = Grossarl, h = 0.1)
summary(bp)
## RSS, BIC, AIC
plot(bp)
plot(0:8, AIC(bp), type = "b")

## probably use 3 (or maybe 4) breakpoints and compare with
## coding of the factors as used by us
##
## politics                   1803      1816 1850
## morals      1736 1753 1771 1803
## nuptiality                 1803 1810 1816      1883
##
## m = 3       1738                     1813      1875
## m = 4       1738      1794           1814      1875
##              2         4              1         3
fm.bp <- lm(marriages ~ breakfactor(breakpoints(bp, breaks = 4)),
            data = Grossarl)

plot(Grossarl$marriages)
lines(fitted(fm.final), col = 3)
lines(fitted(fm.bp), col = 4)




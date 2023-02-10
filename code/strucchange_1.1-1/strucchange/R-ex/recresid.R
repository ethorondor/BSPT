### Name: recresid
### Title: Recursive Residuals
### Aliases: recresid recresid.default recresid.formula recresid.lm
### Keywords: regression

### ** Examples

x <- rnorm(100)
x[51:100] <- x[51:100] + 2
rr <- recresid(x ~ 1)
plot(cumsum(rr), type = "l")

plot(efp(x ~ 1, type = "Rec-CUSUM"))




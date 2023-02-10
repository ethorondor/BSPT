### Name: covHC
### Title: Heteroskedasticity-Consistent Covariance Matrix Estimation
### Aliases: covHC
### Keywords: htest

### ** Examples

## generate linear regression relationship
## with homoskedastic variances
x <- sin(1:100)
y <- 1 + x + rnorm(100)
## compute usual covariance matrix of coefficient estimates
covHC(y~x, type="const")

sigma2 <- sum(residuals(lm(y~x))^2)/98
sigma2 * solve.crossprod(cbind(1,x))




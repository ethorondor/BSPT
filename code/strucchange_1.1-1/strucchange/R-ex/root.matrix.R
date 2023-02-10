### Name: root.matrix
### Title: Root of a Matrix
### Aliases: root.matrix solve.crossprod solve.pd
### Keywords: algebra

### ** Examples

X <- matrix(c(1,2,2,8), ncol=2)
test <- root.matrix(X)
## control results
X
test %*% test




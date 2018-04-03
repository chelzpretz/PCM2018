## Lab 3: Generalized least squares

library(phytools)
library(geiger)

# make phylogeny 

a=tr <- pbtree(n = 100)
dat <- sim.char(tr, matrix(c(1,0,0,1), 2, 2))[,,1]
colnames(dat) <- c('y', 'x')
#layout(matrix(c(1:2, 1, 2)))
par(mfrow = c(1,2))
plot(tr)
plot(dat)

# making a variance-covariance matrix of the tree

tr.v <- vcv.phylo(tr)
round(tr.v[1:10, 1:10], 2) # shows just the upper lefthand corner of your matrix

# Invert the martix it: the inverted matrix is that matrix that multiplied by the original (uninverted) matrix yields a matrix with diagonals of 1 and off-diagonals of 0

round(solve(tr.v), 2)[1:10, 1:10]
round(solve(tr.v) %*% tr.v, 2)[1:10, 1:10]


# Using the inverted matrix, data can be weighted by how phylogenetically distinct the corresponding tips are. 
# To visualize these weights, multiply the inverted matrix by a column-matrix of ones (for example, matrix(1, 10) is a one-column matrix of 10 1s).

one <- matrix(1, length(tr$tip.label), 1) # a matrix of 1s, 1 column
tr.w <- t(one) %*% solve(tr.v)
a=plot(tr, show.tip.label = F, x.lim = nodeheight(tr,1) * 3) # sets the plot window width to be thrice the tree depth
phydataplot(tr.w[1,] / max(tr.w[1,]) * nodeheight(tr, 1), tr)


# GLS estimate of regression coefficients with the OLS (ordinary least squares) estimates, using equations from Symonds and Blomberg Appendix A

b.gls <- function(tr, dat, y = 'y', x = 'x') {
  #print(dat)
  v.solved <- solve(vcv.phylo(tr))
  one <- matrix(1, length(tr$tip.label))
  X <- cbind(origin = one, x = dat[, x])
  out <- solve(t(X) %*% v.solved %*% X) %*% (t(X) %*% v.solved %*% dat[, y])
  out
}

b.ols <- function(dat, y = 'y', x = 'x') {
  one <- matrix(1, dim(dat)[1])
  X <- cbind(origin = one, x = dat[, x])
  out <- solve(t(X) %*% X) %*% (t(X) %*% dat[, y])
  out
}
message('GLS coefficients (intercept and slope)')
b.gls(tr, dat)
message('OLS coefficients (intercept and slope)')
b.ols(dat)

# Revell and Symonds and Blomberg both discuss estimating phylogenetic signal simultaneously with the regression coefficients. This is really simple to do in R. Here is an example:

library(nlme)
out <- gls(y ~ x, as.data.frame(dat), corPagel(1, tr))
summary(out)

gls.r.squared <- function(x) {
  # based on Judge et al. 1985, eq. 2.3.16
  e = x$resid
  V <- corMatrix(x$modelStruct$corStruct)
  Y = x$resid + x$fitted
  one <- matrix(1, length(Y), 1)
  a <- as.numeric(solve(t(one) %*% solve(V) %*% one) %*% (t(one) %*% solve(V) %*% Y))
  r.squared= 1-(t(e) %*% solve(V) %*% e) / (t(Y-a) %*% solve(V) %*% (Y-a))
  return(r.squared[1,1])
}

## example:
gls.r.squared(out) # probably will be pretty low, as we simulated data with no correlation


# Practice questions
# 1. Play around with the inverse of the covariance matrix (tr.w <- t(one) %*% solve(tr.v)) on different trees and see what taxa are upweighted more. What species contribute the most to the regression? Which weightings make sense, and which are surprising? You can use the code above and just run it for a handful of trees. Try your own data, if you have some.

# 2. Wrestle your own data into R or use the maple dataset we started working with on week 1, and use gls and lm to calculate regression coefficients. Try this in the traditional GLS / PIC manner using corBrownian for the correlation structure, and as Revell suggests, using corPagel and simulating lambda. What estimate of lambda do you get? how, if at all, does this affect your inferences?


# 3. Use formula #4 in Rohlf 2006 to calculate the phylogenetic mean (answer at the end of this tutorial). How is it different from the formula for the regression coefficients? Then, use that formula to estimate the mean for a randomly simulated trait on each of 500 randomly generated trees, each of 50 tips, with a root state of 0. Compare with the non-phylogenetic means for these. How do they compare? which gives the smaller variance? the lower bias?

# 4. Use the functions lm and gls OR the functions we generated for regression coefficients to compare the OLS and GLS slope on a Brownian-motion generated pair of traits on each of 500 randomly generated trees of 50 tips each, with a slope and root state of 0. How do they compare? which gives the smaller variance? the lower bias?




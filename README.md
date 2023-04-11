
# sptwed: An R-package to perform inference for spatial Tweedie Compound Poisson-gamma Double generalized linear models

<!-- badges: start -->
![Maintainer](https://img.shields.io/badge/maintainer-arh926-blue)
<!-- badges: end -->

The goal of `sptwed` is to carry out statistical inference for spatial Tweedie Compound Poisson-gamma Double generalized linear models. It leverages a co-ordinate descent algorithm for estimating the coefficients. It contains the following functions

Function | Description
:--------|:-----------
`crossvalPll_sptw.R` | K-fold cross-validation (main callable function)
`pathMM_sptw.R` | Warm-start (supporting function)
`spatial_tweedie.R` | Co-ordinate descent (supporting function)

## Installation

You can install the development version of sptwed from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("arh926/sptwed")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
set.seed(2022)
require(tweedie)

# Generate Data
N = 1e4
L = 1e2

coords = matrix(runif(2*L), nc=2)
par(mfcol=c(1,1))
# plot(coords)
sigma2.true = 1.5
phis.true = 3
Delta = as.matrix(dist(coords))
Sigma = sigma2.true*exp(-phis.true*Delta)
w.true = MASS::mvrnorm(1, rep(0,L), Sigma)

if(N > L) index = sample(1:L, N, replace = T) else if(N == L) index = sample(1:L, N, replace = F)
# Design matrices
x = z = cbind(1, rnorm(N), rnorm(N), rnorm(N))
p = ncol(x)
q = ncol(z)

# bockwise spatial effect
# sp_eff = apply(coords,1, function(x){
#   if(x[2]<0.25) return(-3)
#   if(x[2]>0.25 & x[2]<0.5) return(-1)
#   if(x[2]>0.5 & x[2]<0.75) return(1)
#   else return(3)
# })

theta = rnorm(N, -0.16, 0.02)
mu_sim = 4/theta^2 #* exp(w.true[index])
phi_sim = runif(N,10,15) # change this for increased zeros
beta.true = solve(crossprod(x,x)) %*% crossprod(x,log(mu_sim))
gamma.true = solve(crossprod(z,z)) %*% crossprod(z,log(phi_sim))
mu_sim.sp = 4/theta^2 * exp(w.true[index])

# # covariates
# beta0 = 1
# beta1 = 1.5
# beta2 = 1.1
# beta3 = 1.4
# beta.true = c(beta0, beta1,beta2,beta3)
# mu_sim.sp = exp(x%*%beta.true + w.true[index])
# gamma0 = 1
# gamma1 = 0.5
# gamma2 = 0.1
# gamma3 = 1.1
# gamma.true = c(gamma0, gamma1,gamma2, gamma3)
# phi_sim = exp(z%*%gamma.true)

xi.true = 1.5

y_sim = rtweedie(N, xi = xi.true, mu = mu_sim.sp, phi = phi_sim)
sum(y_sim == 0)/N # proportion of zeros
par(mfcol=c(1,1)); hist(y_sim) # histogram
y.mean_sp = aggregate(y_sim, list(index), sum)[,2]
par(mfrow=c(2,2))
hist(log(y.mean_sp), probability = T, ylim=c(0,0.5),
     main = "",
     xlab = "Log of Spatially aggregated response",
     col="lightblue")
lines(density(log(y.mean_sp)))
hist(w.true, probability = T, ylim=c(0,0.5),
     main = "",
     xlab = "Spatial Effect",
     col="lightblue")
lines(density(w.true))
boxplot(y_sim~round(w.true[index],3), ylab = "Response", xlab = "Spatial Effect")
grid()
plot(w.true,
     y.mean_sp,
     xlab = "Spatial Effect",
     ylab = "Spatially Aggregated Response")
lines(lowess(y.mean_sp~w.true), col="red")
grid()


# spatial plot for w and log(y+1)
mat <- matrix(c(1,2,3,4), nr=1,nc=4, byrow=T)
layout(mat,
       widths = rep(c(3,1.5),2),
       heights = rep(c(3,3),2))
sp_plot(data_frame = cbind(coords,w.true), points.plot = T, contour.plot = T, legend = T)
sp_plot(data_frame = cbind(coords,log(y.mean_sp+1)), points.plot = T, contour.plot = T, legend = T)

cor(w.true, log(y.mean_sp+1))

adjM = apply(Delta, 1, function(s){
  s[s < 0.15] = 1
  s[s > 0.15 & s != 1] = 0
  s
})
diag(adjM) = 0
par(mfcol=c(1,1))
sp_plot(data_frame = cbind(coords,w.true), points.plot = T, contour.plot = T, legend = F)
for(i in 1:L){
  id = which(adjM[i,] == 1)
  for(j in 1:length(id)){
    lines(rbind(coords[i,], coords[id[j],]),
          col="darkgreen",
          lwd = 1.5)
  }
}
degM = diag(as.vector(rowSums(adjM)))

beta.init = rep(0, p)
gamma.init = rep(0,q)
alpha.init = rep(0,nrow(adjM))

p.tw = 1.2
tol = 1e-6
miter = 1e4
l1_seq <- exp(seq(-5,5,length.out = 10))
l2_seq <- exp(seq(-5,5,length.out = 10))
lapMat <- degM - adjM

full_id <- fold_split(K=3,index = index)
fold_1 <- as.numeric(unlist(lapply(full_id, function(x) x[[1]])))
fold_2 <- as.numeric(unlist(lapply(full_id, function(x) x[[2]])))
fold_3 <- as.numeric(unlist(lapply(full_id, function(x) x[[3]])))
id.list <- list(fold1=fold_1, fold2=fold_2, fold3=fold_3)
names(alpha.init) = 1:L
cvM <- crossvalPll_sptw(K=3,
                        y=y_sim,
                        X=x,
                        Z=z,
                        index=index,
                        beta.init = beta.init,
                        gamma.init = gamma.init,
                        alpha.init = alpha.init,
                        id.list=id.list,
                        l1_seq=l1_seq,
                        l2_seq=l2_seq,
                        lapMat=lapMat,
                        miter=miter,
                        tol=tol,
                        p=p.tw,
                        verbose=T)
devM <- (cvM[[1]]$dev+cvM[[2]]$dev+cvM[[3]]$dev)
pM <- (cvM[[1]]$eff.p+cvM[[2]]$eff.p+cvM[[3]]$eff.p)/3
arr.min <- which(devM==min(devM),arr.ind = T)

fit_sptw <- spatial_tweedie(y = y_sim,
                            X = x,
                            Z = z,
                            index = index,
                            index.y.0 = y_sim==0,
                            beta.init = beta.init,
                            gamma.init = gamma.init,
                            alpha.init = alpha.init,
                            pen_mat = l1_seq[arr.min[1]]*diag(nrow(lapMat))+l2_seq[arr.min[2]]*lapMat, #
                            p = pM[arr.min[1],arr.min[2]],
                            tol = tol,
                            miter = miter,
                            inf=T,
                            p.update = F)
beta.est = fit_sptw$optim_pars$beta
gamma.est = fit_sptw$optim_pars$gamma
w.est = fit_sptw$optim_pars$alpha
# sp_plot(data_frame = cbind(coords,w.est), points.plot = T, contour.plot = T, legend = F)

```


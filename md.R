######################### codes for the md datas######################
##################### run it on clusters#################################
##loading required packages and functions
library(gss)
load("quant.RData")
source("ssa.R")
source("esps-high.R")
source("ssanovas.R")

################ASP
nob = dim(dat)[1]
sam.size <- ceiling(50*(nob^{1/4}))
ptm <- proc.time()
esps.re <- esps(y~.+x1:x17+x3:x31+x6:x23+x10:x23+x11:x22+x18:x33+x22:x33,
                dat, sam.size, r=3)
ttt1 <- proc.time()-ptm
r <- 3
p <- esps.re$p
pro.gcv.lamb <- esps.re$lambda*((nob/sam.size)^(-r/(p*r+1)))
lambda.g <- log10(nob*pro.gcv.lamb)
theta <- esps.re$theta
nbasis = round(nob^{2/9})
ptm <- proc.time()
prop.fit <- ssa(y~.+x1:x17+x3:x31+x6:x23+x10:x23+x11:x22+x18:x33+x22:x33,
     data=dat, lambda = lambda.g, theta=theta, nbasis=nbasis, seed=1, alpha=1.0)
proc.time() - ptm
fitted.dt <- fitted(prop.fit)
sqrt(mean((dat$y-fitted.dt)^2))

#################SKIP
nob <- dim(dat)[1]
nbasis <- round(nob^{2/9})
ptm <- proc.time()
fit.skip <- ssanova(y~.+x1:x17+x3:x31+x6:x23+x10:x23+x11:x22+x18:x33+x22:x33,
   data=dat, nbasis=nbasis, seed=1, alpha=1.0, skip.iter=TRUE)
proc.time() - ptm
fitted.skip <- fitted(fit.skip)
sqrt(mean((dat$y-fitted.skip)^2))

#########################five fold cross validation##############
#########ASP
set.seed(1)
inx = sample(dim(dat)[1], dim(dat)[1])
nob = dim(dat)[1]
ys = c(0, 178647, 357294, 535941, 714588, 893238)
mse.cv = c()
for(cv in 1:5){
  dt.inx=inx[(ys[cv]+1):ys[cv+1]]
  cv.dt = dat[setdiff(1:nob, dt.inx), ]
  cv.n = dim(cv.dt)[1]
  pre.dt = dat[dt.inx, -1]
  sam.size <- ceiling(50*((cv.n)^{1/4}))
  esps.re <- esps(y~.+x1:x17+x3:x31+x6:x23+x10:x23+x11:x22+x18:x33+x22:x33,
    cv.dt, sam.size, r=3)
  r = 3
  p = esps.re$p
  cat(p, "\n")
  pro.gcv.lamb = esps.re$lambda*((cv.n/sam.size)^(-r/(p*r+1)))
  lambda.g = log10(cv.n*pro.gcv.lamb)
  theta = esps.re$theta
  nbasis = round(cv.n^{2/9})
  prop.fit = ssa(y~.+x1:x17+x3:x31+x6:x23+x10:x23+x11:x22+x18:x33+x22:x33,
    data=cv.dt, lambda = lambda.g, theta=theta, nbasis=nbasis, seed=1, alpha=1.0)
  pre.cv = predict(prop.fit, newdata=pre.dt)
  mse.cv[cv] = mean((pre.cv - dat[dt.inx, 1])^2)
}
sqrt(mse.cv)

#########SKIP
set.seed(1)
inx <- sample(dim(dat)[1], dim(dat)[1])
nob <- dim(dat)[1]
ys = c(0, 178647, 357294, 535941, 714588, 893238)
mse.cv = c()
for(cv in 1:5){
  dt.inx=inx[(ys[cv]+1):ys[cv+1]]
  cv.dt = dat[setdiff(1:nob, dt.inx), ]
  cv.n = dim(cv.dt)[1]
  pre.dt = dat[dt.inx, -1]
  nbasis = round(cv.n^{2/9})
  gcv.fit = ssanova(y~.+x1:x17+x3:x31+x6:x23+x10:x23+x11:x22+x18:x33+x22:x33,
    data=cv.dt, nbasis=nbasis, seed=1, alpha=1.0, skip.iter=TRUE)
  pre.cv = predict(gcv.fit, newdata=pre.dt)
  mse.cv[cv] = mean((pre.cv - dat[dt.inx, 1])^2)
}
sqrt(mse.cv)

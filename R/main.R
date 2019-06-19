library(gss)
source("ssa.R")
source("curves.R")
source("esps.R")

sample.size <- 20000
set.seed(111)
x1 <- runif(sample.size)
x2 <- runif(sample.size)
x3 <- runif(sample.size)
#univariate function
dat.uni <- uni.func(x1,5)
y.uni <- dat.uni$y
y.unit <- dat.uni$true.y
sam.size <- ceiling(50*(sample.size^{1/4}))
data <- data.frame(y=y.uni,x=x1)
esps.re <- esps(y~x,data,sam.size)
r <- 4
p <- 2
esps.lamb <- esps.re$lambda*(sample.size/sam.size)^(-r/(p*r+1))
lambda.pro <- log10(sample.size*esps.lamb)
theta <- esps.re$theta
uni.fit <- ssa(y~x,data=data,lambda=lambda.pro,theta=theta,seed=111)
mean((fitted(uni.fit)-y.unit)^2)
gss.fit <- ssanova(y~x,data=data,alpha=1.0,seed=111)
mean((fitted(gss.fit)-y.unit)^2)

#multivariate function
data.mult <- multi.func(x1,x2,x3,5)
y.mult <- data.mult$y
y.multt <- data.mult$true.y
sam.size <- ceiling(50*(sample.size^{1/4}))
data <- data.frame(y=y.mult,x1,x2,x3)
esps.re <- esps(y~x1+x2+x3,data,sam.size)
r <- 3
p <- 2
esps.lamb <- esps.re$lambda*(sample.size/sam.size)^(-r/(p*r+1))
lambda.pro <- log10(sample.size*esps.lamb)
theta <- esps.re$theta
mult.fit <- ssa(y~x1+x2+x3,data=data,lambda=lambda.pro,theta=theta,seed=111)
mean((fitted(mult.fit)-y.multt)^2)
gss.fit <- ssanova(y~x1+x2+x3,data=data,alpha=1.0,seed=111)
mean((fitted(gss.fit)-y.multt)^2)


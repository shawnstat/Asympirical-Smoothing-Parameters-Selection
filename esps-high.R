esps = function(form, data, sam.size=NULL, iter.c=10){
##Author: Xiaoxiao Sun
  mf = model.frame(form, data = data)
  resp = model.extract(mf, "response")
  sam.theta <- NULL
  sam.lamb <- NULL
  obs = dim(data)[1]
  for(t in 1:iter.c){
    set.seed(t+10)
    sam.indx <- sample(obs, sam.size)
    sam.dat <- data[sam.indx, ]
    sam.fit <- ssanovas(form, data=sam.dat, seed=t, alpha=1.0)
    sam.theta <- rbind(sam.theta, sam.fit$theta)
    sam.lamb[t] <- sam.fit$nlambda
  }
  me.lamb <- median(sam.lamb)
  me.theta <- apply(sam.theta, 2, median)
  me.lamb.org = 10^me.lamb/sam.size
p.c <- 1
list(lambda=me.lamb.org, theta=me.theta, p=p.c)
}

beta.func <- function(t, p, q){
  gamma(p+q)/(gamma(p)*gamma(q))*t^{p-1}*(1-t)^{q-1}
}

uni.func <- function(x,SNR){
    true.f <- (1/3)*beta.func(x,20,5) + (1/3)*beta.func(x,12,12) + (1/3)*beta.func(x,7,30)
    y <- true.f + rnorm(length(x),0,sd(true.f)/SNR)
    list(y=y,true.y=true.f)
} 


multi.func <- function(x1, x2, x3, SNR){
  n <- length(x1)
  true.f <- 10*sin(pi*x1) + exp(3*x2) + 10^6*x3^11*(1-x3)^6 + 10^4*x3^3*(1-x3)^10
  sd.f <- sd(true.f)
  f <- true.f + rnorm(n,0,sd.f/SNR)
  list(y=f,true.y=true.f)
}


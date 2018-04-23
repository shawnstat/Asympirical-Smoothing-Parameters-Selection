## The ssa function takes the pre-defined smoothing parameter as input.
## Author: Xiaoxiao Sun
## 10/31/2016
##lambda: log10(nlambda), from ssanova function
##theta: log10(theta), from ssanova funciton

ssa <- function(formula,type=NULL,data=list(),lambda,theta,weights,subset,
	                offset,na.action=na.omit,method="v",alpha=1.4,varht=1,
                    id.basis=NULL,nbasis=NULL,seed=NULL)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$type <- mf$method <- mf$varht <- NULL
    mf$alpha <- mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$lambda <- mf$theta <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    wt <- model.weights(mf)
    ## Generate sub-basis
    nobs <- dim(mf)[1]
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=wt)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in ssanova: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Generate terms
    term <- mkterm(mf,type)
    ## Generate s, r, and y
    s <- r <- NULL
    nq <- 0
    for (label in term$labels) {
        if (label=="1") {
            s <- cbind(s,rep(1,len=nobs))
            next
        }
        x <- mf[,term[[label]]$vlist]
        x.basis <- mf[id.basis,term[[label]]$vlist]
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi)
                s <- cbind(s,phi$fun(x,nu=i,env=phi$env))
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),c(nobs,nbasis,nq))
            }
        }
    }
    if (is.null(r))
        stop("gss error in ssanova: use lm for models with only unpenalized terms")
    if (qr(s)$rank<dim(s)[2])
        stop("gss error in ssanova: unpenalized terms are linearly dependent")
    ## Prepare the data
    y <- model.response(mf,"numeric")
    offset <- model.offset(mf)
    if (!is.null(offset)) {
        term$labels <- c(term$labels,"offset")
        term$offset <- list(nphi=0,nrk=0)
        y <- y - offset
    }
    if (!is.null(wt)) wt <- sqrt(wt)
    ## Fit the model
    if (nq==1) {
        r <- r[,,1]
        z <- sreg(s,r,r[id.basis,],y,lambda,theta,wt,method,alpha,varht)
    }
    else z <- mreg(s,r,id.basis,y,lambda,theta,wt,method,alpha,varht)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    desc <- rbind(desc,apply(desc,2,sum))
    rownames(desc) <- c(term$labels,"total")
    colnames(desc) <- c("Unpenalized","Penalized")
    ## Return the results
    obj <- c(list(call=match.call(),mf=mf,terms=term,desc=desc,alpha=alpha,
                  id.basis=id.basis),z)
    class(obj) <- c("ssanova")
    obj
}

## Fit Single Smoothing Parameter (Gaussian) REGression
sreg <- function(s,r,q,y,lambda,theta,wt,method,alpha,varht)
{
    alpha <- abs(alpha)
    ## get dimensions
    nobs <- nrow(r)
    nxi <- ncol(r)
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    nxiz <- nxi
    nn <- nxiz + nnull
    if (!is.null(wt)) {
        y <- wt*y
        s <- wt*s
        r <- wt*r
    }

    ## weighted q using lambda and theta
    q.wk <- 10^(lambda+theta)*q

    z <- .Fortran("reg",
            as.double(cbind(s,10^theta*r)), as.integer(nobs), as.integer(nnull),
            as.double(q.wk), as.integer(nxiz), as.double(y),
            as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
            as.double(alpha), varht=as.double(varht),
            score=double(1), dc=double(nn),
            as.double(.Machine$double.eps),
            chol=double(nn*nn), double(nn),
            jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
            wk=double(3*nobs+nnull), rkv=integer(1), info=integer(1),
            PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
    if (z$info) stop("gss error in ssanova: evaluation of GML score fails")
    assign("fit",z[c(1:5,7)],inherits=TRUE)

    se.q.wk <- 10^theta*q
    se.aux <- seaux(s,10^theta*r,se.q.wk,lambda,fit)
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    c(list(method=method,theta=theta,c=c,d=d,nlambda=lambda),
    	fit[-3],list(se.aux=se.aux))
}

## Fit Multiple Smoothing Parameter (Gaussian) REGression
mreg <- function(s,r,id.basis,y,lambda,theta,wt,method,alpha,varht)
{

    alpha <- abs(alpha)
    ## get dimensions
    nobs <- nrow(r)
    nxi <- ncol(r)
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    nxiz <- nxi
    nn <- nxiz + nnull
    nq <- dim(r)[3]
    ## weighted q using lambda and theta
    r.wk0 <- 0
    for (i in 1:nq) {
        r.wk0 <- r.wk0 + 10^theta[i]*r[,,i]
    }
    qq.wk <- r.wk0[id.basis,]
    q.wk <- 10^lambda*qq.wk

    if (!is.null(wt)) {
        y.wk <- wt*y
        s.wk <- wt*s
        r.wk0 <- wt*r.wk0
    } else {
        y.wk <- y
        s.wk <- s
    }

    z <- .Fortran("reg",
            as.double(cbind(s.wk,r.wk0)), as.integer(nobs), as.integer(nnull),
            as.double(q.wk), as.integer(nxiz), as.double(y.wk),
            as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
            as.double(alpha), varht=as.double(varht),
            score=double(1), dc=double(nn),
            as.double(.Machine$double.eps),
            chol=double(nn*nn), double(nn),
            jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
            wk=double(3*nobs+nnull), rkv=integer(1), info=integer(1),
            PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
    if (z$info) stop("gss error in ssanova: evaluation of GML score fails")
    assign("fit",z[c(1:5,7)],inherits=TRUE)

    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    se.qq.wk <- r.wk[id.basis,]
    se.q.wk <- se.qq.wk

    if (!is.null(wt)) {
        s <- wt*s
        r.wk <- wt*r.wk
    }

    se.aux <- seaux(s,r.wk,se.q.wk,lambda,fit)
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    c(list(method=method,theta=theta[1:nq],c=c,d=d,nlambda=lambda),fit[-3],list(se.aux=se.aux))
}

## Auxiliary Quantities for Standard Error Calculation
seaux <- function(s,r,q,nlambda,fit)
{
    nnull <- dim(s)[2]
    nn <- nnull +  dim(q)[1]
    zzz <- eigen(q,symmetric=TRUE)
    rkq <- min(fit$rkv-nnull,sum(zzz$val/zzz$val[1]>sqrt(.Machine$double.eps)))
    val <- zzz$val[1:rkq]
    vec <- zzz$vec[,1:rkq,drop=FALSE]
    if (nnull) {
        wk1 <- qr(s)
        wk1 <- (qr.qty(wk1,r%*%vec))[-(1:nnull),]
    }
    else wk1 <- r%*%vec
    wk2 <- t(t(wk1)/sqrt(val))
    wk2 <- t(wk2)%*%wk2
    wk2 <- solve(wk2+diag(10^nlambda,dim(wk2)[1]),wk2)
    wk2 <- (wk2+t(wk2))/2
    wk2 <- t(wk2/sqrt(val))/sqrt(val)
    wk2 <- diag(1/val,dim(wk2)[1])-wk2
    z <- .Fortran("regaux",
                  as.double(fit$chol), as.integer(nn),
                  as.integer(fit$jpvt), as.integer(fit$rkv),
                  drcr=as.double(t(cbind(s,r))%*%r%*%vec), as.integer(rkq),
                  sms=double(nnull^2), as.integer(nnull), double(nn*nnull),
                  PACKAGE="gss")[c("drcr","sms")]
    drcr <- matrix(z$drcr,nn,rkq)
    dr <- drcr[1:nnull,,drop=FALSE]
    sms <- 10^nlambda*matrix(z$sms,nnull,nnull)
    wk1 <- matrix(0,nnull+rkq,nnull+rkq)
    wk1[1:nnull,1:nnull] <- sms
    wk1[1:nnull,nnull+(1:rkq)] <- -t(t(dr)/val)
    wk1[nnull+(1:rkq),nnull+(1:rkq)] <- wk2
    z <- chol(wk1,pivot=TRUE)
    wk1 <- z
    rkw <- attr(z,"rank")
    while (wk1[rkw,rkw]<wk1[1,1]*sqrt(.Machine$double.eps)) rkw <- rkw-1
    wk1[row(wk1)>col(wk1)] <- 0
    if (rkw<nnull+rkq)
        wk1[(rkw+1):(nnull+rkq),(rkw+1):(nnull+rkq)] <- diag(0,nnull+rkq-rkw)
    hfac <- wk1
    hfac[,attr(z,"pivot")] <- wk1
    list(vec=vec,hfac=hfac)
}

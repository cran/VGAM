# These functions are
# Copyright (C) 1998-2006 T.W. Yee, University of Auckland. All rights reserved.




VGAM.weights.function = function(w, M, n) {
    ncolw = ncol(as.matrix(w))
    if(ncolw == 1) {
        wz = matrix(w, nrow=n, ncol=M) # w_i * diag(M)
    } else if(ncolw == M) {
        wz = as.matrix(w)
    } else if(ncolw < M && M > 1) {
        stop("ambiguous input for weights")
    } else if(ncolw > M*(M+1)/2) {
        stop("too many columns")
    } else {
        wz = as.matrix(w)
    }
    wz
}








gaussianff = function(dispersion=0, parallel=FALSE, zero=NULL)
{
    if(!is.Numeric(dispersion, allow=1) || dispersion < 0)
        stop("bad input for argument 'dispersion'")
    estimated.dispersion = dispersion==0

    new("vglmff",
    blurb=c("Vector linear/additive model\n",
            "Links:    identity for Y1,...,YM"),
    constraints=eval(substitute(expression({
        constraints = cm.vgam(matrix(1,M,1), x, .parallel, constraints)
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero=zero ))),
    deviance=function(mu, y, w, residuals= FALSE, eta, extra=NULL) {
        M = if(is.matrix(y)) ncol(y) else 1
        n = if(is.matrix(y)) nrow(y) else length(y)
        wz = VGAM.weights.function(w=w, M=M, n=n)
        if(residuals) {
            if(M > 1) {
                U <- vchol(wz, M=M, n=n) 
                temp = mux22(U, y-mu, M=M, upper=TRUE, as.matrix=TRUE)
                dimnames(temp) = dimnames(y)
                temp
            } else (y-mu) * sqrt(wz)
        } else
            rss.vgam(y-mu, wz=wz, M=M)
    },
    initialize=eval(substitute(expression({
        if(is.R())
            assign("CQO.FastAlgorithm", TRUE, envir = VGAMenv) else
            CQO.FastAlgorithm <<- TRUE
        if(any(function.name == c("cqo","cao")) &&
           (length( .zero ) || (is.logical( .parallel ) && .parallel )))
            stop("cannot handle non-default arguments for cqo() and cao()")

        M = if(is.matrix(y)) ncol(y) else 1
        dy = dimnames(y)
        predictors.names = if(!is.null(dy[[2]])) dy[[2]] else
                           paste("Y",1:M,sep="")
        if(!length(etastart)) 
            etastart = 0 * y
    }), list( .parallel=parallel, .zero=zero ))),
    inverse=function(eta, extra=NULL) eta, 
    last=eval(substitute(expression({
        dy = dimnames(y)
        if(!is.null(dy[[2]]))
            dimnames(fit$fitted.values) = dy
        dpar = .dispersion
        if(!dpar) {
                wz = VGAM.weights.function(w=w, M=M, n=n)
                temp = rss.vgam(y-mu, wz=wz, M=M)
                dpar = temp / (length(y) - ncol(xbig.save))
        }
        misc$dispersion = dpar
        misc$default.dispersion = 0
        misc$estimated.dispersion = .estimated.dispersion
        misc$link = rep("identity", length=M)
        names(misc$link) = predictors.names

        if(is.R()) {
            if(exists("CQO.FastAlgorithm", envir = VGAMenv))
                rm("CQO.FastAlgorithm", envir = VGAMenv)
        } else {
            while(exists("CQO.FastAlgorithm"))
                remove("CQO.FastAlgorithm")
        }
    }), list( .dispersion=dispersion,
              .estimated.dispersion=estimated.dispersion ))),
    loglikelihood= function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        M = if(is.matrix(y)) ncol(y) else 1
        n = if(is.matrix(y)) nrow(y) else length(y)
        wz = VGAM.weights.function(w=w, M=M, n=n)
        temp = rss.vgam(y, wz=wz, M=M)
        -0.5 * temp
    },
    link=function(mu, extra=NULL) mu,
    vfamily="gaussianff",
    deriv=expression({
        wz = VGAM.weights.function(w=w, M=M, n=n)
        mux22(cc=t(wz), xmat=y-mu, M=M, as.mat=TRUE)
    }),
    weight= expression({
        wz
    }))
}






dposnorm = function(x, mean=0, sd=1) {
    L = max(length(x), length(mean), length(sd))
    x = rep(x, len=L); mean = rep(mean, len=L); sd = rep(sd, len=L);
    ifelse(x < 0, 0, dnorm(x=x, mean=mean, sd=sd) /
                     (1-pnorm(q=0, mean=mean, sd=sd)))
}

pposnorm = function(q, mean=0, sd=1) {
    L = max(length(q), length(mean), length(sd))
    q = rep(q, len=L); mean = rep(mean, len=L); sd = rep(sd, len=L);
    ifelse(q < 0, 0, (pnorm(q=q, mean=mean, sd=sd) -
                      pnorm(q=0, mean=mean, sd=sd)) /
                     (1-pnorm(q=0, mean=mean, sd=sd)))
}

qposnorm = function(p, mean=0, sd=1) {
    if(!is.Numeric(p, posit=TRUE) || max(p) >= 1)
        stop("bad input for argument \"p\"")
    qnorm(p=p+(1-p)*pnorm(0, mean=mean, sd=sd), mean=mean, sd=sd)
}

rposnorm = function(n, mean=0, sd=1) {
    if(!is.Numeric(n, integ=TRUE, posit=TRUE))
        stop("bad input for argument \"n\"")
    y = rnorm(n, mean=mean, sd=sd)
    mean = rep(mean, length=n)
    sd = rep(sd, length=n)
    repeat {
        index = y < 0
        if(any(index)) {
            y[index] = rnorm(n=sum(index), mean=mean[index], sd=sd[index])
        } else break
    }
    y
}

posnormal1 = function(lmean="identity", lsd="loge",
                      emean=list(), esd=list(),
                      imean=NULL, isd=NULL, zero=NULL)
{
    if(mode(lmean) != "character" && mode(lmean) != "name")
        lmean = as.character(substitute(lmean))
    if(mode(lsd) != "character" && mode(lsd) != "name")
        lsd = as.character(substitute(lsd))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(length(isd) && !is.Numeric(isd, posit=TRUE))
        stop("bad input for argument \"isd\"")
    if(!is.list(emean)) emean = list()
    if(!is.list(esd)) esd = list()

    new("vglmff",
    blurb=c("Positive (univariate) normal distribution\n\n",
            "Links:    ",
            namesof("mean", lmean, earg= emean, tag= TRUE), "; ",
            namesof("sd", lsd, earg= esd, tag= TRUE)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(y <- cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(min(y) <= 0)
            stop("response must be positive")
        predictors.names = c(namesof("mean", .lmean, earg= .emean, tag= FALSE),
                             namesof("sd",   .lsd, earg= .esd, tag= FALSE))
        if(!length(etastart)) {
            init.me = if(length( .imean)) rep( .imean, len=n) else NULL
            init.sd = if(length( .isd  )) rep( .isd  , len=n) else NULL
            if(!length(init.me)) init.me = rep(quantile(y, probs=0.40), len=n)
            if(!length(init.sd)) init.sd = rep(sd(y)*1.2, len=n)
            etastart = cbind(theta2eta(init.me, .lmean, earg= .emean),
                             theta2eta(init.sd, .lsd, earg= .esd))
        }
    }), list( .lmean=lmean, .lsd=lsd, .imean=imean, .isd=isd,
              .emean=emean, .esd=esd ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        mysd = eta2theta(eta[,2], .lsd, earg= .esd)
        mymu + mysd * dnorm(-mymu/mysd) / (1-pnorm(-mymu/mysd))
    }, list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd ))),
    last=eval(substitute(expression({
        misc$link = c("mean"= .lmean, "sd"= .lsd)  # zz mu or mean ?
        misc$earg = list("mean"= .emean, "sd"= .esd )
        misc$expected = TRUE
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        mysd = eta2theta(eta[,2], .lsd, earg= .esd)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
        if(is.R())
            sum(w*(dnorm(y, m=mymu, sd=mysd, log=TRUE) -
                   pnorm(-mymu/mysd, log=TRUE, lower.tail=FALSE))) else
            sum(w*(-log(mysd)-0.5*((y-mymu)/mysd)^2 -log(1-pnorm(-mymu/mysd))))
        }
    }, list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd ))),
    vfamily=c("posnormal1"),
    deriv=eval(substitute(expression({
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        mysd = eta2theta(eta[,2], .lsd, earg= .esd)
        zedd = (y-mymu) / mysd
        temp7 = dnorm(-mymu/mysd)
        temp8 = if(is.R()) pnorm(-mymu/mysd, low=FALSE) else 1-pnorm(-mymu/mysd)
        temp8 = temp8 * mysd
        dl.dmu = zedd / mysd - temp7 / temp8
        dl.dsd = (mymu*temp7/temp8 + zedd^2 - 1) / mysd
        dmu.deta = dtheta.deta(mymu, .lmean, earg= .emean)
        dsd.deta = dtheta.deta(mysd, .lsd, earg= .esd)
        w * cbind(dl.dmu * dmu.deta,
                  dl.dsd * dsd.deta)
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), n, dimm(M))
        ed2l.dmu2 = (1 - temp7*mymu/temp8) / mysd^2  - (temp7/temp8)^2
        ed2l.dmusd = (temp7 /(mysd * temp8)) * (1 + (mymu/mysd)^2 +
                     mymu*temp7 / temp8)
        ed2l.dsd2 = 2 / mysd^2  - (temp7 * mymu /(mysd^2 * temp8)) *
                    (1 + (mymu/mysd)^2 + mymu*temp7/temp8)
        wz[,iam(1,1,M)] = ed2l.dmu2  * dmu.deta^2
        wz[,iam(2,2,M)] = ed2l.dsd2  * dsd.deta^2
        wz[,iam(1,2,M)] = ed2l.dmusd * dsd.deta * dmu.deta
        w * wz
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd ))))
}



dbetanorm = function(x, shape1, shape2, mean=0, sd=1, log.arg=FALSE) {
    if(!is.logical(log.arg) || length(log.arg)!=1)
        stop("bad input for argument \"log.arg\"")
    ans =
    if(is.R() && log.arg) {
        dnorm(x=x, mean=mean, sd=sd, log=TRUE) +
        (shape1-1) * pnorm(q=x, mean=mean, sd=sd, log=TRUE) +
        (shape2-1) * pnorm(q=x, mean=mean, sd=sd, lower=FALSE, log=TRUE) -
        lbeta(shape1, shape2)
    } else {
        dnorm(x=x, mean=mean, sd=sd) *
        pnorm(q=x, mean=mean, sd=sd)^(shape1-1) *
        pnorm(q=x, mean=mean, sd=sd, lower=FALSE)^(shape2-1) /
        beta(shape1, shape2)
    }
    if(!is.R() && log.arg) ans = log(ans)
    ans
}

pbetanorm = function(q, shape1, shape2, mean=0, sd=1,
    lower.tail=TRUE, log.p=FALSE) {
    pbeta(q=pnorm(q=q, mean=mean, sd=sd), shape1=shape1, shape2=shape2,
          lower.tail = lower.tail, log.p = log.p)
}

qbetanorm = function(p, shape1, shape2, mean=0, sd=1) {
    if(!is.Numeric(p, posit=TRUE) || max(p) >= 1)
        stop("bad input for argument \"p\"")
    qnorm(p=qbeta(p=p, shape1=shape1, shape2=shape2), mean=mean, sd=sd)
}

rbetanorm = function(n, shape1, shape2, mean=0, sd=1) {
    if(!is.Numeric(n, integ=TRUE, posit=TRUE))
        stop("bad input for argument \"n\"")
    qnorm(p=qbeta(p=runif(n), shape1=shape1, shape2=shape2), mean=mean, sd=sd)
}



tikuv = function(d, lmean="identity", lsigma="loge",
                 emean=list(), esigma=list(),
                 isigma=NULL, zero=2)
{
    if(mode(lmean) != "character" && mode(lmean) != "name")
        lmean = as.character(substitute(lmean))
    if(mode(lsigma) != "character" && mode(lsigma) != "name")
        lsigma = as.character(substitute(lsigma))
    if(length(zero) && (!is.Numeric(zero, integer=TRUE, posit=TRUE) ||
       max(zero) > 2))
        stop("bad input for argument \"zero\"")
    if(!is.Numeric(d, allow=1) || max(d) >= 2)
        stop("bad input for argument \"d\"")
    if(!is.list(emean)) emean = list()
    if(!is.list(esigma)) esigma = list()

    new("vglmff",
    blurb=c("Short-tailed symmetric [Tiku and Vaughan (1999)] distribution\n",
            "Link:     ",
            namesof("mean", lmean, earg= emean), ", ",
            namesof("sigma", lsigma, earg= esigma),
            "\n",
            "\n",
            "Mean:     mean"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("the response must be a vector or one-column matrix")
        predictors.names = 
            c(namesof("mean", .lmean, earg= .emean, tag= FALSE),
              namesof("sigma", .lsigma, earg= .esigma, tag= FALSE))
        if(!length(etastart)) {
            sigma.init = if(length(.isigma)) rep(.isigma, length=n) else {
                hh = 2 - .d
                KK = 1 / (1 + 1/hh + 0.75/hh^2)
                K2 = 1 + 3/hh + 15/(4*hh^2)
                rep(sqrt(var(y) / (KK*K2)), len=n)
            }
            mean.init = rep(weighted.mean(y, w), len=n) 
            etastart = cbind(theta2eta(mean.init,  .lmean, earg= .emean),
                             theta2eta(sigma.init, .lsigma, earg= .esigma))
        }
    }),list( .lmean=lmean, .lsigma=lsigma, .isigma=isigma, .d=d,
             .emean=emean, .esigma=esigma ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], .lmean, earg= .emean)
    }, list( .lmean=lmean,
             .emean=emean, .esigma=esigma ))),
    last=eval(substitute(expression({
        misc$link = c("mean"= .lmean, "sigma"= .lsigma)
        misc$earg = list("mean"= .emean, "sigma"= .esigma )
        misc$expected = TRUE
        misc$d = .d 
    }), list( .lmean=lmean, .lsigma=lsigma, .d=d,
             .emean=emean, .esigma=esigma ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        sigma = eta2theta(eta[,2], .lsigma, earg= .esigma)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            zedd = (y - mymu) / sigma
            hh = 2 - .d
            sum(w * (-log(sigma) + 2 * log(1 + 0.5*zedd^2 / hh) - 0.5*zedd^2))
        }
    }, list( .lmean=lmean, .lsigma=lsigma, .d=d,
             .emean=emean, .esigma=esigma ))),
    vfamily=c("tikuv"),
    deriv=eval(substitute(expression({
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        sigma = eta2theta(eta[,2], .lsigma, earg= .esigma)
        dmu.deta = dtheta.deta(mymu, .lmean, earg= .emean)
        dsigma.deta = dtheta.deta(sigma, .lsigma, earg= .esigma)
        zedd = (y - mymu) / sigma
        hh = 2 - .d 
        gzedd = zedd / (1 + 0.5*zedd^2 / hh)
        dl.dmu = zedd / sigma - 2 * gzedd / (hh*sigma)
        dl.dsigma = (zedd^2 - 1 - 2 * zedd * gzedd / hh) / sigma
        w * cbind(dl.dmu * dmu.deta, dl.dsigma * dsigma.deta)
    }), list( .lmean=lmean, .lsigma=lsigma, .d=d,
             .emean=emean, .esigma=esigma ))),
    weight=eval(substitute(expression({
        ayy = 1 / (2*hh)
        Dnos = 1 - (2/hh) * (1 - ayy) / (1 + 2*ayy + 3*ayy^2)
        Dstar = -1 + 3 * (1 + 2*ayy + 11*ayy^2) / (1 + 2*ayy + 3*ayy^2)
        ed2l.dmymu2 = Dnos / sigma^2
        ed2l.dnu2   = Dstar / sigma^2
        wz = matrix(as.numeric(NA), n, M)  # diagonal matrix
        wz[,iam(1,1,M)] = ed2l.dmymu2 * dmu.deta^2
        wz[,iam(2,2,M)] = ed2l.dnu2 * dsigma.deta^2
        w * wz
    }), list( .lmean=lmean, .lsigma=lsigma,
             .emean=emean, .esigma=esigma ))))
}


dtikuv = function(x, d, mean=0, sigma=1) {
    if(!is.Numeric(d, allow=1) || max(d) >= 2)
        stop("bad input for argument \"d\"")
    L = max(length(x), length(mean), length(sigma))
    x = rep(x, len=L); mean = rep(mean, len=L); sigma = rep(sigma, len=L);
    hh = 2 - d
    KK = 1 / (1 + 1/hh + 0.75/hh^2)
    dnorm(x=x, mean=mean, sd=sigma) * KK * (1 + ((x-mean)/sigma)^2 / (2*hh))^2
}


ptikuv = function(q, d, mean=0, sigma=1) {
    if(!is.Numeric(d, allow=1) || max(d) >= 2)
        stop("bad input for argument \"d\"")
    L = max(length(q), length(mean), length(sigma))
    q = rep(q, len=L); mean = rep(mean, len=L); sigma = rep(sigma, len=L);
    zedd1 = 0.5 * ((q - mean) / sigma)^2
    ans = q*0 + 0.5
    hh = 2 - d
    KK = 1 / (1 + 1/hh + 0.75/hh^2)
    if(any(lhs <- q < mean)) {
        ans[lhs] = ( KK/(2*sqrt(pi))) * (
        gamma(0.5) * (1 - pgamma(zedd1[lhs], 0.5)) +
        2 * gamma(1.5) * (1 - pgamma(zedd1[lhs], 1.5)) / hh +
        gamma(2.5) * (1 - pgamma(zedd1[lhs], 2.5)) / hh^2)
    }
    if(any(rhs <- q > mean)) {
        ans[rhs] = 1.0 - Recall(q=(2*mean[rhs]-q[rhs]), d=d,
                   mean=mean[rhs], sigma=sigma[rhs])
    }
    ans
}


qtikuv = function(p, d, mean=0, sigma=1, ...) {
    if(!is.Numeric(p, posit=TRUE) || max(p) >= 1)
        stop("bad input for argument \"p\"")
    if(!is.Numeric(d, allow=1) || max(d) >= 2)
        stop("bad input for argument \"d\"")
    if(!is.Numeric(mean))
        stop("bad input for argument \"mean\"")
    if(!is.Numeric(sigma))
        stop("bad input for argument \"sigma\"")
    L = max(length(p), length(mean), length(sigma))
    p = rep(p, len=L); mean = rep(mean, len=L); sigma = rep(sigma, len=L);
    ans = rep(0.0, len=L)
    myfun = function(x, d, mean=0, sigma=1, p)
        ptikuv(q=x, d=d, mean=mean, sigma=sigma) - p
    for(i in 1:L) {
        Lower = ifelse(p[i] <= 0.5, mean[i] - 3 * sigma[i], mean[i])
        while(ptikuv(q=Lower, d=d, mean=mean[i], sigma=sigma[i]) > p[i])
            Lower = Lower - sigma[i]
        Upper = ifelse(p[i] >= 0.5, mean[i] + 3 * sigma[i], mean[i])
        while(ptikuv(q=Upper, d=d, mean=mean[i], sigma=sigma[i]) < p[i])
            Upper = Upper + sigma[i]
        ans[i] = uniroot(f=myfun, lower=Lower, upper=Upper,
                         d=d, mean=mean[i], sigma=sigma[i], p=p[i], ...)$root
    }
    ans
}


rtikuv = function(n, d, mean=0, sigma=1, Smallno=1.0e-6) {
    if(!is.Numeric(n, posit=TRUE, integ=TRUE))
        stop("bad input for argument \"n\"")
    if(!is.Numeric(d, allow=1) || max(d) >= 2)
        stop("bad input for argument \"d\"")
    if(!is.Numeric(mean, allow=1))
        stop("bad input for argument \"mean\"")
    if(!is.Numeric(sigma, allow=1))
        stop("bad input for argument \"sigma\"")
    if(!is.Numeric(Smallno, posit=TRUE, allow=1) || Smallno > 0.01 ||
       Smallno < 2 * .Machine$double.eps)
        stop("bad input for argument \"Smallno\"")
    ans = rep(0.0, len=n)

    ptr1 = 1; ptr2 = 0
    hh = 2 - d
    KK = 1 / (1 + 1/hh + 0.75/hh^2)
    ymax = ifelse(hh < 2,
                  dtikuv(x=mean + sigma*sqrt(4 - 2*hh), d=d, m=mean, s=sigma),
                  KK / (sqrt(2 * pi) * sigma))
    while(ptr2 < n) {
        Lower = mean - 5 * sigma
        while(ptikuv(q=Lower, d=d, mean=mean, sigma=sigma) > Smallno)
            Lower = Lower - sigma
        Upper = mean + 5 * sigma
        while(ptikuv(q=Upper, d=d, mean=mean, sigma=sigma) < 1-Smallno)
            Upper = Upper + sigma
        x = runif(2*n, min=Lower, max=Upper)
        index = runif(2*n, max=ymax) < dtikuv(x,d=d,m=mean,s=sigma)
        sindex = sum(index)
        if(sindex) {
            ptr2 = min(n, ptr1 + sindex - 1)
            ans[ptr1:ptr2] = (x[index])[1:(1+ptr2-ptr1)]
            ptr1 = ptr2 + 1
        }
    }
    ans
}




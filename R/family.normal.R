# These functions are
# Copyright (C) 1998-2010 T.W. Yee, University of Auckland. All rights reserved.




VGAM.weights.function = function(w, M, n) {
    ncolw = ncol(as.matrix(w))
    if (ncolw == 1) {
        wz = matrix(w, nrow=n, ncol=M) # w_i * diag(M)
    } else if (ncolw == M) {
        wz = as.matrix(w)
    } else if (ncolw < M && M > 1) {
        stop("ambiguous input for weights")
    } else if (ncolw > M*(M+1)/2) {
        stop("too many columns")
    } else {
        wz = as.matrix(w)
    }
    wz
}








 gaussianff = function(dispersion=0, parallel=FALSE, zero=NULL)
{
    if (!is.Numeric(dispersion, allow=1) || dispersion < 0)
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
        M = if (is.matrix(y)) ncol(y) else 1
        n = if (is.matrix(y)) nrow(y) else length(y)
        wz = VGAM.weights.function(w=w, M=M, n=n)
        if (residuals) {
            if (M > 1) {
                U <- vchol(wz, M=M, n=n) 
                temp = mux22(U, y-mu, M=M, upper=TRUE, as.matrix=TRUE)
                dimnames(temp) = dimnames(y)
                temp
            } else (y-mu) * sqrt(wz)
        } else
            rss.vgam(y-mu, wz=wz, M=M)
    },
    initialize=eval(substitute(expression({
        if (is.R())
            assign("CQO.FastAlgorithm", TRUE, envir = VGAMenv) else
            CQO.FastAlgorithm <<- TRUE
        if (any(function.name == c("cqo","cao")) &&
           (length( .zero ) || (is.logical( .parallel ) && .parallel )))
            stop("cannot handle non-default arguments for cqo() and cao()")

        M = if (is.matrix(y)) ncol(y) else 1
        dy = dimnames(y)
        predictors.names = if (!is.null(dy[[2]])) dy[[2]] else
                           paste("Y",1:M,sep="")
        if (!length(etastart)) 
            etastart = 0 * y
    }), list( .parallel=parallel, .zero=zero ))),
    inverse=function(eta, extra=NULL) eta, 
    last=eval(substitute(expression({
        dy = dimnames(y)
        if (!is.null(dy[[2]]))
            dimnames(fit$fitted.values) = dy
        dpar = .dispersion
        if (!dpar) {
                wz = VGAM.weights.function(w=w, M=M, n=n)
                temp = rss.vgam(y-mu, wz=wz, M=M)
                dpar = temp / (length(y) -
                       (if(is.numeric(ncol(X_vlm_save))) ncol(X_vlm_save) else 0))
        }
        misc$dispersion = dpar
        misc$default.dispersion = 0
        misc$estimated.dispersion = .estimated.dispersion
        misc$link = rep("identity", length=M)
        names(misc$link) = predictors.names

        if (is.R()) {
            if (exists("CQO.FastAlgorithm", envir = VGAMenv))
                rm("CQO.FastAlgorithm", envir = VGAMenv)
        } else {
            while(exists("CQO.FastAlgorithm"))
                remove("CQO.FastAlgorithm")
        }
    }), list( .dispersion=dispersion,
              .estimated.dispersion=estimated.dispersion ))),
    loglikelihood= function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        M = if (is.matrix(y)) ncol(y) else 1
        n = if (is.matrix(y)) nrow(y) else length(y)
        wz = VGAM.weights.function(w=w, M=M, n=n)
        temp = rss.vgam(y-mu, wz=wz, M=M)
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







dposnorm = function(x, mean=0, sd=1, log=FALSE) {
    log.arg = log
    rm(log)
    if (!is.logical(log.arg) || length(log.arg)!=1)
        stop("bad input for argument 'log'")
    L = max(length(x), length(mean), length(sd))
    x = rep(x, len=L); mean = rep(mean, len=L); sd = rep(sd, len=L);

    if (log.arg) {
        ifelse(x < 0, log(0), dnorm(x, m=mean, sd=sd, log=TRUE) -
               pnorm(mean/sd, log=TRUE))
    } else {
        ifelse(x < 0, 0, dnorm(x=x, me=mean, sd=sd) / pnorm(mean/sd))
    }
}

pposnorm = function(q, mean=0, sd=1) {
    L = max(length(q), length(mean), length(sd))
    q = rep(q, len=L); mean = rep(mean, len=L); sd = rep(sd, len=L);
    ifelse(q < 0, 0, (pnorm(q=q, mean=mean, sd=sd) -
                      pnorm(q=0, mean=mean, sd=sd)) / pnorm(q= mean/sd))
}

qposnorm = function(p, mean=0, sd=1) {
    if (!is.Numeric(p, posit=TRUE) || max(p) >= 1)
        stop("bad input for argument 'p'")
    qnorm(p=p+(1-p)*pnorm(0, mean=mean, sd=sd), mean=mean, sd=sd)
}

rposnorm = function(n, mean=0, sd=1) {
    if (!is.Numeric(n, integ=TRUE, posit=TRUE))
        stop("bad input for argument 'n'")
    mean = rep(mean, length=n)
    sd = rep(sd, length=n)
    qnorm(p=runif(n, min=pnorm(0, m=mean, sd=sd)), m=mean, sd=sd)
}



 posnormal1.control <- function(save.weight=TRUE, ...) {
    list(save.weight=save.weight)
}




 posnormal1 = function(lmean="identity", lsd="loge",
                       emean=list(), esd=list(),
                       imean=NULL, isd=NULL,
                       nsimEIM=100, zero=NULL)
{
 warning("this VGAM family function is not working properly yet")

    if (mode(lmean) != "character" && mode(lmean) != "name")
        lmean = as.character(substitute(lmean))
    if (mode(lsd) != "character" && mode(lsd) != "name")
        lsd = as.character(substitute(lsd))
    if (length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument 'zero'")
    if (length(isd) && !is.Numeric(isd, posit=TRUE))
        stop("bad input for argument 'isd'")
    if (!is.list(emean)) emean = list()
    if (!is.list(esd)) esd = list()
    if (length(nsimEIM))
        if (!is.Numeric(nsimEIM, allow=1, integ=TRUE) || nsimEIM <= 10)
            stop("'nsimEIM' should be an integer greater than 10")

    new("vglmff",
    blurb=c("Positive (univariate) normal distribution\n\n",
            "Links:    ",
            namesof("mean", lmean, earg= emean, tag= TRUE), "; ",
            namesof("sd", lsd, earg= esd, tag= TRUE)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if (ncol(y <- cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if (min(y) <= 0)
            stop("response must be positive")
        predictors.names = c(namesof("mean", .lmean, earg= .emean, tag= FALSE),
                             namesof("sd",   .lsd, earg= .esd, tag= FALSE))
        if (!length(etastart)) {
            init.me = if (length( .imean)) rep( .imean, len=n) else NULL
            init.sd = if (length( .isd  )) rep( .isd  , len=n) else NULL
            if (!length(init.me)) init.me = rep(quantile(y, probs=0.40), len=n)
            if (!length(init.sd)) init.sd = rep(sd(y)*1.2, len=n)
            etastart = cbind(theta2eta(init.me, .lmean, earg= .emean),
                             theta2eta(init.sd, .lsd, earg= .esd))
        }
    }), list( .lmean=lmean, .lsd=lsd, .imean=imean, .isd=isd,
              .emean=emean, .esd=esd ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        mysd = eta2theta(eta[,2], .lsd, earg= .esd)
        mymu + mysd * dnorm(-mymu/mysd) / pnorm(mymu/mysd)
    }, list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd ))),
    last=eval(substitute(expression({
        misc$link = c("mean"= .lmean, "sd"= .lsd)
        misc$earg = list("mean"= .emean, "sd"= .esd )
        misc$expected = TRUE
        misc$nsimEIM = .nsimEIM
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd,
              .nsimEIM=nsimEIM ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        mysd = eta2theta(eta[,2], .lsd, earg= .esd)
        if (residuals) stop("loglikelihood residuals not implemented yet") else {

            sum(w * dposnorm(x=y, m=mymu, sd=mysd, log=TRUE))
        }
    }, list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd ))),
    vfamily=c("posnormal1"),
    deriv=eval(substitute(expression({
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        mysd = eta2theta(eta[,2], .lsd, earg= .esd)
        zedd = (y-mymu) / mysd
        temp7 = dnorm(-mymu/mysd)
        temp8 = pnorm(mymu/mysd) * mysd
        dl.dmu = zedd / mysd^2 - temp7 / temp8
        dl.dsd = (mymu*temp7/temp8 + zedd^3 / mysd - 1) / mysd
        dmu.deta = dtheta.deta(mymu, .lmean, earg= .emean)
        dsd.deta = dtheta.deta(mysd, .lsd, earg= .esd)
        dthetas.detas = cbind(dmu.deta, dsd.deta)
        w * dthetas.detas * cbind(dl.dmu, dl.dsd)
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd ))),
    weight=eval(substitute(expression({
        run.varcov = 0
        ind1 = iam(NA, NA, M=M, both=TRUE, diag=TRUE)
        if (length( .nsimEIM )) {
            for(ii in 1:( .nsimEIM )) {
                ysim <- rposnorm(n, m=mymu, sd=mysd)
                zedd = (ysim-mymu) / mysd
                temp7 = dnorm(-mymu/mysd)
                temp8 = pnorm(mymu/mysd) * mysd
                dl.dmu = zedd / mysd^2 - temp7 / temp8
                dl.dsd = (mymu*temp7/temp8 + zedd^3 / mysd - 1) / mysd

                rm(ysim)
                temp3 = matrix(c(dl.dmu, dl.dsd), n, 2)
                run.varcov = ((ii-1) * run.varcov +
                           temp3[,ind1$row.index]*temp3[,ind1$col.index]) / ii
            }
            wz = if (intercept.only)
                matrix(colMeans(run.varcov),
                       n, ncol(run.varcov), byrow=TRUE) else run.varcov

            wz = wz * dthetas.detas[,ind1$row] * dthetas.detas[,ind1$col]
            wz = w * matrix(wz, n, dimm(M))
        } else {
            wz = matrix(as.numeric(NA), n, dimm(M))
            ed2l.dmu2 = (1 - temp7*mymu/temp8) / mysd^2  - (temp7/temp8)^2
            ed2l.dmusd = (temp7 /(mysd * temp8)) * (1 + (mymu/mysd)^2 +
                         mymu*temp7 / temp8)
            ed2l.dsd2 = 2 / mysd^2  - (temp7 * mymu /(mysd^2 * temp8)) *
                        (1 + (mymu/mysd)^2 + mymu*temp7/temp8)
            wz[,iam(1,1,M)] = ed2l.dmu2  * dmu.deta^2
            wz[,iam(2,2,M)] = ed2l.dsd2  * dsd.deta^2
            wz[,iam(1,2,M)] = ed2l.dmusd * dsd.deta * dmu.deta
            wz = w * wz
        }
        wz
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd,
              .nsimEIM=nsimEIM ))))
}




dbetanorm = function(x, shape1, shape2, mean=0, sd=1, log=FALSE) {
    log.arg = log
    rm(log)
    if (!is.logical(log.arg) || length(log.arg)!=1)
        stop("bad input for argument 'log'")
    ans =
    if (is.R() && log.arg) {
        dnorm(x=x, mean=mean, sd=sd, log=TRUE) +
        (shape1-1) * pnorm(q=x, mean=mean, sd=sd, log=TRUE) +
        (shape2-1) * pnorm(q=x, mean=mean, sd=sd, lower=FALSE, log=TRUE) -
        lbeta(shape1, shape2)
    } else {
        dnorm(x=x, mean=mean, sd=sd) *
        pnorm(q=x, mean=mean, sd=sd)^(shape1-1) *
    pnorm(q=x, mean=mean, sd=sd, lower=FALSE)^(shape2-1) / beta(shape1, shape2)
    }
    if (!is.R() && log.arg) ans = log(ans)
    ans
}


pbetanorm = function(q, shape1, shape2, mean=0, sd=1,
    lower.tail=TRUE, log.p=FALSE) {
    pbeta(q=pnorm(q=q, mean=mean, sd=sd), shape1=shape1, shape2=shape2,
          lower.tail = lower.tail, log.p = log.p)
}


qbetanorm = function(p, shape1, shape2, mean=0, sd=1) {
    if (!is.Numeric(p, posit=TRUE) || max(p) >= 1)
        stop("bad input for argument 'p'")
    qnorm(p=qbeta(p=p, shape1=shape1, shape2=shape2), mean=mean, sd=sd)
}


rbetanorm = function(n, shape1, shape2, mean=0, sd=1) {
    if (!is.Numeric(n, integ=TRUE, posit=TRUE))
        stop("bad input for argument 'n'")
    qnorm(p=qbeta(p=runif(n), shape1=shape1, shape2=shape2), mean=mean, sd=sd)
}




dtikuv = function(x, d, mean=0, sigma=1, log = FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    if (!is.Numeric(d, allow=1) || max(d) >= 2)
        stop("bad input for argument 'd'")
    L = max(length(x), length(mean), length(sigma))
    x = rep(x, len=L); mean = rep(mean, len=L); sigma = rep(sigma, len=L);
    hh = 2 - d
    KK = 1 / (1 + 1/hh + 0.75/hh^2)
    if (log.arg) {
        dnorm(x=x, mean=mean, sd=sigma, log=TRUE) + log(KK) +
        2 * log1p(((x-mean)/sigma)^2 / (2*hh))
    } else {
        dnorm(x=x, mean=mean, sd=sigma) * KK *
        (1 + ((x-mean)/sigma)^2 / (2*hh))^2
    }
}


ptikuv = function(q, d, mean=0, sigma=1) {
    if (!is.Numeric(d, allow=1) || max(d) >= 2)
        stop("bad input for argument 'd'")
    L = max(length(q), length(mean), length(sigma))
    q = rep(q, len=L); mean = rep(mean, len=L); sigma = rep(sigma, len=L);
    zedd1 = 0.5 * ((q - mean) / sigma)^2
    ans = q*0 + 0.5
    hh = 2 - d
    KK = 1 / (1 + 1/hh + 0.75/hh^2)
    if (any(lhs <- q < mean)) {
        ans[lhs] = ( KK/(2*sqrt(pi))) * (
        gamma(0.5) * (1 - pgamma(zedd1[lhs], 0.5)) +
        2 * gamma(1.5) * (1 - pgamma(zedd1[lhs], 1.5)) / hh +
        gamma(2.5) * (1 - pgamma(zedd1[lhs], 2.5)) / hh^2)
    }
    if (any(rhs <- q > mean)) {
        ans[rhs] = 1.0 - Recall(q=(2*mean[rhs]-q[rhs]), d=d,
                   mean=mean[rhs], sigma=sigma[rhs])
    }
    ans
}


qtikuv = function(p, d, mean=0, sigma=1, ...) {
    if (!is.Numeric(p, posit=TRUE) || max(p) >= 1)
        stop("bad input for argument 'p'")
    if (!is.Numeric(d, allow=1) || max(d) >= 2)
        stop("bad input for argument 'd'")
    if (!is.Numeric(mean))
        stop("bad input for argument 'mean'")
    if (!is.Numeric(sigma))
        stop("bad input for argument 'sigma'")
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
    if (!is.Numeric(n, posit=TRUE, integ=TRUE))
        stop("bad input for argument 'n'")
    if (!is.Numeric(d, allow=1) || max(d) >= 2)
        stop("bad input for argument 'd'")
    if (!is.Numeric(mean, allow=1))
        stop("bad input for argument 'mean'")
    if (!is.Numeric(sigma, allow=1))
        stop("bad input for argument 'sigma'")
    if (!is.Numeric(Smallno, posit=TRUE, allow=1) || Smallno > 0.01 ||
       Smallno < 2 * .Machine$double.eps)
        stop("bad input for argument 'Smallno'")
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
        if (sindex) {
            ptr2 = min(n, ptr1 + sindex - 1)
            ans[ptr1:ptr2] = (x[index])[1:(1+ptr2-ptr1)]
            ptr1 = ptr2 + 1
        }
    }
    ans
}




 tikuv = function(d, lmean="identity", lsigma="loge",
                  emean=list(), esigma=list(),
                  isigma=NULL, zero=2)
{
    if (mode(lmean) != "character" && mode(lmean) != "name")
        lmean = as.character(substitute(lmean))
    if (mode(lsigma) != "character" && mode(lsigma) != "name")
        lsigma = as.character(substitute(lsigma))
    if (length(zero) && (!is.Numeric(zero, integer=TRUE, posit=TRUE) ||
       max(zero) > 2))
        stop("bad input for argument 'zero'")
    if (!is.Numeric(d, allow=1) || max(d) >= 2)
        stop("bad input for argument 'd'")
    if (!is.list(emean)) emean = list()
    if (!is.list(esigma)) esigma = list()

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
        if (ncol(cbind(y)) != 1)
            stop("the response must be a vector or one-column matrix")
        predictors.names = 
            c(namesof("mean", .lmean, earg= .emean, tag= FALSE),
              namesof("sigma", .lsigma, earg= .esigma, tag= FALSE))
        if (!length(etastart)) {
            sigma.init = if (length(.isigma)) rep(.isigma, length=n) else {
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
        if (residuals) stop("loglikelihood residuals not implemented yet") else {
            sum(w * dtikuv(x=y, d= .d, mean=mymu, sigma=sigma, log = TRUE))
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



dfnorm = function(x, mean=0, sd=1, a1=1, a2=1) {
    if (!is.Numeric(a1, posit=TRUE) || !is.Numeric(a2, posit=TRUE))
        stop("bad input for arguments 'a1' and 'a2'")
    if (any(a1 <= 0 | a2 <= 0))
        stop("arguments 'a1' and 'a2' must have positive values only")
    ans = dnorm(x=x/(a1*sd) - mean/sd)/(a1*sd) +
          dnorm(x=x/(a2*sd) + mean/sd)/(a2*sd)
    ans[x < 0] = 0
    ans[a1 <= 0 | a2 <= 0 | is.na(a1) | is.na(a2)] = NA
    ans
}

pfnorm = function(q, mean=0, sd=1, a1=1, a2=1) {
    if (!is.Numeric(a1, posit=TRUE) || !is.Numeric(a2, posit=TRUE))
        stop("bad input for arguments 'a1' and 'a2'")
    if (any(a1 <= 0 | a2 <= 0))
        stop("arguments 'a1' and 'a2' must have positive values only")
    L = max(length(q), length(mean), length(sd))
    q = rep(q, len=L); mean = rep(mean, len=L); sd = rep(sd, len=L);
    ifelse(q < 0, 0, pnorm(q=q/(a1*sd) - mean/sd) - pnorm(q=-q/(a2*sd) - mean/sd))
}

qfnorm = function(p, mean=0, sd=1, a1=1, a2=1, ...) {
    if (!is.Numeric(p, posit=TRUE) || max(p) >= 1)
        stop("bad input for argument 'p'")
    if (!is.Numeric(a1, posit=TRUE) || !is.Numeric(a2, posit=TRUE))
        stop("bad input for arguments 'a1' and 'a2'")
    if (any(a1 <= 0 | a2 <= 0))
        stop("arguments 'a1' and 'a2' must have positive values only")

    L = max(length(p), length(mean), length(sd), length(a1), length(a2))
    p = rep(p, len=L); mean = rep(mean, len=L); sd = rep(sd, len=L);
    a1 = rep(a1, len=L); a2 = rep(a2, len=L);
    ans = rep(0.0, len=L)
    myfun = function(x, mean=0, sd=1, a1=1, a2=2, p)
        pfnorm(q=x, mean=mean, sd=sd, a1=a1, a2=a2) - p
    for(i in 1:L) {
        mytheta = mean[i]/sd[i]
        EY = sd[i] * ((a1[i]+a2[i]) * (mytheta * pnorm(mytheta) + dnorm(mytheta)) -
             a2[i] * mytheta)
        Upper = 2 * EY
        while(pfnorm(q=Upper, mean=mean[i], sd=sd[i], a1=a1[i], a2=a2[i]) < p[i])
            Upper = Upper + sd[i]
        ans[i] = uniroot(f=myfun, lower=0, upper=Upper,
                         mean=mean[i], sd=sd[i], a1=a1[i], a2=a2[i], p=p[i], ...)$root
    }
    ans
}

rfnorm = function(n, mean=0, sd=1, a1=1, a2=1) {
    if (!is.Numeric(n, integ=TRUE, posit=TRUE))
        stop("bad input for argument 'n'")
    if (!is.Numeric(a1, posit=TRUE) || !is.Numeric(a2, posit=TRUE))
        stop("bad input for arguments 'a1' and 'a2'")
    if (any(a1 <= 0 | a2 <= 0))
        stop("arguments 'a1' and 'a2' must have positive values only")
    X = rnorm(n, mean=mean, sd=sd)
    pmax(a1 * X, -a2*X)
}


 fnormal1 =  function(lmean="identity", lsd="loge", emean=list(), esd=list(),
                      imean=NULL, isd=NULL, a1=1, a2=1, nsimEIM=500,
                      method.init=1, zero=NULL)
{
    if (!is.Numeric(a1, posit=TRUE, allow=1) ||
       !is.Numeric(a2, posit=TRUE, allow=1))
        stop("bad input for arguments 'a1' and 'a2'")
    if (any(a1 <= 0 | a2 <= 0))
        stop("arguments 'a1' and 'a2' must each be a positive value")
    if (!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 2)
        stop("'method.init' must be 1 or 2")

    if (mode(lmean) != "character" && mode(lmean) != "name")
        lmean = as.character(substitute(lmean))
    if (mode(lsd) != "character" && mode(lsd) != "name")
        lsd = as.character(substitute(lsd))
    if (length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument 'zero'")
    if (!is.list(emean)) emean = list()
    if (!is.list(esd)) esd = list()
    if (!is.Numeric(nsimEIM, allow=1, integ=TRUE) || nsimEIM <= 10)
        stop("'nsimEIM' should be an integer greater than 10")
    if (length(imean) && !is.Numeric(imean))
        stop("bad input for 'imean'")
    if (length(isd) && !is.Numeric(isd, posit=TRUE))
        stop("bad input for 'isd'")

    new("vglmff",
    blurb=c("(Generalized) folded univariate normal distribution\n\n",
            "Link:     ",
            namesof("mean", lmean, earg=emean, tag= TRUE), "; ",
            namesof("sd", lsd, earg=esd, tag= TRUE)),
    initialize=eval(substitute(expression({
        predictors.names = c(namesof("mean", .lmean, earg=.emean, tag=FALSE),
                             namesof("sd",   .lsd, earg=.esd, tag=FALSE))
        if ((ncol(y <- cbind(y)) != 1) || any(y <= 0))
 stop("response must be a vector or a one-column matrix with positive values")
        if (!length(etastart)) {
            junk = if (is.R()) lm.wfit(x=x, y=y, w=w) else
                              lm.wfit(x=x, y=y, w=w, method="qr")


 if (FALSE) {
        if ((ncol(cbind(w)) != 1) || any(w != round(w)))
   stop("'weights' must be a vector or a one-column matrix with integer values")
            m1d = meany = weighted.mean(y, w)
            m2d = weighted.mean(y^2, w)
            stddev = sqrt( sum(w * junk$resid^2) / junk$df.residual )
            Ahat = m1d^2 / m2d
            thetahat = sqrt(max(1/Ahat -1, 0.1))
            mean.init = rep(if(length( .imean)) .imean else
                thetahat * sqrt((stddev^2 + meany^2) * Ahat), len=n)
            sd.init = rep(if(length( .isd)) .isd else
                sqrt((stddev^2 + meany^2) * Ahat), len=n)
}


            stddev = sqrt( sum(w * junk$resid^2) / junk$df.residual )
            meany = weighted.mean(y, w)
            mean.init = rep(if(length( .imean)) .imean else
                {if( .method.init == 1) median(y) else meany}, len=n)
            sd.init = rep(if(length( .isd)) .isd else
                {if( .method.init == 1)  stddev else 1.2*sd(y)}, len=n)
            etastart = cbind(theta2eta(mean.init, .lmean, earg= .emean),
                             theta2eta(sd.init, .lsd, earg= .esd))
        }
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd,
              .imean=imean, .isd=isd, .a1=a1, .a2=a2,
              .method.init=method.init ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        mysd = eta2theta(eta[,2], .lsd, earg= .esd)
        mytheta = mymu/mysd
        mysd * (( .a1+ .a2) * (mytheta * pnorm(mytheta) +
                dnorm(mytheta)) - .a2 * mytheta)
    }, list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd, .a1=a1, .a2=a2 ))),
    last=eval(substitute(expression({
        misc$link = c("mu"= .lmean, "sd"= .lsd)
        misc$earg = list("mu"= .emean, "sd"= .esd)
        misc$expected = TRUE
        misc$nsimEIM = .nsimEIM
        misc$simEIM = TRUE
        misc$method.init = .method.init
        misc$a1 = .a1
        misc$a2 = .a2
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd,
              .method.init=method.init, .nsimEIM=nsimEIM, .a1=a1, .a2=a2 ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        mysd = eta2theta(eta[,2], .lsd, earg= .esd)
        a1vec = .a1
        a2vec = .a2
        if (residuals) stop("loglikelihood residuals not implemented yet") else {
            sum(w*log(dnorm(x=y/(a1vec*mysd) - mymu/mysd)/(a1vec*mysd) +
                      dnorm(x=y/(a2vec*mysd) + mymu/mysd)/(a2vec*mysd)))
        }
    }, list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd, .a1=a1, .a2=a2 ))),
    vfamily=c("fnormal1"),
    deriv=eval(substitute(expression({
        mymu = eta2theta(eta[,1], .lmean, earg= .emean)
        mysd = eta2theta(eta[,2], .lsd, earg= .esd)
        dmu.deta = dtheta.deta(mymu, .lmean, earg= .emean)
        dsd.deta = dtheta.deta(mysd, .lsd, earg= .esd)
        a1vec = .a1
        a2vec = .a2
        d3 = deriv3(~ log((exp(-0.5*(y/(a1vec*mysd) - mymu/mysd)^2)/a1vec +
                           exp(-0.5*(y/(a2vec*mysd) +
                               mymu/mysd)^2)/a2vec)/(mysd*sqrt(2*pi))),
                    name=c("mymu","mysd"), hessian= FALSE)
        eval.d3 = eval(d3)
        dl.dthetas =  attr(eval.d3, "gradient")  # == cbind(dl.dmu, dl.dsd)
        dtheta.detas = cbind(dmu.deta, dsd.deta)
        w * dtheta.detas * dl.dthetas
    }), list( .lmean=lmean, .lsd=lsd, .emean=emean, .esd=esd, .a1=a1, .a2=a2 ))),
    weight=eval(substitute(expression({
        de3 = deriv3(~ log((exp(-0.5*(ysim/(a1vec*mysd) - mymu/mysd)^2)/a1vec +
                            exp(-0.5*(ysim/(a2vec*mysd) + mymu/mysd)^2)/a2vec)/(mysd*sqrt(2*pi))),
                     name=c("mymu","mysd"), hessian= TRUE)
        run.mean = 0
        for(ii in 1:( .nsimEIM )) {
            ysim = rfnorm(n=n, mean=mymu, sd=mysd, a1= a1vec, a2= a2vec)
            eval.de3 = eval(de3)
            d2l.dthetas2 =  attr(eval.de3, "hessian")
            rm(ysim)

            temp3 = matrix(0, n, dimm(M))
            for(ss in 1:M)
                for(tt in ss:M)
                    temp3[,iam(ss,tt,M)] =  -d2l.dthetas2[,ss,tt]

            run.mean = ((ii-1) * run.mean + temp3) / ii
        }

        wz = if (intercept.only)
            matrix(colMeans(run.mean), n, dimm(M), byrow=TRUE) else run.mean

        index0 = iam(NA, NA, M=M, both=TRUE, diag=TRUE)
        wz = wz * dtheta.detas[,index0$row] * dtheta.detas[,index0$col]
        w * wz
    }), list( .nsimEIM=nsimEIM, .a1=a1, .a2=a2 ))))
}





lqnorm.control = function(trace=TRUE, ...)
{
    list(trace=trace)
}


lqnorm = function(qpower=2, link="identity", earg=list(),
                  method.init=1, imu=NULL, shrinkage.init=0.95)
{
    if (mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if (!is.list(earg)) eerg = list()
    if (!is.Numeric(qpower, allow=1) || qpower <= 1)
        stop("bad input for argument 'qpower'")
    if (!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 3)
        stop("'method.init' must be 1 or 2 or 3")
    if (!is.Numeric(shrinkage.init, allow=1) || shrinkage.init < 0 ||
       shrinkage.init > 1) stop("bad input for argument 'shrinkage.init'")

    new("vglmff",
    blurb=c("Minimizing the q-norm of residuals\n",
            "Links:    ",
            namesof("Y1", link, earg=earg, tag= TRUE)),
    initialize=eval(substitute(expression({
        M = if (is.matrix(y)) ncol(y) else 1
        if (M != 1)
            stop("response must be a vector or a one-column matrix")
        dy = dimnames(y)
        predictors.names = if (!is.null(dy[[2]])) dy[[2]] else
                           paste("mu", 1:M, sep="")
        predictors.names = namesof(predictors.names, link= .link,
                                   earg= .earg, short=TRUE)
        if (!length(etastart))  {
            meany = weighted.mean(y, w)
            mean.init = rep(if(length( .imu)) .imu else
                {if( .method.init == 2) median(y) else 
                 if ( .method.init == 1) meany else
                 .sinit * meany + (1 - .sinit) * y
                }, len=n)
            etastart = theta2eta(mean.init, link= .link, earg= .earg)
        }
    }), list( .method.init=method.init, .imu=imu, .sinit=shrinkage.init,
              .link=link, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mu = eta2theta(eta, link= .link, earg= .earg)
        mu
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        dy = dimnames(y)
        if (!is.null(dy[[2]]))
            dimnames(fit$fitted.values) = dy
        misc$link = rep( .link, length=M)
        names(misc$link) = predictors.names
        misc$earg = list(mu = .earg)
        misc$qpower = .qpower
        misc$method.init = .method.init
        misc$objectiveFunction = sum( w * (abs(y - mu))^(.qpower) )
    }), list( .qpower=qpower,
              .link=link, .earg=earg,
              .method.init=method.init ))),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, link= .link, earg=.earg)
    }, list( .link=link, .earg=earg ))),
    vfamily="lqnorm",
    deriv=eval(substitute(expression({
        dmu.deta = dtheta.deta(theta=mu, link=.link, earg= .earg )
        myresid = y - mu
        signresid = sign(myresid)
        temp2 = (abs(myresid))^(.qpower-1)
        .qpower * w * temp2 * signresid * dmu.deta
    }), list( .qpower=qpower, .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        temp3 = (abs(myresid))^(.qpower-2)
        wz = .qpower * (.qpower - 1) * w * temp3 * dmu.deta^2
        wz
    }), list( .qpower=qpower, .link=link, .earg=earg ))))
}






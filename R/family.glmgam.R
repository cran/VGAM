# These functions are
# Copyright (C) 1998-2010 T.W. Yee, University of Auckland. All rights reserved.









 binomialff = function(link="logit", earg=list(),
                       dispersion=1, mv=FALSE, onedpar=!mv,
                       parallel = FALSE, zero=NULL)

{


    estimated.dispersion <- dispersion==0
    if (mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb=if(mv) c("Multivariate Binomial model\n\n", 
           "Link:     ", namesof("mu[,j]", link, earg= earg), "\n",
           "Variance: mu[,j]*(1-mu[,j])") else
           c("Binomial model\n\n", 
           "Link:     ", namesof("mu", link, earg= earg), "\n",
           "Variance: mu*(1-mu)"),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(1,M,1), x, .parallel, constraints)
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero=zero ))),
    deviance=function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        Deviance.categorical.data.vgam(mu=cbind(mu,1-mu), y=cbind(y,1-y),
                                       w=w, residuals = residuals,
                                       eta=eta, extra=extra)
    },
    initialize=eval(substitute(expression({
        if (is.R()) {
            assign("CQO.FastAlgorithm", ( .link=="logit" || .link=="cloglog"),
                   envir=VGAMenv)
            assign("modelno", if ( .link=="logit") 1 else
                   if ( .link=="cloglog") 4 else NULL, envir=VGAMenv)
        }  else {
            CQO.FastAlgorithm <<- ( .link == "logit" || .link=="cloglog")
          modelno <<- if ( .link=="logit") 1 else if ( .link=="cloglog") 4 else NULL
        }
        if (.mv) {
            y = as.matrix(y)
            M = ncol(y)
            if (!all(y == 0 | y == 1))
                stop("response must contain 0's and 1's only")
            dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
            dn2 = if (length(dn2)) {
                paste("E[", dn2, "]", sep="") 
            } else {
                paste("mu", 1:M, sep="") 
            }
            predictors.names = namesof(if(M>1) dn2 else
                "mu", .link, earg= .earg, short=TRUE)

            mustart = (0.5 + w * y) / (1 + w)
        } else {
            NCOL = function (x) 
                if (is.array(x) && length(dim(x)) > 1 ||
                is.data.frame(x)) ncol(x) else as.integer(1)

            if (NCOL(y) == 1) {
                if (is.factor(y)) y = y != levels(y)[1]
                nn = rep(1, n)
                if (!all(y >= 0 & y <= 1))
                    stop("response values must be in [0, 1]")
                mustart = (0.5 + w * y) / (1 + w)
                no.successes = w * y
                if (any(abs(no.successes - round(no.successes)) > 0.001))
                    stop("Number of successes must be integer-valued")
            } else if (NCOL(y) == 2) {
                if (any(abs(y - round(y)) > 0.001))
                    stop("Count data must be integer-valued")
                nn = y[,1] + y[,2]
                y = ifelse(nn > 0, y[,1]/nn, 0)
                w = w * nn
                mustart = (0.5 + nn * y) / (1 + nn)
            } else 
                 stop("Response not of the right form")
            predictors.names = namesof("mu", .link, earg= .earg, short=TRUE)
        }
    }), list( .link=link, .mv=mv, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mu =  eta2theta(eta, link= .link, earg = .earg)
        mu
    }, list( .link=link, .earg = earg  ))),
    last=eval(substitute(expression({
        if (is.R()) {
            if (exists("CQO.FastAlgorithm", envir = VGAMenv))
                rm("CQO.FastAlgorithm", envir = VGAMenv)
            if (exists("modelno", envir = VGAMenv))
                rm("modelno", envir = VGAMenv)
        } else {
            while(exists("CQO.FastAlgorithm"))
                remove("CQO.FastAlgorithm")
            while(exists("modelno"))
                remove("modelno")
        }
        dpar <- .dispersion
        if (!dpar) {
            temp87 = (y-mu)^2 * wz /
                     (dtheta.deta(mu, link= .link, earg = .earg )^2) # w cancel
            if (.mv && ! .onedpar) {
                dpar = rep(as.numeric(NA), len=M)
                temp87 = cbind(temp87)
                nrow.mu = if (is.matrix(mu)) nrow(mu) else length(mu)
                for(ii in 1:M)
                    dpar[ii] = sum(temp87[,ii]) / (nrow.mu - ncol(x))
                if (is.matrix(y) && length(dimnames(y)[[2]])==length(dpar))
                    names(dpar) = dimnames(y)[[2]]
            } else 
                dpar = sum(temp87) / (length(mu) - ncol(x))
        }
        misc$mv = .mv
        misc$dispersion <- dpar
        misc$default.dispersion <- 1
        misc$estimated.dispersion <- .estimated.dispersion
        misc$link = rep( .link, length=M)
        names(misc$link) = if (M>1) dn2 else "mu"

        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg

        misc$expected = TRUE
    }), list( .dispersion=dispersion, .estimated.dispersion=estimated.dispersion,
              .onedpar=onedpar, .link=link, .mv=mv, .earg = earg ))),
    link=eval(substitute(function(mu, extra=NULL)
        theta2eta(mu, .link, earg = .earg )
    , list( .link=link, .earg = earg ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        if (residuals) w*(y/mu - (1-y)/(1-mu)) else {
            sum(dbinom(x=w*y, size=w, prob=mu, log=TRUE))
        }
    },
    vfamily=c("binomialff", "vcategorical"),
    deriv=eval(substitute(expression({
        if ( .link == "logit") {
            w * (y - mu)
        } else if ( .link == "cloglog") {
            mu.use = mu
            smallno = 100 * .Machine$double.eps
            mu.use[mu.use < smallno] = smallno
            mu.use[mu.use > 1 - smallno] = 1 - smallno
            -w * (y - mu) * log1p(-mu.use) / mu.use
        } else
            w * dtheta.deta(mu, link= .link, earg = .earg )* (y/mu - 1)/(1-mu)
    }), list( .link=link, .earg = earg ))),
    weight=eval(substitute(expression({
        tmp100 = mu*(1-mu)

        tmp200 = if ( .link == "logit") {
            cbind(w * tmp100)
        } else if ( .link == "cloglog") {
            cbind(w * (1-mu.use) * (log1p(-mu.use))^2 / mu.use )
        } else {
            cbind(w * dtheta.deta(mu, link= .link, earg = .earg)^2 / tmp100)
        }
        for(ii in 1:M) {
            index200 = !is.finite(tmp200[,ii]) |
                       (abs(tmp200[,ii]) < .Machine$double.eps)
            if (any(index200)) {
                tmp200[index200,ii] = .Machine$double.eps # Diagonal 0's are bad 
            }
        }
        tmp200
    }), list( .link=link, .earg = earg ))))
}



 gammaff = function(link="nreciprocal", earg=list(), dispersion=0)
{
    estimated.dispersion <- dispersion==0
    if (mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Gamma distribution\n\n",
           "Link:     ", namesof("mu", link, earg=earg), "\n",
           "Variance: mu^2 / k"),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        devi <- -2 * w * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu)
        if (residuals) {
            sign(y - mu) * sqrt(abs(devi) * w)
        } else sum(w * devi)
    },
    initialize=eval(substitute(expression({
        mustart <- y + 0.167 * (y == 0)
            M = if (is.matrix(y)) ncol(y) else 1
            dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
            dn2 = if (length(dn2)) {
                paste("E[", dn2, "]", sep="") 
            } else {
                paste("mu", 1:M, sep="") 
            }
            predictors.names = namesof(if(M>1) dn2 else "mu", .link,
                 earg=.earg, short=TRUE)
        if (!length(etastart))
            etastart <- theta2eta(mustart, link= .link, earg=.earg)
    }), list( .link=link, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta, link= .link, earg=.earg)
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        dpar <- .dispersion
        if (!dpar) {
            if (M == 1) {
                temp = w * dmu.deta^2
                dpar = sum(w * (y-mu)^2 * wz / temp) / (length(mu) - ncol(x))
            } else {
                dpar = rep(0, len=M)
                for(spp in 1:M) {
                    temp = w * dmu.deta[,spp]^2
                    dpar[spp] = sum(w * (y[,spp]-mu[,spp])^2 * wz[,spp]/temp) /
                                (length(mu[,spp]) - ncol(x))
                }
            }
        }
        misc$dispersion <- dpar
        misc$default.dispersion <- 0
        misc$estimated.dispersion <- .estimated.dispersion
        misc$link = rep( .link, length=M)
        names(misc$link) = if (M>1) paste("mu", 1:M, sep="") else "mu"

        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg

        misc$expected = TRUE
    }), list( .dispersion=dispersion, .earg=earg,
              .estimated.dispersion=estimated.dispersion,
              .link=link ))),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, link= .link, earg=.earg)
    }, list( .link=link, .earg=earg ))),
    vfamily="gammaff",
    deriv=eval(substitute(expression({
        dl.dmu = (y-mu) / mu^2
        dmu.deta = dtheta.deta(theta=mu, link= .link, earg=.earg)
        w * dl.dmu * dmu.deta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        d2l.dmu2 = 1 / mu^2
        w * dmu.deta^2 * d2l.dmu2
    }), list( .link=link, .earg=earg ))))
}



 inverse.gaussianff = function(link="natural.ig", dispersion=0)
{
    estimated.dispersion <- dispersion==0
    warning("@deviance() not finished")
    warning("needs checking, but I'm sure it works")

    if (mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Inverse Gaussian distribution\n\n",
           "Link:     ", namesof("mu", link), "\n",
           "Variance: mu^3 /k"),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        pow <- 3  # Use Quasi()$deviance with pow==3
        devy  <- y^(2-pow) / (1-pow) - y^(2-pow) / (2-pow)
        devmu <- y * mu^(1-pow) / (1-pow) - mu^(2-pow) / (2-pow)
        devi <- 2 * (devy - devmu)
        if (residuals) {
            sign(y - mu) * sqrt(abs(devi) * w)
        } else sum(w * devi)
    },
    initialize=eval(substitute(expression({
        mu <- y + 0.167 * (y == 0)
        if (!length(etastart))
            etastart <- theta2eta(mu, link= .link)
    }), list( .link=link ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta, link= .link)
    }, list( .link=link ))),
    last=eval(substitute(expression({
        dpar <- .dispersion
        if (!dpar) {
            temp <- w * dmu.deta^2
            dpar <- sum( w * (y-mu)^2 * wz / temp ) / (length(mu) - ncol(x))
        }
        misc$dispersion <- dpar
        misc$default.dispersion <- 0
        misc$estimated.dispersion <- .estimated.dispersion
        misc$link = rep( .link, length=M)
        names(misc$link) = if (M>1) paste("mu", 1:M, sep="") else "mu"
    }), list( .dispersion=dispersion,
              .estimated.dispersion=estimated.dispersion,
              .link=link ))),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, link= .link)
    }, list( .link=link ))),
    vfamily="inverse.gaussianff",
    deriv=eval(substitute(expression({
        dl.dmu <- (y-mu) / mu^3
        dmu.deta <- dtheta.deta(theta=mu, link= .link)
        w * dl.dmu * dmu.deta
    }), list( .link=link ))),
    weight=eval(substitute(expression({
        d2l.dmu2 <- 1 / mu^3
        w * dmu.deta^2 * d2l.dmu2
    }), list( .link=link ))))
}




dinv.gaussian = function(x, mu, lambda, log=FALSE) {
    if (!is.logical(log.arg <- log))
        stop("bad input for argument 'log'")
    rm(log)

    LLL = max(length(x), length(mu), length(lambda))
    x = rep(x, len=LLL); mu = rep(mu, len=LLL); lambda = rep(lambda, len=LLL)
    logdensity = rep(log(0), len=LLL)
    xok = (x > 0)
    logdensity[xok] = 0.5 * log(lambda[xok] / (2 * pi * x[xok]^3)) -
                      lambda[xok] * (x[xok]-mu[xok])^2 / (2*mu[xok]^2 * x[xok])
    logdensity[mu     <= 0] = NaN
    logdensity[lambda <= 0] = NaN
    if (log.arg) logdensity else exp(logdensity)
}


pinv.gaussian = function(q, mu, lambda) {
    if (any(mu <=0)) stop("mu must be positive")
    if (any(lambda <=0)) stop("lambda must be positive")
    ans = q
    mu = rep(mu, len=length(q))
    lambda = rep(lambda, len=length(q))
    ans[q <= 0] = 0
    bb = q > 0
    ans[bb] = pnorm(sqrt(lambda[bb]/q[bb])*(q[bb]/mu[bb]-1)) +
              exp(2*lambda[bb]/mu[bb]) *
              pnorm(-sqrt(lambda[bb]/q[bb])*(q[bb]/mu[bb]+1))
    ans
}


rinv.gaussian = function(n, mu, lambda) {
    use.n = if ((length.n <- length(n)) > 1) length.n else
            if (!is.Numeric(n, integ=TRUE, allow=1, posit=TRUE))
                stop("bad input for argument 'n'") else n
    mu = rep(mu, len=use.n); lambda = rep(lambda, len=use.n)
    u = runif(use.n)
    Z = rnorm(use.n)^2 # rchisq(use.n, df=1)
    phi = lambda / mu
    y1 = 1 - 0.5 * (sqrt(Z^2 + 4*phi*Z) - Z) / phi
    ans = mu * ifelse((1+y1)*u > 1, 1/y1, y1)
    ans[mu     <= 0] = NaN
    ans[lambda <= 0] = NaN
    ans
}



 inv.gaussianff = function(lmu="loge", llambda="loge",
                           emu=list(), elambda=list(),
                           ilambda=1,
                           zero=NULL)
{
    if (mode(lmu) != "character" && mode(lmu) != "name")
        lmu <- as.character(substitute(lmu))
    if (mode(llambda) != "character" && mode(llambda) != "name")
        llambda <- as.character(substitute(llambda))
    if (!is.list(emu)) emu = list()
    if (!is.list(elambda)) elambda = list()

    new("vglmff",
    blurb=c("Inverse Gaussian distribution\n\n",
           "f(y) = sqrt(lambda/(2*pi*y^3)) * exp(-lambda*(y-mu)^2/(2*mu^2*y)), y&lambda>0",
           "Link:     ", namesof("mu", lmu, earg= emu), ", ",
                         namesof("lambda", llambda, earg= elambda), "\n",
           "Mean:     ", "mu\n",
           "Variance: mu^3 / lambda"),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if (any(y <= 0)) stop("Require the response to have positive values")
        predictors.names =
        c(namesof("mu", .lmu, earg= .emu, short= TRUE),
          namesof("lambda", .llambda, earg= .elambda, short= TRUE))
        if (!length(etastart)) {
            initmu = y + 1/8
            initlambda = rep(if(length( .ilambda)) .ilambda else 1, len=n)
            etastart = cbind(
                theta2eta(initmu, link=.lmu, earg= .emu), 
                theta2eta(initlambda, link=.llambda, earg= .elambda))
        }
    }), list( .lmu=lmu, .llambda=llambda,
              .emu=emu, .elambda=elambda,
              .ilambda=ilambda ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], link=.lmu, earg= .emu)
    }, list( .lmu=lmu, .emu=emu, .elambda=elambda ))),
    last=eval(substitute(expression({
        misc$link = c(mu = .lmu, lambda = .llambda)
        misc$earg = list(mu = .emu, lambda = .elambda)
    }), list( .lmu=lmu, .llambda=llambda, .emu=emu, .elambda=elambda ))),
    loglikelihood=eval(substitute(
             function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        lambda <- eta2theta(eta[,2], link=.llambda, earg= .elambda)
        if (residuals) stop("loglikelihood residuals not implemented yet") else {
            sum(w * dinv.gaussian(x=y, mu=mu, lambda=lambda, log=TRUE))
        }
    }, list( .llambda=llambda, .emu=emu, .elambda=elambda ))),
    vfamily="inv.gaussianff",
    deriv=eval(substitute(expression({
        lambda <- eta2theta(eta[,2], link=.llambda, earg= .elambda)
        dl.dmu = lambda * (y-mu) / mu^3
        dl.dlambda <- 0.5 / lambda - (y-mu)^2 / (2 * mu^2 * y)
        dmu.deta <- dtheta.deta(theta=mu, link=.lmu, earg= .emu)
        dlambda.deta <- dtheta.deta(theta=lambda, link=.llambda, earg= .elambda)
        w * cbind(dl.dmu * dmu.deta, dl.dlambda * dlambda.deta)
    }), list( .lmu=lmu, .llambda=llambda, .emu=emu, .elambda=elambda ))),
    weight=eval(substitute(expression({
        d2l.dmu2 = lambda / mu^3
        d2l.dlambda2 = 0.5 / (lambda^2)
        w * cbind(dmu.deta^2 * d2l.dmu2, dlambda.deta^2 * d2l.dlambda2)
    }), list( .lmu=lmu, .llambda=llambda, .emu=emu, .elambda=elambda ))))
}



 poissonff = function(link="loge", earg=list(),
                      dispersion=1, onedpar=FALSE,
                      imu=NULL, method.init=1,
                      parallel=FALSE, zero=NULL)
{

    estimated.dispersion <- dispersion==0
    if (mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if (!is.list(earg)) earg = list()
    if (!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 3)
        stop("'method.init' must be 1 or 2 or 3")
    if (length(imu) && !is.Numeric(imu, posit=TRUE))
        stop("bad input for argument 'imu'")

    new("vglmff",
    blurb=c("Poisson distribution\n\n",
           "Link:     ", namesof("mu", link, earg= earg), "\n",
           "Variance: mu"),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(1,M,1), x, .parallel, constraints)
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero=zero ))),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        nz = y > 0
        devi =  -(y - mu)
        devi[nz] = devi[nz] + y[nz] * log(y[nz]/mu[nz])
        if (residuals) sign(y - mu) * sqrt(2 * abs(devi) * w) else
            2 * sum(w * devi)
    },
    initialize=eval(substitute(expression({
        y = as.matrix(y)
        M = ncoly = ncol(y)

        if (is.R()) assign("CQO.FastAlgorithm",
            ( .link == "loge"), envir = VGAMenv) else
            CQO.FastAlgorithm <<- ( .link == "loge")
        dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
        dn2 = if (length(dn2)) {
            paste("E[", dn2, "]", sep="") 
        } else {
            paste("mu", 1:M, sep="") 
        }
        predictors.names = namesof(if(M>1) dn2 else "mu", .link,
            earg= .earg, short=TRUE)

        if (!length(etastart)) {
            mu.init = pmax(y, 1/8)
            for(iii in 1:ncol(y)) {
                if ( .method.init == 2) {
                    mu.init[,iii] = weighted.mean(y[,iii], w) + 1/8
                } else if ( .method.init == 3) {
                    mu.init[,iii] = median(y[,iii]) + 1/8
                }
            }
            if (length(.imu))
                mu.init = matrix( .imu, n, ncoly, byrow=TRUE)
            etastart <- theta2eta(mu.init, link= .link, earg= .earg)
        }
    }), list( .link=link, .estimated.dispersion=estimated.dispersion,
              .method.init=method.init, .imu=imu, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mu = eta2theta(eta, link= .link, earg= .earg)
        mu
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
       if (is.R()) {
            if (exists("CQO.FastAlgorithm", envir = VGAMenv))
                rm("CQO.FastAlgorithm", envir = VGAMenv)
        } else {
            while(exists("CQO.FastAlgorithm"))
                remove("CQO.FastAlgorithm")
        }
        dpar <- .dispersion
        if (!dpar) {
            temp87 = (y-mu)^2 *
                wz / (dtheta.deta(mu, link= .link, earg= .earg)^2) # w cancel
            if (M > 1 && ! .onedpar) {
                dpar = rep(as.numeric(NA), len=M)
                temp87 = cbind(temp87)
                nrow.mu = if (is.matrix(mu)) nrow(mu) else length(mu)
                for(ii in 1:M)
                    dpar[ii] = sum(temp87[,ii]) / (nrow.mu - ncol(x))
                if (is.matrix(y) && length(dimnames(y)[[2]])==length(dpar))
                    names(dpar) = dimnames(y)[[2]]
            } else 
                dpar = sum(temp87) / (length(mu) - ncol(x))
        }
        misc$dispersion <- dpar
        misc$default.dispersion <- 1
        misc$estimated.dispersion <- .estimated.dispersion
        misc$expected = TRUE
        misc$link = rep( .link, length=M)
        names(misc$link) = if (M>1) dn2 else "mu"
        misc$method.init = .method.init

        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg
    }), list( .dispersion=dispersion, .method.init=method.init,
              .estimated.dispersion=estimated.dispersion,
              .onedpar=onedpar, .link=link, .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, link= .link, earg= .earg)
    }, list( .link=link, .earg=earg ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        if (residuals) w*(y/mu - 1) else {
            sum(w * dpois(x=y, lambda=mu, log=TRUE))
        }
    },
    vfamily="poissonff",
    deriv=eval(substitute(expression({
        if ( .link == "loge" && (any(mu < .Machine$double.eps))) {
            w * (y - mu)
        } else {
            lambda <- mu
            dl.dlambda <- (y-lambda) / lambda
            dlambda.deta <- dtheta.deta(theta=lambda, link= .link, earg= .earg)
            w * dl.dlambda * dlambda.deta
        }
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        if ( .link == "loge" && (any(mu < .Machine$double.eps))) {
            tmp600 = mu
            tmp600[tmp600 < .Machine$double.eps] = .Machine$double.eps
            w * tmp600
        } else {
            d2l.dlambda2 = 1 / lambda
            d2lambda.deta2=d2theta.deta2(theta=lambda,link= .link,earg= .earg)
            w * dlambda.deta^2 * d2l.dlambda2
        }
    }), list( .link=link, .earg=earg ))))
}


 quasibinomialff = function(link = "logit", mv = FALSE, onedpar = !mv, 
                            parallel = FALSE, zero = NULL) {
    dispersion = 0 # Estimated; this is the only difference with binomialff()
    ans =
    binomialff(link = link, dispersion=dispersion, mv=mv, onedpar=onedpar,
               parallel=parallel, zero=zero) 
    ans@vfamily = "quasibinomialff"
    ans
}

 quasipoissonff = function(link = "loge", onedpar = FALSE, parallel = FALSE,
                           zero = NULL) {
    dispersion = 0 # Estimated; this is the only difference with poissonff()
    ans =
    poissonff(link = link, dispersion=dispersion, onedpar=onedpar,
               parallel=parallel, zero=zero) 
    ans@vfamily = "quasipoissonff"
    ans
}




poissonqn.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}


 poissonqn = function(link="loge", earg=list(),
                      dispersion=1, onedpar=FALSE,
                      parallel=FALSE, zero=NULL,
                      wwts=c("expected","observed","qn"))
{
    estimated.dispersion <- dispersion==0
    if (mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if (mode(wwts) != "character" && mode(wwts) != "name")
        wwts <- as.character(substitute(wwts))
    wwts <- match.arg(wwts, c("expected","observed","qn"))[1]
    if (!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Poisson distribution\n\n",
           "Link:     ", namesof("mu", link, earg= earg), "\n",
           "Variance: mu"),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(1,M,1), x, .parallel, constraints)
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero=zero ))),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        nz = y > 0
        devi =  -(y - mu)
        devi[nz] = devi[nz] + y[nz] * log(y[nz]/mu[nz])
        if (residuals) sign(y - mu) * sqrt(2 * abs(devi) * w) else
            2 * sum(w * devi)
    },
    initialize=eval(substitute(expression({
        M = if (is.matrix(y)) ncol(y) else 1
        dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
        dn2 = if (length(dn2)) {
            paste("E[", dn2, "]", sep="") 
        } else {
            paste("mu", 1:M, sep="") 
        }
        predictors.names = namesof(if(M>1) dn2 else "mu", .link,
            earg= .earg, short=TRUE)
        mu = pmax(y, 0.167)  # y + 0.167 * (y == 0)
        if (!length(etastart))
            etastart <- theta2eta(mu, link= .link, earg= .earg)
    }), list( .link=link, .estimated.dispersion=estimated.dispersion,
              .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta, link= .link, earg= .earg)
    }, list( .link=link,
              .earg=earg ))),
    last=eval(substitute(expression({
        dpar <- .dispersion
        if (!dpar) {
            temp87= (y-mu)^2 * wz/(dtheta.deta(mu, link= .link, earg= .earg)^2)
            if (M > 1 && ! .onedpar) {
                dpar = rep(as.numeric(NA), len=M)
                temp87 = cbind(temp87)
                nrow.mu = if (is.matrix(mu)) nrow(mu) else length(mu)
                for(i in 1:M)
                    dpar[i] = sum(temp87[,i]) / (nrow.mu - ncol(x))
                if (is.matrix(y) && length(dimnames(y)[[2]])==length(dpar))
                    names(dpar) = dimnames(y)[[2]]
            } else 
                dpar = sum(temp87) / (length(mu) - ncol(x))
        }
        misc$BFGS = TRUE
        misc$dispersion <- dpar
        misc$default.dispersion <- 1
        misc$estimated.dispersion <- .estimated.dispersion
        misc$expected = FALSE
        misc$link = rep( .link, length=M)
        names(misc$link) = if (M>1) dn2 else "mu"

        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg
    }), list( .dispersion=dispersion,
              .earg=earg, 
              .estimated.dispersion=estimated.dispersion,
              .onedpar=onedpar, .link=link ))),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, link= .link, earg= .earg)
    }, list( .link=link,
              .earg=earg ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        if (residuals) w*(y/mu - 1) else {
            sum(w * dpois(x=y, lambda=mu, log=TRUE))
        }
    },
    vfamily="poissonqn",
    deriv=eval(substitute(expression({
        if (iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }

        derivnew =
        if ( .link == "loge" && (any(mu < .Machine$double.eps))) {
            w * (y - mu)
        } else {
            lambda <- mu
            dl.dlambda <- (y-lambda) / lambda
            dlambda.deta <- dtheta.deta(theta=lambda, link= .link, earg= .earg)
            w * dl.dlambda * dlambda.deta
        }
        derivnew
    }), list( .link=link,
              .earg=earg ))),
    weight=eval(substitute(expression({
        if ( .wwts == "qn") {
            if (iter == 1) {
                wznew = cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
            } else {
                wzold = wznew
                wznew = qnupdate(w=w, wzold=wzold, dderiv=(derivold-derivnew),
                                 deta=etanew-etaold, M=M,
                                 trace=trace)  # weights incorporated in args
            }
        } else if ( .wwts == "expected") {
            wznew = if ( .link == "loge") {
                tmp600 = mu
                tmp600[tmp600 < .Machine$double.eps] = .Machine$double.eps
                w * tmp600
            } else {
                d2l.dlambda2 = 1 / lambda
                w * dlambda.deta^2 * d2l.dlambda2
            }
        } else {
            wznew = if ( .link == "loge") {
                tmp600 = y
                tmp600[y < .Machine$double.eps] = sqrt(.Machine$double.eps)
                w * tmp600
            } else {
                stop("this is not programmed in yet")
            }
        }
        wznew
    }), list( .wwts=wwts, .link=link,
              .earg=earg ))))
}




 dexppoisson = function(lmean="loge", emean=list(),
                        ldispersion="logit", edispersion=list(),
                        idispersion=0.8,
                        zero=NULL)
{
    if (mode(lmean)!= "character" && mode(lmean)!= "name")
        lmean = as.character(substitute(lmean))
    if (mode(ldispersion)!= "character" && mode(ldispersion)!= "name")
        ldispersion = as.character(substitute(ldispersion))
    if (!is.Numeric(idispersion, posit=TRUE))
        stop("bad input for 'idispersion'")
    if (!is.list(emean)) emean = list()
    if (!is.list(edispersion)) edispersion = list()

    new("vglmff",
    blurb=c("Double Exponential Poisson distribution\n\n",
           "Link:     ",
           namesof("mean", lmean, earg= emean), ", ",
           namesof("dispersion", lmean, earg= edispersion), "\n",
           "Mean:     ", "mean\n",
           "Variance: mean / dispersion"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if (ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        M = if (is.matrix(y)) ncol(y) else 1
        dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
        dn2 = if (length(dn2)) {
            paste("E[", dn2, "]", sep="") 
        } else {
            "mu"
        }
        predictors.names =
            c(namesof(dn2, link= .lmean, earg= .emean, short=TRUE),
              namesof("dispersion", link= .ldispersion,
                                    earg= .edispersion, short=TRUE))
        init.mu = pmax(y, 1/8)
        if (!length(etastart))
            etastart = cbind(theta2eta(init.mu, link= .lmean,earg= .emean),
                             theta2eta(rep( .idispersion, len=n),
                                       link= .ldispersion, earg= .edispersion))
    }), list( .lmean=lmean, .emean=emean,
              .ldispersion=ldispersion, .edispersion=edispersion,
              .idispersion=idispersion ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], link= .lmean, earg= .emean)
    }, list( .lmean=lmean, .emean=emean,
             .ldispersion=ldispersion, .edispersion=edispersion ))),
    last=eval(substitute(expression({
        misc$expected = TRUE
        misc$link = c("mean"= .lmean, "dispersion"= .ldispersion)
        misc$earg = list(mean= .emean, dispersion= .edispersion)
    }), list( .lmean=lmean, .emean=emean,
              .ldispersion=ldispersion, .edispersion=edispersion ))),
    loglikelihood=eval(substitute(
                      function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        lambda = eta2theta(eta[,1], link= .lmean, earg= .emean)
        Disper = eta2theta(eta[,2], link= .ldispersion, earg= .edispersion)
        if (residuals) stop("loglikelihood residuals not implemented yet") else {
            sum(w*(0.5*log(Disper) + Disper*(y-lambda) + Disper*y*log(lambda)))
        }
    }, list( .lmean=lmean, .emean=emean,
             .ldispersion=ldispersion, .edispersion=edispersion ))),
    vfamily="dexppoisson",
    deriv=eval(substitute(expression({
        lambda = eta2theta(eta[,1], link= .lmean, earg= .emean)
        Disper = eta2theta(eta[,2], link= .ldispersion, earg= .edispersion)
        dl.dlambda = Disper * (y / lambda - 1)
        dl.dDisper = y * log(lambda) + y - lambda + 0.5 / Disper
        dlambda.deta = dtheta.deta(theta=lambda, link= .lmean, earg= .emean)
        dDisper.deta = dtheta.deta(theta=Disper, link= .ldispersion,
                                   earg= .edispersion)
        w * cbind(dl.dlambda * dlambda.deta,
                  dl.dDisper * dDisper.deta)
    }), list( .lmean=lmean, .emean=emean,
              .ldispersion=ldispersion, .edispersion=edispersion ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), nrow=n, ncol=2) # diagonal
        usethis.lambda = pmax(lambda, .Machine$double.eps / 10000)
        wz[,iam(1,1,M)] = (Disper / usethis.lambda) * dlambda.deta^2
        wz[,iam(2,2,M)] = (0.5 / Disper^2) * dDisper.deta^2
        w * wz
    }), list( .lmean=lmean, .emean=emean,
              .ldispersion=ldispersion, .edispersion=edispersion ))))
}



 dexpbinomial = function(lmean="logit", ldispersion="logit",
                         emean=list(), edispersion=list(),
                         idispersion=0.25,
                         zero=2)
{
    if (mode(lmean)!= "character" && mode(lmean)!= "name")
        lmean = as.character(substitute(lmean))
    if (mode(ldispersion)!= "character" && mode(ldispersion)!= "name")
        ldispersion = as.character(substitute(ldispersion))
    if (!is.Numeric(idispersion, posit=TRUE))
        stop("bad input for 'idispersion'")
    if (!is.list(emean)) emean = list()
    if (!is.list(edispersion)) edispersion = list()

    new("vglmff",
    blurb=c("Double Exponential Binomial distribution\n\n",
           "Link:     ",
           namesof("mean", lmean, earg= emean), ", ",
           namesof("dispersion", lmean, earg= edispersion), "\n",
           "Mean:     ", "mean\n"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if (ncol(cbind(w)) != 1)
            stop("'weights' must be a vector or a one-column matrix")

            NCOL = function (x)
                if (is.array(x) && length(dim(x)) > 1 ||
                is.data.frame(x)) ncol(x) else as.integer(1)

            if (NCOL(y) == 1) {
                if (is.factor(y)) y = y != levels(y)[1]
                nn = rep(1, n)
                if (!all(y >= 0 & y <= 1))
                    stop("response values must be in [0, 1]")
                init.mu = (0.5 + w * y) / (1 + w)
                no.successes = w * y
                if (any(abs(no.successes - round(no.successes)) > 0.001))
                    stop("Number of successes must be integer-valued")
            } else if (NCOL(y) == 2) {
                if (any(abs(y - round(y)) > 0.001))
                    stop("Count data must be integer-valued")
                nn = y[,1] + y[,2]
                y = ifelse(nn > 0, y[,1]/nn, 0)
                w = w * nn
                init.mu = (0.5 + nn * y) / (1 + nn)
            } else
                 stop("Response not of the right form")


        dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
        dn2 = if (length(dn2)) {
            paste("E[", dn2, "]", sep="") 
        } else {
            "mu"
        }
        predictors.names =
            c(namesof(dn2, link= .lmean, earg= .emean, short=TRUE),
              namesof("dispersion", link= .ldispersion,
                                    earg= .edispersion, short=TRUE))
        if (!length(etastart))
            etastart = cbind(theta2eta(init.mu, link= .lmean,earg= .emean),
                             theta2eta(rep( .idispersion, len=n),
                                       link= .ldispersion, earg= .edispersion))
    }), list( .lmean=lmean, .emean=emean,
              .ldispersion=ldispersion, .edispersion=edispersion,
              .idispersion=idispersion ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], link= .lmean, earg= .emean)
    }, list( .lmean=lmean, .emean=emean,
             .ldispersion=ldispersion, .edispersion=edispersion ))),
    last=eval(substitute(expression({
        misc$expected = TRUE
        misc$link = c("mean"= .lmean, "dispersion"= .ldispersion)
        misc$earg = list(mean= .emean, dispersion= .edispersion)
    }), list( .lmean=lmean, .emean=emean,
              .ldispersion=ldispersion, .edispersion=edispersion ))),
    loglikelihood=eval(substitute(
                      function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        prob = eta2theta(eta[,1], link= .lmean, earg= .emean)
        Disper = eta2theta(eta[,2], link= .ldispersion, earg= .edispersion)
        if (residuals) stop("loglikelihood residuals not implemented yet") else {
            temp1 = y * log(ifelse(y > 0, y, 1)) # y*log(y)
            temp2 = (1.0-y) * log1p(ifelse(y < 1, -y, 0)) # (1-y)*log(1-y)
            sum(0.5*log(Disper) + w*(y*Disper*log(prob) +
                   (1-y)*Disper*log1p(-prob) +
                   temp1*(1-Disper) + temp2*(1-Disper)))
        }
    }, list( .lmean=lmean, .emean=emean,
             .ldispersion=ldispersion, .edispersion=edispersion ))),
    vfamily="dexpbinomial",
    deriv=eval(substitute(expression({
        prob = eta2theta(eta[,1], link= .lmean, earg= .emean)
        Disper = eta2theta(eta[,2], link= .ldispersion, earg= .edispersion)
        temp1 = y * log(ifelse(y > 0, y, 1)) # y*log(y)
        temp2 = (1.0-y) * log1p(ifelse(y < 1, -y, 0)) # (1-y)*log(1-y)
        temp3 = prob * (1.0-prob)
        temp3 = pmax(temp3, .Machine$double.eps * 10000)
        dl.dprob = w * Disper * (y - prob) / temp3
        dl.dDisper = 0.5 / Disper + w * (y * log(prob) + 
                     (1-y)*log1p(-prob) - temp1 - temp2)
        dprob.deta = dtheta.deta(theta=prob, link= .lmean, earg= .emean)
        dDisper.deta = dtheta.deta(theta=Disper, link= .ldispersion,
                                   earg= .edispersion)
        cbind(dl.dprob * dprob.deta,
              dl.dDisper * dDisper.deta)
    }), list( .lmean=lmean, .emean=emean,
              .ldispersion=ldispersion, .edispersion=edispersion ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), nrow=n, ncol=2) # diagonal
        wz[,iam(1,1,M)] = w * (Disper / temp3) * dprob.deta^2
        wz[,iam(2,2,M)] = (0.5 / Disper^2) * dDisper.deta^2
        wz
    }), list( .lmean=lmean, .emean=emean,
              .ldispersion=ldispersion, .edispersion=edispersion ))))
}




 mbinomial = function(mvar=NULL, link="logit", earg=list(),
                      parallel = TRUE, smallno = .Machine$double.eps^(3/4))
{
    if (mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if (!is.list(earg)) earg = list()
    if (!is.Numeric(smallno, positive=TRUE, allow=1) || smallno > 1e-4)
        stop("bad input for 'smallno'")
    if (is.logical(parallel) && !parallel)
        stop("'parallel' must be TRUE")

    temp = terms(mvar)
    mvar = attr(temp,"term.labels")
    if (length(mvar) != 1) stop("cannot obtain the matching variable")
    if (!is.character(mvar) || length(mvar) != 1) {
        stop("bad input for 'mvar'")
    }

    new("vglmff",
    blurb= c("Matched binomial model (intercepts fitted)\n\n", 
           "Link:     ", namesof("mu[,j]", link, earg= earg)),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(1,M,1), x, .parallel, constraints,
                               intercept.apply=TRUE)
        constraints[[extra$mvar]] <- diag(M)

        specialCM = list(a = vector("list", M-1))
        for(ii in 1:(M-1)) {
            specialCM[[1]][[ii]] = (constraints[[extra$mvar]])[,1+ii,drop=FALSE]
        }
        names(specialCM) = extra$mvar
    }), list( .parallel=parallel ))),
    initialize=eval(substitute(expression({
        mvar = .mvar

        NCOL = function (x) 
            if (is.array(x) && length(dim(x)) > 1 ||
            is.data.frame(x)) ncol(x) else as.integer(1)

        if (NCOL(y) == 1) {
            if (is.factor(y)) y = y != levels(y)[1]
            nn = rep(1, n)
            if (!all(y >= 0 & y <= 1))
                stop("response values must be in [0, 1]")
            mustart = (0.5 + w * y) / (1 + w)
            no.successes = w * y
            if (any(abs(no.successes - round(no.successes)) > 0.001))
                stop("Number of successes must be integer-valued")
        } else if (NCOL(y) == 2) {
            if (any(abs(y - round(y)) > 0.001))
                stop("Count data must be integer-valued")
            nn = y[,1] + y[,2]
            y = ifelse(nn > 0, y[,1]/nn, 0)
            w = w * nn
            mustart = (0.5 + nn * y) / (1 + nn)
        } else 
             stop("Response not of the right form")

        temp1 = attr(x, "assign")
        if (colnames(x)[1] != "(Intercept)") stop("x must have an intercept")
        M = CCC = length(temp1[[mvar]]) + (colnames(x)[1] == "(Intercept)")
        temp9 = x[,temp1[[mvar]],drop=FALSE]
        temp9 = temp9 * matrix(2:CCC, n, CCC-1, byrow=TRUE)
        temp9 = apply(temp9, 1, max)
        temp9[temp9 == 0] = 1
        extra$NoMatchedSets = CCC
        extra$n = n
        extra$M = M
        extra$mvar = mvar
        extra$index9 = temp9

        predictors.names = namesof("mu", .link, earg= .earg, short=TRUE)
        predictors.names = rep(predictors.names, len=M)
    }), list( .link=link, .earg=earg, .mvar=mvar ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mu = eta2theta(eta, link= .link, earg = .earg)
        mu[cbind(1:extra$n, extra$index9)]
    }, list( .link=link, .earg = earg  ))),
    last=eval(substitute(expression({
        misc$link = rep( .link, length=M)
        names(misc$link) = if (M>1) paste("mu(matched set ",
            1:M, ")", sep="") else "mu"
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg

        misc$expected = TRUE
    }), list( .link=link, .earg = earg ))),
    link=eval(substitute(function(mu, extra=NULL) {
        temp = theta2eta(mu, .link, earg = .earg )
        matrix(temp, extra$n, extra$M)
    }, list( .link=link, .earg = earg ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        if (residuals) w*(y/mu - (1-y)/(1-mu)) else {
            sum(w*(y*log(mu) + (1-y)*log1p(-mu)))
        }
    },
    vfamily=c("mbinomial", "vcategorical"),
    deriv=eval(substitute(expression({
        answer =
        if ( .link == "logit") {
            w * (y - mu)
        } else if ( .link == "cloglog") {
            mu.use = mu
            smallno = 100 * .Machine$double.eps
            mu.use[mu.use < smallno] = smallno
            mu.use[mu.use > 1 - smallno] = 1 - smallno
            -w * (y - mu) * log1p(-mu.use) / mu.use
        } else
            w * dtheta.deta(mu, link= .link, earg = .earg )* (y/mu - 1)/(1-mu)
        result = matrix(0, n, M)
        result[cbind(1:n, extra$index9)] = answer
        result
    }), list( .link=link, .earg = earg ))),
    weight=eval(substitute(expression({
        tmp100 = mu*(1-mu)
        answer = if ( .link == "logit") {
            cbind(w * tmp100)
        } else if ( .link == "cloglog") {
            cbind(w * (1-mu.use) * (log1p(-mu.use))^2 / mu.use )
        } else {
            cbind(w * dtheta.deta(mu, link= .link, earg = .earg)^2 / tmp100)
        }

        result = matrix( .smallno, n, M)
        result[cbind(1:n, extra$index9)] = answer
        result
    }), list( .link=link, .earg = earg, .smallno=smallno ))))
}




mypool = function(x, index) {
    answer = x
    uindex = unique(index)
    for(ii in uindex) {
        ind0 = index == ii
        answer[ind0] = sum(x[ind0])
    }
    answer
}


 if (FALSE)
 mbino = function()
{
    link = "logit"
    earg = list()
    parallel = TRUE

    if (mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if (!is.list(earg)) earg = list()
    if (is.logical(parallel) && !parallel)
        stop("'parallel' must be TRUE")


    new("vglmff",
    blurb= c("Matched binomial model (intercepts not fitted)\n\n", 
           "Link:     ", namesof("mu[,j]", link, earg= earg)),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(1,M,1), x, .parallel, constraints,
                               intercept.apply=FALSE)
    }), list( .parallel=parallel ))),
    initialize=eval(substitute(expression({
        if (colnames(x)[1] == "(Intercept)")
            stop("the model matrix must not have an intercept")

        NCOL = function (x) 
            if (is.array(x) && length(dim(x)) > 1 ||
            is.data.frame(x)) ncol(x) else as.integer(1)

        if (NCOL(y) == 1) {
            if (is.factor(y)) y = y != levels(y)[1]
            nn = rep(1, n)
            if (!all(y >= 0 & y <= 1))
                stop("response values must be in [0, 1]")
            mustart = (0.5 + w * y) / (1 + w)
            no.successes = w * y
            if (any(abs(no.successes - round(no.successes)) > 0.001))
                stop("Number of successes must be integer-valued")
        } else if (NCOL(y) == 2) {
            if (any(abs(y - round(y)) > 0.001))
                stop("Count data must be integer-valued")
            nn = y[,1] + y[,2]
            y = ifelse(nn > 0, y[,1]/nn, 0)
            w = w * nn
            mustart = (0.5 + nn * y) / (1 + nn)
        } else 
             stop("Response not of the right form")

        if (!length(etastart))
            etastart <- theta2eta(mustart, link= "logit", earg= list())

        temp1 = attr(x, "assign")
        mvar = extra$mvar
        if (length(mvar) != n) stop("input extra$mvar doesn't look right")

        if (any(y != 0 & y != 1))
            stop("response vector must have 0 or 1 values only")
        xrle = rle(mvar)
        if (length(unique(mvar)) != length(xrel$length))
            stop("extra$mvar must take on contiguous values")

        temp9 = factor(mvar)
        extra$NoMatchedSets = levels(temp9)
        extra$n = n
        extra$M = M
        extra$rlex = xrle
        extra$index9 = temp9
        predictors.names = namesof("mu", .link, earg= .earg, short=TRUE)
    }), list( .link=link, .earg=earg, .mvar=mvar ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        denominator = exp(eta)
        numerator = mypool(denominator, extra$mvar)
        numerator / denominator
    }, list( .link=link, .earg = earg  ))),
    last=eval(substitute(expression({
        misc$link = c(mu = .link)
        misc$earg = list( mu = .earg )
        misc$expected = TRUE
    }), list( .link=link, .earg = earg ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        if (residuals) w*(y/mu - (1-y)/(1-mu)) else {
            sum(w*(y*log(mu) + (1-y)*log1p(-mu)))
        }
    },
    vfamily=c("mbin", "vcategorical"),
    deriv=eval(substitute(expression({
        answer =
        if ( .link == "logit") {
            w * (y - mu)
        } else stop("can only handle the logit link")
        answer
    }), list( .link=link, .earg = earg ))),
    weight=eval(substitute(expression({
        tmp100 = mu*(1-mu)
        answer = if ( .link == "logit") {
            cbind(w * tmp100)
        } else stop("can only handle the logit link")

        result = matrix( .smallno, n, M)
        result[cbind(1:n, extra$index9)] = answer
        result
    }), list( .link=link, .earg = earg, .smallno=smallno ))))
}







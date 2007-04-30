# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.





if(FALSE)
quasiff = function(link="polw",
                   earg=if(link=="powl") list(power=1) else list(),
                   dispersion=0)
{
    warning("link=powl doesn't work yet")
    estimated.dispersion <- dispersion==0

    if(mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()



    result <- 
    new("vglmff",
    blurb=c("Quasi family\n\n",
           "Link:     ", namesof("mu", link, earg=earg), "\n",
           "Variance: ", ifelse(power.variance==1, variance,
           paste(variance, "^", power.variance, sep=""))),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        pow <- extra$power.variance
        thing <- extra$variance
        if(thing=="mu" && (pow==1 || pow==2))
            stop("this deviance function not right")

        devy  <- y^(2-pow) / (1-pow) - y^(2-pow) / (2-pow)
        devmu <- y * mu^(1-pow) / (1-pow) - mu^(2-pow) / (2-pow)
        devi <- 2 * (devy - devmu)
        if(residuals) {
            sign(y - mu) * sqrt(abs(devi) * w)
        } else sum(w * devi)
    },
    initialize=eval(substitute(expression({
        extra$link <- .link
        extra$variance <- .variance
        extra$power.variance <- .power.variance

        if(.variance=="mu(1-mu)") {
            delete.zero.colns <- TRUE 
            eval(process.categorical.data.vgam)
    
            mustart <- mustart[,1]
            y <- y[,1]
        } else {
            mustart <- y + 0.167 * (y == 0)
        }


    }), list( .link=link, .variance=variance,
              .earg=earg, .power.variance=power.variance ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta, link= .link, earg=.earg)
    }, list( .link=link,
              .earg=earg ))),
    last=eval(substitute(expression({
        dpar <- .dispersion
        if(!dpar) {
            temp <- w * dmu.deta^2
            dpar <- sum( w * (y-mu)^2 * wz / temp ) / (length(mu) - ncol(x))
        }
        misc$dispersion <- dpar
        misc$default.dispersion <- 0
        misc$estimated.dispersion <- .estimated.dispersion
        misc$power.variance <- .power.variance
        misc$link = c("mu" = .link )
    }), list( .dispersion=dispersion,
              .earg=earg, .estimated.dispersion=estimated.dispersion,
            .link=link, .power.variance=power.variance ))),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, link= .link, earg=.earg)
    }, list( .link=link,
              .earg=earg ))),
    vfamily="quasiff",
    deriv=eval(substitute(expression({
        pow <- extra$power.variance
        thing <- extra$variance
        dQ.dmu <- if(thing=="mu") (y-mu)/mu^pow else (y-mu)/(mu*(1-mu))^pow
        dmu.deta <- dtheta.deta(theta=mu, link= .link, earg=.earg)
        w * dQ.dmu * dmu.deta
    }), list( .link=link, .power.variance=power.variance,
              .earg=earg ))),
    weight=eval(substitute(expression({
        d2Q.dmu2 <- if(thing=="mu") 1 / mu^pow else 
            1 / (mu*(1-mu))^pow
        w * dmu.deta^2 * d2Q.dmu2
    }), list( .link=link, .power.variance=power.variance,
              .earg=earg ))))

    if(variance=="mu") {
        if(power.variance==1)
            result@deviance <- poissonff()@deviance
        if(power.variance==2)
            result@deviance <- gammaff()@deviance
    } else {
        result@deviance <- if(power.variance==1) binomialff()@deviance else NULL 
    }

    result
}



binomialff <- function(link="logit", earg=list(),
                       dispersion=1, mv=FALSE, onedpar=!mv,
                       parallel = FALSE,
                       zero=NULL)

{


    estimated.dispersion <- dispersion==0
    if(mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

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
        devy <- y
        nz <- y != 0
        devy[nz] <- y[nz] * log(y[nz])
        nz <- (1 - y) != 0
        devy[nz] <- devy[nz] + (1 - y[nz]) * log(1 - y[nz])
        devmu <- y * log(mu) + (1 - y) * log(1 - mu)
        if(any(small <- mu * (1 - mu) < .Machine$double.eps)) {
            warning("fitted values close to 0 or 1")
            smu <- mu[small]
            sy <- y[small]
            smu <- ifelse(smu < .Machine$double.eps, 
                        .Machine$double.eps, smu)
            onemsmu <- ifelse((1 - smu) < .Machine$
                        double.eps, .Machine$double.eps, 1 - smu)
            devmu[small] <- sy * log(smu) + (1 - sy) * log(onemsmu)
        }
        devi <- 2 * (devy - devmu)
        if(residuals) {
            sign(y - mu) * sqrt(abs(devi) * w)
        } else sum(w * devi)
    },
    initialize=eval(substitute(expression({
        if(is.R()) {
            assign("CQO.FastAlgorithm", ( .link=="logit" || .link=="cloglog"),
                   envir=VGAMenv)
            assign("modelno", if( .link=="logit") 1 else
                   if( .link=="cloglog") 4 else NULL, envir=VGAMenv)
        }  else {
            CQO.FastAlgorithm <<- ( .link == "logit" || .link=="cloglog")
          modelno <<- if( .link=="logit") 1 else if( .link=="cloglog") 4 else NULL
        }
        if(.mv) {
            y = as.matrix(y)
            M = ncol(y)
            if(!all(y == 0 | y == 1))
                stop("response must contain 0's and 1's only")
            dn2 = if(is.matrix(y)) dimnames(y)[[2]] else NULL
            dn2 = if(length(dn2)) {
                paste("E[", dn2, "]", sep="") 
            } else {
                paste("mu", 1:M, sep="") 
            }
            predictors.names = namesof(if(M>1) dn2 else
                "mu", .link, earg= .earg, short=TRUE)

            mustart = (0.5 + w * y) / (1 + w)
        } else {
            NCOL = function (x) 
                if(is.array(x) && length(dim(x)) > 1 ||
                is.data.frame(x)) ncol(x) else as.integer(1)

            if(NCOL(y) == 1) {
                if(is.factor(y)) y = y != levels(y)[1]
                nn = rep(1, n)
                if(!all(y >= 0 & y <= 1))
                    stop("response values must be in [0, 1]")
                mustart = (0.5 + w * y) / (1 + w)
                no.successes = w * y
                if(any(abs(no.successes - round(no.successes)) > 0.001))
                    stop("Number of successes must be integer-valued")
            } else if(NCOL(y) == 2) {
                if(any(abs(y - round(y)) > 0.001))
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
        if(is.R()) {
            if(exists("CQO.FastAlgorithm", envir = VGAMenv))
                rm("CQO.FastAlgorithm", envir = VGAMenv)
            if(exists("modelno", envir = VGAMenv))
                rm("modelno", envir = VGAMenv)
        } else {
            while(exists("CQO.FastAlgorithm"))
                remove("CQO.FastAlgorithm")
            while(exists("modelno"))
                remove("modelno")
        }
        dpar <- .dispersion
        if(!dpar) {
            temp87 = (y-mu)^2 * wz /
                     (dtheta.deta(mu, link= .link, earg = .earg )^2) # w cancel
            if(.mv && ! .onedpar) {
                dpar = rep(as.numeric(NA), len=M)
                temp87 = cbind(temp87)
                nrow.mu = if(is.matrix(mu)) nrow(mu) else length(mu)
                for(i in 1:M)
                    dpar[i] = sum(temp87[,i]) / (nrow.mu - ncol(x))
                if(is.matrix(y) && length(dimnames(y)[[2]])==length(dpar))
                    names(dpar) = dimnames(y)[[2]]
            } else 
                dpar = sum(temp87) / (length(mu) - ncol(x))
        }
        misc$mv = .mv
        misc$dispersion <- dpar
        misc$default.dispersion <- 1
        misc$estimated.dispersion <- .estimated.dispersion
        misc$link = rep( .link, length=M)
        names(misc$link) = if(M>1) dn2 else "mu"

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
        if(residuals) w*(y/mu - (1-y)/(1-mu)) else
            sum(w*(y*log(mu) + (1-y)*log(1-mu)))
    },
    vfamily=c("binomialff", "vcategorical"),
    deriv=eval(substitute(expression({
        if( .link == "logit") {
            w * (y - mu)
        } else if( .link == "cloglog") {
            mu.use = mu
            smallno = 100 * .Machine$double.eps
            mu.use[mu.use < smallno] = smallno
            mu.use[mu.use > 1 - smallno] = 1 - smallno
            -w * (y - mu) * log(1 - mu.use) / mu.use
        } else
            w * dtheta.deta(mu, link= .link, earg = .earg )* (y/mu - 1)/(1-mu)
    }), list( .link=link, .earg = earg ))),
    weight=eval(substitute(expression({
        tmp100 = mu*(1-mu)

        tmp200 = if( .link == "logit") {
            cbind(w * tmp100)
        } else if( .link == "cloglog") {
            cbind(w * (1-mu.use) * (log(1-mu.use))^2 / mu.use )
        } else {
            cbind(w * dtheta.deta(mu, link= .link, earg = .earg)^2 / tmp100)
        }
        for(ii in 1:M) {
            index200 = !is.finite(tmp200[,ii]) |
                       (abs(tmp200[,ii]) < .Machine$double.eps)
            if(any(index200)) {
                tmp200[index200,ii] = .Machine$double.eps # Diagonal 0's are bad 
            }
        }
        tmp200
    }), list( .link=link, .earg = earg ))))
}



gammaff <- function(link="nreciprocal", earg=list(),
                    dispersion=0)
{
    estimated.dispersion <- dispersion==0
    if(mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Gamma distribution\n\n",
           "Link:     ", namesof("mu", link, earg=earg), "\n",
           "Variance: mu^2 / k"),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        devi <- -2 * w * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu)
        if(residuals) {
            sign(y - mu) * sqrt(abs(devi) * w)
        } else sum(w * devi)
    },
    initialize=eval(substitute(expression({
        mustart <- y + 0.167 * (y == 0)
            M = if(is.matrix(y)) ncol(y) else 1
            dn2 = if(is.matrix(y)) dimnames(y)[[2]] else NULL
            dn2 = if(length(dn2)) {
                paste("E[", dn2, "]", sep="") 
            } else {
                paste("mu", 1:M, sep="") 
            }
            predictors.names = namesof(if(M>1) dn2 else "mu", .link,
                 earg=.earg, short=TRUE)
        if(!length(etastart))
            etastart <- theta2eta(mustart, link= .link, earg=.earg)
    }), list( .link=link, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta, link= .link, earg=.earg)
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        dpar <- .dispersion
        if(!dpar) {
            if(M == 1) {
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
        names(misc$link) = if(M>1) paste("mu", 1:M, sep="") else "mu"

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



inverse.gaussianff <- function(link="natural.ig", dispersion=0)
{
    estimated.dispersion <- dispersion==0
    warning("@deviance() not finished")
    warning("needs checking, but I'm sure it works")

    if(mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Inverse Gaussian distribution\n\n",
           "Link:     ", namesof("mu", link), "\n",
           "Variance: mu^3 /k"),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        pow <- 3  # Use Quasi()$deviance with pow==3
        devy  <- y^(2-pow) / (1-pow) - y^(2-pow) / (2-pow)
        devmu <- y * mu^(1-pow) / (1-pow) - mu^(2-pow) / (2-pow)
        devi <- 2 * (devy - devmu)
        if(residuals) {
            sign(y - mu) * sqrt(abs(devi) * w)
        } else sum(w * devi)
    },
    initialize=eval(substitute(expression({
        mu <- y + 0.167 * (y == 0)
        if(!length(etastart))
            etastart <- theta2eta(mu, link= .link)
    }), list( .link=link ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta, link= .link)
    }, list( .link=link ))),
    last=eval(substitute(expression({
        dpar <- .dispersion
        if(!dpar) {
            temp <- w * dmu.deta^2
            dpar <- sum( w * (y-mu)^2 * wz / temp ) / (length(mu) - ncol(x))
        }
        misc$dispersion <- dpar
        misc$default.dispersion <- 0
        misc$estimated.dispersion <- .estimated.dispersion
        misc$link = rep( .link, length=M)
        names(misc$link) = if(M>1) paste("mu", 1:M, sep="") else "mu"
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




dinv.gaussian = function(x, mu, lambda) {
    if(any(mu <=0)) stop("mu must be positive")
    if(any(lambda <=0)) stop("lambda must be positive")
    ans = x
    mu = rep(mu, len=length(x))
    lambda = rep(lambda, len=length(x))
    ans[x <= 0] = 0
    bb = x > 0
    ans[bb] = sqrt(lambda[bb]/(2*pi*x[bb]^3)) *
              exp(-lambda[bb]*(x[bb]-mu[bb])^2/(2*mu[bb]^2*x[bb]))
    ans
}


pinv.gaussian = function(q, mu, lambda) {
    if(any(mu <=0)) stop("mu must be positive")
    if(any(lambda <=0)) stop("lambda must be positive")
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
    if(!is.Numeric(n, positive=TRUE, integer=TRUE, allow=1))
        stop("'n' must be a single positive integer")
    if(!is.Numeric(mu, positive=TRUE))
        stop("'mu' must have positive values only")
    if(!is.Numeric(lambda, positive=TRUE))
        stop("'lambda' must have positive values only")
    mu = rep(mu, len=n)
    lambda = rep(lambda, len=n)
    u = runif(n)
    z = rnorm(n)^2
    phi = lambda / mu
    y1 = 1 - 0.5 * (sqrt(z^2 + 4*phi*z) - z) / phi
    mu * ifelse((1+y1)*u > 1, 1/y1, y1)
}



inv.gaussianff <- function(lmu="loge", llambda="loge",
                           emu=list(), elambda=list(),
                           ilambda=1,
                           zero=NULL)
{
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu <- as.character(substitute(lmu))
    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda <- as.character(substitute(llambda))
    if(!is.list(emu)) emu = list()
    if(!is.list(elambda)) elambda = list()

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
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(any(y <= 0)) stop("Require the response to have positive values")
        predictors.names =
        c(namesof("mu", .lmu, earg= .emu, short= TRUE),
          namesof("lambda", .llambda, earg= .elambda, short= TRUE))
        if(!length(etastart)) {
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
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w*(0.5 * log(lambda / (2 * pi * y^3)) -
                lambda *(y-mu)^2 / (2*mu^2 * y)))
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



poissonff <- function(link="loge", earg=list(),
                      dispersion=1, onedpar=FALSE,
                      parallel=FALSE, zero=NULL)
{

    estimated.dispersion <- dispersion==0
    if(mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Poisson distribution\n\n",
           "Link:     ", namesof("mu", link, earg= earg), "\n",
           "Variance: mu"),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(1,M,1), x, .parallel, constraints)
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero=zero ))),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        nz <- y > 0
        devi <-  - (y - mu)
        devi[nz] <- devi[nz] + y[nz] * log(y[nz]/mu[nz])
        if(residuals) sign(y - mu) * sqrt(2 * abs(devi) * w) else
            2 * sum(w * devi)
    },
    initialize=eval(substitute(expression({
        if(is.R()) assign("CQO.FastAlgorithm",
            ( .link == "loge"), envir = VGAMenv) else
            CQO.FastAlgorithm <<- ( .link == "loge")
        M = if(is.matrix(y)) ncol(y) else 1
        dn2 = if(is.matrix(y)) dimnames(y)[[2]] else NULL
        dn2 = if(length(dn2)) {
            paste("E[", dn2, "]", sep="") 
        } else {
            paste("mu", 1:M, sep="") 
        }
        predictors.names = namesof(if(M>1) dn2 else "mu", .link,
            earg= .earg, short=TRUE)
        mu = pmax(y, 1/8) # y + 0.167 * (y == 0)
        if(!length(etastart))
            etastart <- theta2eta(mu, link= .link, earg= .earg)
    }), list( .link=link, .estimated.dispersion=estimated.dispersion,
              .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mu = eta2theta(eta, link= .link, earg= .earg)
        mu
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
       if(is.R()) {
            if(exists("CQO.FastAlgorithm", envir = VGAMenv))
                rm("CQO.FastAlgorithm", envir = VGAMenv)
        } else {
            while(exists("CQO.FastAlgorithm"))
                remove("CQO.FastAlgorithm")
        }
        dpar <- .dispersion
        if(!dpar) {
            temp87 = (y-mu)^2 *
                wz / (dtheta.deta(mu, link= .link, earg= .earg)^2) # w cancel
            if(M > 1 && ! .onedpar) {
                dpar = rep(as.numeric(NA), len=M)
                temp87 = cbind(temp87)
                nrow.mu = if(is.matrix(mu)) nrow(mu) else length(mu)
                for(i in 1:M)
                    dpar[i] = sum(temp87[,i]) / (nrow.mu - ncol(x))
                if(is.matrix(y) && length(dimnames(y)[[2]])==length(dpar))
                    names(dpar) = dimnames(y)[[2]]
            } else 
                dpar = sum(temp87) / (length(mu) - ncol(x))
        }
        misc$dispersion <- dpar
        misc$default.dispersion <- 1
        misc$estimated.dispersion <- .estimated.dispersion
        misc$expected = TRUE
        misc$link = rep( .link, length=M)
        names(misc$link) = if(M>1) dn2 else "mu"

        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M) misc$earg[[ii]] = .earg
    }), list( .dispersion=dispersion, .estimated.dispersion=estimated.dispersion,
            .onedpar=onedpar, .link=link, .earg=earg ))),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, link= .link, earg= .earg)
    }, list( .link=link, .earg=earg ))),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        if(residuals) w*(y/mu - 1) else sum(w*(-mu + y*log(mu) - lgamma(y+1)))
    },
    vfamily="poissonff",
    deriv=eval(substitute(expression({
        if( .link == "loge" && (any(mu < .Machine$double.eps))) {
            w * (y - mu)
        } else {
            lambda <- mu
            dl.dlambda <- (y-lambda) / lambda
            dlambda.deta <- dtheta.deta(theta=lambda, link= .link, earg= .earg)
            w * dl.dlambda * dlambda.deta
        }
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        if( .link == "loge" && (any(mu < .Machine$double.eps))) {
            tmp600 = mu
            tmp600[tmp600 < .Machine$double.eps] = .Machine$double.eps
            w * tmp600
        } else {
            d2l.dlambda2 = 1 / lambda
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


poissonqn <- function(link="loge", earg=list(),
                      dispersion=1, onedpar=FALSE,
                      parallel=FALSE, zero=NULL,
                      wwts=c("expected","observed","qn"))
{
    estimated.dispersion <- dispersion==0
    if(mode(link )!= "character" && mode(link )!= "name")
        link <- as.character(substitute(link))
    if(mode(wwts) != "character" && mode(wwts) != "name")
        wwts <- as.character(substitute(wwts))
    wwts <- match.arg(wwts, c("expected","observed","qn"))[1]
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Poisson distribution\n\n",
           "Link:     ", namesof("mu", link, earg= earg), "\n",
           "Variance: mu"),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(matrix(1,M,1), x, .parallel, constraints)
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .parallel=parallel, .zero=zero ))),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        nz <- y > 0
        devi <-  - (y - mu)
        devi[nz] <- devi[nz] + y[nz] * log(y[nz]/mu[nz])
        if(residuals) sign(y - mu) * sqrt(2 * abs(devi) * w) else
            2 * sum(w * devi)
    },
    initialize=eval(substitute(expression({
        M = if(is.matrix(y)) ncol(y) else 1
        dn2 = if(is.matrix(y)) dimnames(y)[[2]] else NULL
        dn2 = if(length(dn2)) {
            paste("E[", dn2, "]", sep="") 
        } else {
            paste("mu", 1:M, sep="") 
        }
        predictors.names = namesof(if(M>1) dn2 else "mu", .link,
            earg= .earg, short=TRUE)
        mu = pmax(y, 0.167)  # y + 0.167 * (y == 0)
        if(!length(etastart))
            etastart <- theta2eta(mu, link= .link, earg= .earg)
    }), list( .link=link, .estimated.dispersion=estimated.dispersion,
              .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta, link= .link, earg= .earg)
    }, list( .link=link,
              .earg=earg ))),
    last=eval(substitute(expression({
        dpar <- .dispersion
        if(!dpar) {
            temp87 = (y-mu)^2 * wz / (dtheta.deta(mu, link= .link, earg= .earg)^2) # w cancel
            if(M > 1 && ! .onedpar) {
                dpar = rep(as.numeric(NA), len=M)
                temp87 = cbind(temp87)
                nrow.mu = if(is.matrix(mu)) nrow(mu) else length(mu)
                for(i in 1:M)
                    dpar[i] = sum(temp87[,i]) / (nrow.mu - ncol(x))
                if(is.matrix(y) && length(dimnames(y)[[2]])==length(dpar))
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
        names(misc$link) = if(M>1) dn2 else "mu"

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
        if(residuals) w*(y/mu - 1) else sum(w*(-mu + y*log(mu) - lgamma(y+1)))
    },
    vfamily="poissonqn",
    deriv=eval(substitute(expression({
        if(iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }

        derivnew =
        if( .link == "loge" && (any(mu < .Machine$double.eps))) {
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
        if( .wwts == "qn") {
            if(iter == 1) {
                wznew = cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
            } else {
                wzold = wznew
                wznew = qnupdate(w=w, wzold=wzold, dderiv=(derivold-derivnew),
                                 deta=etanew-etaold, M=M,
                                 trace=trace)  # weights incorporated in args
            }
        } else if( .wwts == "expected") {
            wznew = if( .link == "loge") {
                tmp600 = mu
                tmp600[tmp600 < .Machine$double.eps] = .Machine$double.eps
                w * tmp600
            } else {
                d2l.dlambda2 = 1 / lambda
                w * dlambda.deta^2 * d2l.dlambda2
            }
        } else {
            wznew = if( .link == "loge") {
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


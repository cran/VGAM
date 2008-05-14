# These functions are
# Copyright (C) 1998-2008 T.W. Yee, University of Auckland. All rights reserved.





gev <- function(llocation="identity",
                lscale="loge",
                lshape="logoff",
                elocation = list(),
                escale = list(),
                eshape = if(lshape=="logoff") list(offset=0.5) else 
                if(lshape=="elogit") list(min=-0.5, max=0.5) else list(),
                percentiles=c(95,99),
                iscale=NULL, ishape=NULL,
                method.init=1, gshape=c(-0.45, 0.45), tshape0=0.001,
                zero=3)
{


    mean = FALSE
    if(length(iscale) && !is.Numeric(iscale, posit=TRUE))
        stop("bad input for argument \"iscale\"")
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(!mean &&  length(percentiles) &&
       (!is.Numeric(percentiles, posit=TRUE) || max(percentiles) >= 100))
        stop("bad input for argument \"percentiles\"")
    if(!is.Numeric(method.init, allow=1, posit=TRUE, integer=TRUE) ||
       method.init > 2.5)
        stop("argument \"method.init\" must be 1 or 2")
    if(length(ishape) && !is.Numeric(ishape))
        stop("bad input for argument \"ishape\"")
    if(!is.Numeric(tshape0, allow=1, posit=TRUE) || tshape0 > 0.1)
        stop("bad input for argument \"tshape0\"")
    if(!is.Numeric(gshape, allow=2) || gshape[1] >= gshape[2])
        stop("bad input for argument \"gshape\"")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(elocation)) elocation = list()
    if(!is.list(escale)) escale = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("Generalized Extreme Value Distribution\n",
            "Links:    ",
            namesof("location", link=llocation, earg= elocation), ", ", 
            namesof("scale", link=lscale, earg= escale), ", ",
            namesof("shape", link=lshape, earg= eshape)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names = 
        c(namesof("location", .llocation, earg= .elocation, short=TRUE),
          namesof("scale", .lscale, earg= .escale, short=TRUE),
          namesof("shape", .lshape, earg= .eshape, short=TRUE))
        y = as.matrix(y)
        if(ncol(y) > 1)
            y = -t(apply(-y, 1, sort, na.last=TRUE))

        sum1 = function(mat) {
            if(!is.matrix(mat)) mat = as.matrix(mat)
            apply(mat, 1, sum, na.rm=TRUE)
        }
        r.vec = sum1(!is.na(y))
        if(any(r.vec == 0))
            stop("A row contains all missing values")

        extra$percentiles = .percentiles
        if(!length(etastart)) {
            siginit= if(length( .iscale)) rep( .iscale, len=nrow(y)) else NULL
            xiinit = if(length( .ishape)) rep( .ishape, len=nrow(y)) else NULL
            eshape = .eshape
            if( .lshape=="elogit" && length(xiinit) &&
                (any(xiinit <= eshape$min | xiinit >= eshape$max)))
                stop("bad input for argument \"eshape\"")
            if( .method.init == 1) {
                nvector = 4:10   # Arbitrary; could be made an argument
                ynvector = quantile(y[,1], probs = 1-1/nvector)
                objecFunction = -Inf   # Actually the log-likelihood
                est.sigma = !length(siginit)
                gshape = .gshape
                temp234 = if(length(xiinit)) xiinit[1] else
                          seq(gshape[1], gshape[2], len=12)
                for(xi.try in temp234) {
                    xvec = if(abs(xi.try) < .tshape0) log(nvector) else
                           (nvector^xi.try - 1) / xi.try
                    fit0 = lsfit(x=xvec, y=ynvector, intercept=TRUE)
                    if(est.sigma) {
                        sigmaTry = rep(fit0$coef["X"], len=nrow(y))
                    } else { 
                        sigmaTry = siginit
                    }
                    muTry = rep(fit0$coef["Intercept"], len=nrow(y))
                    llTry = egev()@loglikelihood(mu=NULL,y=y[,1],w=w,
                            residuals=FALSE, extra=list(giveWarning=FALSE),
                            eta=cbind(theta2eta(muTry,  .llocation, earg= .elocation),
                                      theta2eta(sigmaTry, .lscale, earg= .escale), 
                                      theta2eta(xi.try, link= .lshape,
                                                earg= .eshape)))
                    if(llTry >= objecFunction) {
                        if(est.sigma)
                            siginit = sigmaTry
                        muinit = rep(muTry, len=nrow(y))
                        objecFunction = llTry
                        bestxi = xi.try
                    }
                }
                if(!length(xiinit))
                    xiinit = rep(bestxi, len=nrow(y))
            } else {
                xiinit = rep(0.05, len=nrow(y))
                if(!length(siginit))
                    siginit = rep(sqrt(6 * var(y[,1]))/pi, len=nrow(y))
                EulerM = -digamma(1)
                muinit = rep(median(y[,1]) - EulerM*siginit, len=nrow(y))
            }

            bad = ((1 + xiinit*(y-muinit)/siginit) <= 0)
            if(fred <- sum(sum1(bad))) {
                warning(paste(fred, "observations violating boundary",
                "constraints while initializing. Taking corrective action."))
                xiinit[bad] = ifelse(y[bad] > muinit[bad], 0.1, -0.1)
            }

            etastart = cbind(theta2eta(muinit,  .llocation, earg= .elocation),
                             theta2eta(siginit, .lscale, earg= .escale),
                             theta2eta(xiinit,  .lshape, earg= .eshape))
        }
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .percentiles=percentiles,
              .elocation = elocation, .escale = escale,
              .eshape= eshape, .tshape0=tshape0,
              .method.init=method.init,
              .iscale=iscale, .ishape=ishape, .gshape=gshape ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        loc = eta2theta(eta[,1], .llocation, earg= .elocation)
        sigma = eta2theta(eta[,2], .lscale, earg= .escale)
        xi = eta2theta(eta[,3], .lshape, earg= .eshape)
        iszero = (abs(xi) < .tshape0)
        cent = extra$percentiles
        lp = length(cent)
        fv = matrix(as.numeric(NA), nrow(eta), lp)
        if(lp) {
            for(i in 1:lp) {
                yp = -log(cent[i]/100)
                fv[!iszero,i] = loc[!iszero] - sigma[!iszero] *
                                (1 - yp^(-xi[!iszero])) / xi[!iszero]
                fv[iszero,i] = loc[iszero] - sigma[iszero] * log(yp)
            }
            dimnames(fv) = list(dimnames(eta)[[1]],
                                paste(as.character(cent), "%", sep=""))
        } else {
            EulerM = -digamma(1)
            fv = loc + sigma * EulerM  # When xi=0, is Gumbel
            fv[!iszero] = loc[!iszero] + sigma[!iszero] *
                          (gamma(1-xi[!iszero])-1) / xi[!iszero]
            fv[xi >= 1] = NA  # Mean exists only if xi < 1.
        }
        fv
    }, list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
             .eshape= eshape, .tshape0=tshape0 ))),
    last=eval(substitute(expression({
        misc$links = c(location = .llocation, scale = .lscale, shape =.lshape)
        misc$true.mu = !length( .percentiles) # @fitted is not a true mu
        misc$percentiles = .percentiles
        misc$earg= list(location= .elocation, scale= .escale, shape= .eshape)
        misc$expected = TRUE
        misc$tshape0 = .tshape0
        if(ncol(y)==1)
            y = as.vector(y)
        if(any(xi < -0.5))
            warning("some values of the shape parameter are less than -0.5")
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .elocation = elocation, .escale = escale, .eshape= eshape,
              .tshape0=tshape0, .percentiles=percentiles ))),
    loglikelihood=eval(substitute(
      function(mu,y,w,residuals=FALSE,eta,extra=NULL) {
            mmu = eta2theta(eta[,1], .llocation, earg= .elocation)
            sigma = eta2theta(eta[,2], .lscale, earg= .escale)
            xi = eta2theta(eta[,3], .lshape, earg= .eshape)
            iszero = (abs(xi) < .tshape0)
            sum1 = function(mat) {
                if(!is.matrix(mat)) mat = as.matrix(mat)
                apply(mat, 1, sum, na.rm=TRUE)
            }
            zedd = (y-mmu) / sigma
            r.vec = sum1(!is.na(y))
            A = 1 + xi * (y-mmu)/sigma
            ii = 1:nrow(eta)
            A1 = A[cbind(ii, r.vec)]
            mytolerance = 0  # .Machine$double.eps
            if(any(bad <- (A1 <= mytolerance), na.rm=TRUE)) {
                cat("There are",sum(bad),"range violations in @loglikelihood\n")
                if(exists("flush.console")) flush.console()
            }
            igev = !iszero &  !bad
            igum =  iszero &  !bad
            pow = 1 + 1/xi[igev]
        if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(bad) * (-1.0e10) +
            sum(w[igum] * (-r.vec[igum]*log(sigma[igum]) -
                           exp(-zedd[igum,r.vec]) - sum1(zedd))) +
            sum(w[igev] * (-r.vec[igev]*log(sigma[igev]) -
                           pow*sum1(log(A[igev])) -
                           A1[igev]^(-1/xi[igev])))
        }, list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
                 .elocation = elocation, .escale = escale, .eshape= eshape,
                 .tshape0=tshape0 ))),
    vfamily=c("gev", "vextremes"),
    deriv=eval(substitute(expression({
        sum1 = function(mat) {
            if(!is.matrix(mat)) mat = as.matrix(mat)
            apply(mat, 1, sum, na.rm=TRUE)
        }
        r.vec = sum1(!is.na(y))
        mmu = eta2theta(eta[,1], .llocation, earg= .elocation)
        sigma = eta2theta(eta[,2], .lscale, earg= .escale)
        xi = eta2theta(eta[,3], .lshape, earg= .eshape)
        iszero = (abs(xi) < .tshape0)
        ii = 1:nrow(eta)
        zedd = (y-mmu) / sigma
        A = 1 + xi * zedd
        dA.dxi = zedd                # matrix
        dA.dmu = -xi/sigma           # vector
        dA.dsigma = -xi*zedd/sigma   # matrix
        pow = 1 + 1/xi
        A1 = A[cbind(ii, r.vec)]

        AAr1 = dA.dmu/(xi * A1^pow) - pow * sum1(dA.dmu/A)
        AAr2 = dA.dsigma[cbind(ii,r.vec)] / (xi * A1^pow) -
               pow * sum1(dA.dsigma/A)
        AAr3 = 1/(xi * A1^pow) - pow * sum1(dA.dsigma/A)
        dl.dmu = AAr1
        dl.dsi = AAr2 - r.vec/sigma
        dl.dxi = sum1(log(A))/xi^2 - pow * sum1(dA.dxi/A) - (log(A1) / 
                xi^2 - dA.dxi[cbind(ii,r.vec)] /(xi*A1)) * A1^(-1/xi)

        if(any(iszero)) {
            zorro = c(zedd[cbind(1:n,r.vec)])
            zorro = zorro[iszero]
            ezedd = exp(-zorro)
            dl.dmu[iszero] = (1-ezedd) / sigma[iszero]
            dl.dsi[iszero] = (zorro * (1-ezedd) - 1) / sigma[iszero]
            dl.dxi[iszero] = zorro * ((1 - ezedd) * zorro / 2 - 1)
        }
        dmu.deta = dtheta.deta(mmu, .llocation, earg= .elocation)
        dsi.deta = dtheta.deta(sigma, .lscale, earg= .escale)
        dxi.deta = dtheta.deta(xi, .lshape, earg= .eshape)
        w * cbind(dl.dmu * dmu.deta, dl.dsi * dsi.deta, dl.dxi * dxi.deta)
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .elocation = elocation, .escale = escale, .eshape= eshape,
              .tshape0=tshape0 ))),
    weight=eval(substitute(expression({
        kay = -xi
        dd = digamma(r.vec-kay+1)
        ddd = digamma(r.vec+1) # Unnecessarily evaluated at each iteration
        temp13 = -kay * dd + (kay^2 - kay + 1) / (1-kay)
        temp33 = 1 - 2 * kay * ddd + kay^2 * (1 + trigamma(r.vec+1) + ddd^2)
        temp23 = -kay * dd + (1+(1-kay)^2) / (1-kay)
        GR.gev = function(j, ri, kay) gamma(ri - j*kay + 1) /  gamma(ri)
        tmp2 = (1-kay)^2 * GR.gev(2, r.vec, kay)  # Latter is GR2
        tmp1 = (1-2*kay) * GR.gev(1, r.vec, kay)  # Latter is GR1
        k0 = (1-2*kay)
        k1 = k0 * kay
        k2 = k1 * kay
        k3 = k2 * kay   # kay^3 * (1-2*kay)
        wz = matrix(as.numeric(NA), n, 6)
        wz[,iam(1,1,M)] = tmp2 / (sigma^2 * k0)
        wz[,iam(1,2,M)] = (tmp2 - tmp1) / (sigma^2 * k1)
        wz[,iam(1,3,M)] = (tmp1 * temp13 - tmp2) / (sigma * k2)
        wz[,iam(2,2,M)] = (r.vec*k0 - 2*tmp1 + tmp2) / (sigma^2 * k2)
        wz[,iam(2,3,M)] = (r.vec*k1*ddd + tmp1 *
                           temp23 - tmp2 - r.vec*k0) / (sigma * k3)
        wz[,iam(3,3,M)] = (2*tmp1*(-temp13) + tmp2 + r.vec*k0*temp33)/(k3*kay)

        if(any(iszero)) {
            if(ncol(y) > 1)
                stop("cannot handle xi==0 with a multivariate response")

            EulerM = -digamma(1)
            wz[iszero,iam(2,2,M)] = (pi^2/6 + (1-EulerM)^2) / sigma^2
            wz[iszero,iam(3,3,M)] = 2.4236
            wz[iszero,iam(1,2,M)] = (digamma(2) + 2*(EulerM-1)) / sigma^2
            wz[iszero,iam(1,3,M)]= -(trigamma(1)/2 + digamma(1)*
                                    (digamma(1)/2+1))/sigma
            wz[iszero,iam(2,3,M)] = (-dgammadx(2,3)/6 + dgammadx(1,1) +
                                    2*dgammadx(1,2) + 2*dgammadx(1,3)/3) / sigma

            if(FALSE ) {
            wz[,iam(1,2,M)] = 2 * r.vec / sigma^2
            wz[,iam(2,2,M)] = -4 * r.vec * digamma(r.vec+1) + 2 * r.vec +
                (4 * dgammadx(r.vec+1, der=1) - 
                 3 * dgammadx(r.vec+1, der=2)) / gamma(r.vec) # Not checked
            }
        }
        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] * dmu.deta^2
        wz[,iam(2,2,M)] = wz[,iam(2,2,M)] * dsi.deta^2
        wz[,iam(3,3,M)] = wz[,iam(3,3,M)] * dxi.deta^2
        wz[,iam(1,2,M)] = wz[,iam(1,2,M)] * dmu.deta * dsi.deta
        wz[,iam(1,3,M)] = wz[,iam(1,3,M)] * dmu.deta * (-dxi.deta)
        wz[,iam(2,3,M)] = wz[,iam(2,3,M)] * dsi.deta * (-dxi.deta)
        w * wz
    }), list( .eshape = eshape ))))
}





dgammadx = function(x, deriv.arg=1) {
    if(deriv.arg==0) {
        gamma(x)
    } else if(deriv.arg == 1) {
        digamma(x) * gamma(x)
    } else if(deriv.arg == 2) {
        gamma(x) * (trigamma(x) + digamma(x)^2)
    } else if(deriv.arg == 3) {
        gamma(x) * (psigamma(x, der=2) + 2 * digamma(x) * trigamma(x)) +
        dgammadx(x, der=1) * (trigamma(x) + digamma(x)^2)
    } else if(deriv.arg == 4) {
        dgammadx(x, der=2) * (trigamma(x) + digamma(x)^2) +
    2 * dgammadx(x, der=1) * (psigamma(x, der=2) + 2*digamma(x) * trigamma(x)) +
        gamma(x) * (psigamma(x, der=3) + 2*trigamma(x)^2 +
                    2 * digamma(x) * psigamma(x, der=2))
    } else stop("cannot handle deriv>4")
}



egev <- function(llocation="identity",
                 lscale="loge",
                 lshape="logoff",
                elocation = list(),
                escale = list(),
                eshape = if(lshape=="logoff") list(offset=0.5) else 
                if(lshape=="elogit") list(min=-0.5, max=0.5) else list(),
                 percentiles=c(95,99),
                 iscale=NULL, ishape=NULL,
                 method.init=1, gshape=c(-0.45, 0.45),
                 tshape0=0.001,
                 zero=3)
{
    if(length(iscale) && !is.Numeric(iscale, posit=TRUE))
        stop("bad input for argument \"iscale\"")
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale <- as.character(substitute(lscale))
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation <- as.character(substitute(llocation))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape <- as.character(substitute(lshape))
    if(!is.Numeric(gshape, allow=2) || gshape[1] >= gshape[2])
        stop("bad input for argument \"gshape\"")
    if(length(percentiles) && 
      (!is.Numeric(percentiles, posit=TRUE) || max(percentiles) >= 100))
        stop("bad input for argument \"percentiles\"")
    if(!is.Numeric(method.init, allow=1, posit=TRUE, integer=TRUE) ||
       method.init > 2.5)
        stop("argument \"method.init\" must be 1 or 2")
    if(length(ishape) && !is.Numeric(ishape))
        stop("bad input for argument \"ishape\"")
    if(!is.Numeric(tshape0, allow=1, posit=TRUE) || tshape0 > 0.1)
        stop("bad input for argument \"tshape0\"")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(elocation)) elocation = list()
    if(!is.list(escale)) escale = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("Generalized Extreme Value Distribution\n",
            "Links:    ",
            namesof("location", link=llocation, earg= elocation), ", ", 
            namesof("scale", link=lscale, earg= escale), ", ",
            namesof("shape", link=lshape, earg= eshape)),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names =
           c(namesof("location", .llocation, earg= .elocation,  short=TRUE),
             namesof("scale", .lscale, earg= .escale,  short=TRUE),
             namesof("shape",  .lshape, earg= .eshape, short=TRUE))
        if(ncol(as.matrix(y)) != 1)
            stop("response must be a vector or one-column matrix")
        if(!length(etastart)) {
            siginit= if(length( .iscale)) rep( .iscale, len=length(y)) else NULL
            xiinit = if(length( .ishape)) rep( .ishape, len=length(y)) else NULL
            eshape = .eshape
            if( .lshape == "elogit" && length(xiinit) && 
               (any(xiinit <= eshape$min | xiinit >= eshape$max)))
                stop("bad input for argument \"eshape\"")
            if( .method.init == 1) {
                nvector = 4:10   # Arbitrary; could be made an argument
                ynvector = quantile(y, probs = 1-1/nvector)
                objecFunction = -Inf   # Actually the log-likelihood
                est.sigma = !length(siginit)
                gshape = .gshape
                temp234 = if(length(xiinit)) xiinit[1] else
                          seq(gshape[1], gshape[2], len=12)
                for(xi.try in temp234) {
                    xvec = if(abs(xi.try) < .tshape0) log(nvector) else
                           (nvector^xi.try - 1) / xi.try
                    fit0 = lsfit(x=xvec, y=ynvector, intercept=TRUE)
                    if(est.sigma) {
                        sigmaTry = rep(fit0$coef["X"], len=length(y))
                    } else { 
                        sigmaTry = siginit
                    }
                    muTry = rep(fit0$coef["Intercept"], len=length(y))
                    llTry = egev()@loglikelihood(mu=NULL,y=y,w=w,
                            residuals=FALSE, extra=list(giveWarning=FALSE),
                        eta=cbind(theta2eta(muTry, .llocation,earg= .elocation),
                                  theta2eta(sigmaTry, .lscale, earg= .escale), 
                                  theta2eta(xi.try,  .lshape, earg= .eshape)))
                    if(llTry >= objecFunction) {
                        if(est.sigma)
                            siginit = sigmaTry
                        muinit = rep(muTry, len=length(y))
                        objecFunction = llTry
                        bestxi = xi.try
                    }
                }
                if(!length(xiinit))
                    xiinit = rep(bestxi, len=length(y))

            } else {
                xiinit = rep(if(length(xiinit)) xiinit else 0.05, len=length(y))
                if(!length(siginit))
                    siginit = rep(sqrt(6*var(y))/pi, len=length(y))
                EulerM = -digamma(1)
                muinit = rep(median(y) - EulerM * siginit, len=length(y))
            }
            bad <- (1 + xiinit*(y-muinit)/siginit <= 0)
            if(fred <- sum(bad, na.rm=TRUE)) {
                warning(paste(fred, "observations violating boundary",
                "constraints while initializing. Taking corrective action."))
                xiinit[bad] = ifelse(y[bad] > muinit[bad], 0.01, -0.01)
            }

            extra$percentiles = .percentiles

            etastart = cbind(theta2eta(muinit,  .llocation, earg= .elocation),
                             theta2eta(siginit, .lscale, earg= .escale), 
                             theta2eta(xiinit,  .lshape, earg= .eshape))
        }
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .percentiles=percentiles, .tshape0=tshape0,
              .elocation = elocation, .escale = escale, .eshape= eshape,
              .method.init=method.init,
              .iscale=iscale, .ishape=ishape, .gshape=gshape ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        loc <- eta2theta(eta[,1], .llocation, earg= .elocation)
        sigma <- eta2theta(eta[,2], .lscale, earg= .escale)
        xi <- eta2theta(eta[,3], .lshape, earg= .eshape)
        iszero <- (abs(xi) < .tshape0)
        cent = extra$percentiles
        lp <- length(cent)
        fv <- matrix(as.numeric(NA), nrow(eta), lp)
        if(lp) {
            for(i in 1:lp) {
                yp = -log(cent[i]/100)
                fv[!iszero,i] = loc[!iszero] - sigma[!iszero] *
                                (1 - yp^(-xi[!iszero])) / xi[!iszero]
                fv[iszero,i] = loc[iszero] - sigma[iszero] * log(yp)
            }
            dimnames(fv) = list(dimnames(eta)[[1]],
                                paste(as.character(cent), "%", sep=""))
        } else {
            EulerM = -digamma(1)
            fv = loc + sigma * EulerM  # When xi=0, is Gumbel
            fv[!iszero] = loc[!iszero] + sigma[!iszero] *
                          (gamma(1-xi[!iszero])-1) / xi[!iszero]
            fv[xi >= 1] = NA  # Mean exists only if xi < 1.
        }
        fv
    }, list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
             .elocation = elocation, .escale = escale, .eshape= eshape,
             .tshape0=tshape0 ))),
    last=eval(substitute(expression({
        misc$links <- c(location = .llocation, scale = .lscale, shape = .lshape)
        misc$true.mu = !length( .percentiles) # @fitted is not a true mu
        misc$percentiles <- .percentiles
        misc$earg= list(location= .elocation, scale= .escale, shape= .eshape)
        misc$tshape0 = .tshape0
        misc$expected = TRUE 
        if(any(xi < -0.5))
            warning("some values of the shape parameter are less than -0.5")
      }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
                .elocation = elocation, .escale = escale, .eshape= eshape,
                .tshape0=tshape0,  .percentiles=percentiles ))),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals=FALSE, eta, extra=NULL) {
            giveWarning = if(length(extra$giveWarning)) extra$giveWarning else
                          FALSE
            mmu <- eta2theta(eta[,1], .llocation, earg= .elocation )
            sigma <- eta2theta(eta[,2], .lscale, earg= .escale )
            xi <- eta2theta(eta[,3], .lshape, earg= .eshape )
            zedd <- (y-mmu) / sigma
            A <- 1 + xi * zedd
            if(any(iszero <- (abs(xi) < .tshape0))) {
            }
            mytolerance = 0  # .Machine$double.eps
            bad <- (A<=mytolerance)   # Range violation
            if(any(bad) && giveWarning) {
                cat("There are",sum(bad),"range violations in @loglikelihood\n")
                if(exists("flush.console")) flush.console()
            }
            igev = !iszero &  !bad
            igum =  iszero &  !bad
        if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(bad) * (-1.0e10) +
            sum(w[igum] * (-log(sigma[igum]) - zedd[igum] - exp(-zedd[igum]))) +
            sum(w[igev] * (-log(sigma[igev]) - (1 + 1/xi[igev])*log(A[igev]) -
                           A[igev]^(-1/xi[igev])))
        }, list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
                 .elocation = elocation, .escale = escale, .eshape= eshape,
                 .tshape0=tshape0 ))),
    vfamily=c("egev", "vextremes"),
    deriv=eval(substitute(expression({
        mmu = eta2theta(eta[,1], .llocation, earg= .elocation)
        sigma = eta2theta(eta[,2], .lscale, earg= .escale )
        xi = eta2theta(eta[,3], .lshape, earg= .eshape)
        iszero <- (abs(xi) < .tshape0)
        zedd = (y-mmu) / sigma
        A = 1 + xi * zedd
        dA.dxi = zedd
        dA.dmu = -xi / sigma
        dA.dsigma = -xi * zedd / sigma
        pow = 1 + 1/xi
        if(any(bad <- A<=0, na.rm=TRUE)) stop(paste(sum(bad, na.rm=TRUE),
            "observations violating boundary constraints in @deriv"))
        AA = 1/(xi*A^pow)- pow/A 
        dl.dmu = dA.dmu * AA
        dl.dsi = dA.dsigma * AA - 1/sigma
        dl.dxi =  log(A)/xi^2 - pow * dA.dxi / A -
               (log(A)/xi^2 - dA.dxi /(xi*A)) * A^(-1/xi)
        if(any(iszero)) {
            ezedd = exp(-zedd[iszero])
            dl.dmu[iszero] = (1-ezedd) / sigma[iszero]
            dl.dsi[iszero] = (zedd[iszero] * (1-ezedd) - 1) / sigma[iszero]
            dl.dxi[iszero] = zedd[iszero] * ((1 - ezedd) * zedd[iszero] / 2 - 1)
        }
        dmu.deta = dtheta.deta(mmu, .llocation, earg= .elocation)
        dsi.deta = dtheta.deta(sigma, .lscale, earg= .escale )
        dxi.deta = dtheta.deta(xi, .lshape, earg= .eshape)
        w * cbind(dl.dmu * dmu.deta, dl.dsi * dsi.deta, dl.dxi*dxi.deta)
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
              .elocation = elocation, .escale = escale, .eshape= eshape,
              .tshape0=tshape0 ))),
    weight=eval(substitute(expression({
        bad <- A <= 0
        if(any(bad, na.rm = TRUE)) stop(paste(sum(bad, na.rm = TRUE),
             "observations violating boundary constraints in @weight"))
        kay = -xi  # for the formulae 
        kay[abs(kay-0.5) < .tshape0] = 0.501
        temp100 = gamma(2-kay)
        pp = (1-kay)^2 * gamma(1-2*kay)    # gamma(0) is undefined so kay != 0.5
        qq = temp100 * (digamma(1-kay) - (1-kay)/kay)
        wz = matrix(as.numeric(NA), n, 6)
        wz[,iam(1,1,M)] = pp / sigma^2
        wz[,iam(2,2,M)] = (1-2*temp100 + pp) / (sigma * kay)^2
        EulerM = -digamma(1)
        wz[,iam(3,3,M)] = (pi^2 / 6 + (1-EulerM-1/kay)^2 +
                           (2*qq + pp/kay)/kay) / kay^2 
        wz[,iam(1,2,M)] = (pp - temp100) / (sigma^2 * kay)
        wz[,iam(1,3,M)] = -(qq + pp/kay) / (sigma * kay)
        wz[,iam(2,3,M)] = (1-EulerM - (1-temp100)/kay - qq -
                            pp/kay) / (sigma * kay^2)
        if(any(iszero)) {
            wz[iszero,iam(2,2,M)] = (pi^2/6 + (1-EulerM)^2) / sigma^2
            wz[iszero,iam(3,3,M)] = 2.4236
            wz[iszero,iam(1,2,M)] = (digamma(2) + 2*(EulerM-1)) / sigma^2
            wz[iszero,iam(1,3,M)]= -(trigamma(1)/2 + digamma(1)*
                                    (digamma(1)/2+1))/sigma
            wz[iszero,iam(2,3,M)] = (-dgammadx(2,3)/6 + dgammadx(1,1) +
                                    2*dgammadx(1,2) + 2*dgammadx(1,3)/3) / sigma
        }
        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] * dmu.deta^2
        wz[,iam(2,2,M)] = wz[,iam(2,2,M)] * dsi.deta^2
        wz[,iam(3,3,M)] = wz[,iam(3,3,M)] * dxi.deta^2
        wz[,iam(1,2,M)] = wz[,iam(1,2,M)] * dmu.deta * dsi.deta
        wz[,iam(1,3,M)] = wz[,iam(1,3,M)] * dmu.deta * (-dxi.deta)
        wz[,iam(2,3,M)] = wz[,iam(2,3,M)] * dsi.deta * (-dxi.deta)
        w * wz
    }), list( .eshape= eshape, .tshape0=tshape0 ))))
}






rgumbel = function(n, location=0, scale=1) {
    if(!is.Numeric(n, posit=TRUE, allow=1, integ=TRUE)) 
        stop("bad input for argument \"n\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    location - scale * log(-log(runif(n)))
}

dgumbel = function(x, location=0, scale=1) {
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    temp = exp(-(x-location) / scale)
    temp * exp(-temp) / scale
}

qgumbel = function(p, location=0, scale=1) {
    if(!is.Numeric(p, posit=TRUE) || any(p >= 1)) stop("0 < p < 1 is required")
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    location - scale * log(-log(p))
}

pgumbel = function(q, location=0, scale=1) {
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    exp(-exp(-(q-location) / scale))
}


gumbel <- function(llocation="identity",
                   lscale="loge",
                   elocation = list(),
                   escale = list(),
                   iscale=NULL,
                   R=NA, percentiles=c(95,99),
                   mpv=FALSE, zero=NULL)
{
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(!is.logical(mpv) || length(mpv) != 1)
        stop("bad input for argument \"mpv\"")
    if(length(percentiles) &&
       (!is.Numeric(percentiles, posit=TRUE) || max(percentiles) >= 100))
        stop("bad input for argument \"percentiles\"")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(length(iscale) && !is.Numeric(iscale, posit=TRUE))
        stop("bad input for argument \"iscale\"")
    if(!is.list(elocation)) elocation = list()
    if(!is.list(escale)) escale = list()

    new("vglmff",
    blurb=c("Gumbel Distribution for Extreme Value Regression\n",
           "Links:    ",
           namesof("location", link=llocation, earg= elocation), ", ",
           namesof("scale", link=lscale, earg= escale )),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names = 
        c(namesof("location", .llocation, earg= .elocation, short=TRUE),
          namesof("scale", .lscale, earg= .escale , short=TRUE))
        y = as.matrix(y)
        if(ncol(y) > 1)
        y = if(is.R()) -t(apply(-y,1,sort, na.last = TRUE)) else {
            y[is.na(y)] = -Inf
            y = -t(apply(-y,1,sort, na.last = TRUE))
            y[y == -Inf] = NA 
            as.matrix(y)
        }
        r.vec = as.vector((!is.na(y)) %*% rep(1,ncol(y))) # r_i vector 
        if(any(r.vec==0))
            stop("There is at least one row of the response containing all NAs")
        if(ncol(y) > 1) {
            yiri = y[cbind(1:nrow(y), r.vec)]
            sc.init = if(is.Numeric( .iscale, posit=TRUE))
                        .iscale else {3*(apply(y, 1, mean, na.rm=TRUE) - yiri)}
            sc.init = rep(sc.init, length=nrow(y))
            sc.init[sc.init <= 0.0001] = 1 # Used to be .iscale
            loc.init = yiri + sc.init * log(r.vec)
        } else {
            sc.init =  if(is.Numeric( .iscale, posit=TRUE))
                           .iscale else 1.1 * (0.01+sqrt(var(y)*6)) / pi
            sc.init = rep(sc.init, len=n)
            EulerM = -digamma(1)
            loc.init = (y - sc.init * EulerM)
            loc.init[loc.init <= 0] = min(y)
        }

        extra$R = .R
        extra$mpv = .mpv
        extra$percentiles = .percentiles

        if(!length(etastart)) 
            etastart = cbind(theta2eta(loc.init, .llocation, earg= .elocation),
                             theta2eta(sc.init, .lscale, earg= .escale ))
    }), list( .llocation=llocation, .lscale=lscale, .iscale=iscale,
              .elocation = elocation, .escale = escale,
              .R=R, .mpv=mpv, .percentiles=percentiles ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        loc = eta2theta(eta[,1], .llocation, earg= .elocation)
        sigma = eta2theta(eta[,2], .lscale, earg= .escale )  # sigma
        Percentiles = extra$percentiles
        lp = length(Percentiles)  # may be 0
        if(lp) {
            mpv = extra$mpv
            mu = matrix(as.numeric(NA), nrow(eta), lp + mpv) # lp could be 0
            Rvec = extra$R
            if(1 <= lp)
            for(i in 1:lp) {
                ci = if(is.Numeric(Rvec)) Rvec * (1 - Percentiles[i] / 100) else
                    -log( Percentiles[i] / 100)
                mu[,i] = loc - sigma * log(ci)
            }
            if(mpv) 
                mu[,ncol(mu)] = loc - sigma * log(log(2))
            dmn2 = if(lp>=1) paste(as.character(Percentiles), "%", sep="") else NULL
            if(mpv) 
                dmn2 = c(dmn2, "MPV")
            dimnames(mu) = list(dimnames(eta)[[1]], dmn2)
        } else {
            EulerM = -digamma(1)
            mu = loc + sigma * EulerM
        }
        mu
    }, list( .llocation=llocation, .lscale=lscale,
             .elocation = elocation, .escale = escale ))),
    last=eval(substitute(expression({
        misc$R = .R
        misc$links = c(location = .llocation, scale = .lscale)
        misc$earg= list(location= .elocation, scale= .escale )
        misc$mpv = .mpv
        misc$true.mu = !length( .percentiles) # @fitted is not a true mu
        misc$percentiles = .percentiles
    }), list( .llocation=llocation, .lscale=lscale, .percentiles=percentiles,
              .elocation = elocation, .escale = escale,
              .mpv=mpv, .R=R ))),
    vfamily=c("gumbel", "vextremes"),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals=FALSE, eta, extra=NULL) {
        loc = eta2theta(eta[,1], .llocation, earg= .elocation)
        sigma = eta2theta(eta[,2], .lscale, earg= .escale )
        r.vec = as.vector((!is.na(y)) %*% rep(1,ncol(y))) # r_i vector 
        yiri = y[cbind(1:nrow(y),r.vec)]
        ans = -r.vec * log(sigma) - exp( -(yiri-loc)/sigma )
        max.r.vec = max(r.vec)
        for(jay in 1:max.r.vec) {
            index = (jay <= r.vec)
            ans[index] = ans[index] - (y[index,jay]-loc[index]) / sigma[index]
        }
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * ans)
        }, list( .llocation=llocation, .lscale=lscale,
                 .elocation = elocation, .escale = escale ))),
    deriv=eval(substitute(expression({
        loc = eta2theta(eta[,1], .llocation, earg= .elocation)
        sigma = eta2theta(eta[,2], .lscale, earg= .escale )
        r.vec = as.vector((!is.na(y)) %*% rep(1,ncol(y))) # r_i vector 
        yiri = y[cbind(1:nrow(y),r.vec)]
        yi.bar = apply(y, 1, mean, na.rm=TRUE)
        temp2 = (yiri - loc) / sigma
        term2 = exp(-temp2)
        dloc.deta = dtheta.deta(loc, .llocation, earg= .elocation)
        dsigma.deta = dtheta.deta(sigma, .lscale, earg= .escale )
        dl.dloc = (r.vec - term2) / sigma
        dl.dsigma = (apply((y - loc) / sigma, 1, sum, na.rm=TRUE) - r.vec -
             temp2 * term2) / sigma
        w * cbind(dl.dloc * dloc.deta, dl.dsigma * dsigma.deta)
    }), list( .llocation=llocation, .lscale=lscale,
              .elocation = elocation, .escale = escale ))),
    weight=eval(substitute(expression({
        temp6 = digamma(r.vec)  # , integer=T
        temp5 = digamma(1:max(r.vec))  # , integer=T
        temp5 = matrix(temp5, n, max(r.vec), byrow=TRUE)
        temp5[col(temp5) > r.vec] = 0
        temp5 = temp5 %*% rep(1, ncol(temp5))
        wz = matrix(as.numeric(NA), n, dimm(M=2))  # 3=dimm(M=2)
        wz[,iam(1,1,M)] = r.vec / sigma^2
        wz[,iam(2,1,M)] = -(1 + r.vec * temp6) / sigma^2
        wz[,iam(2,2,M)] = (2*(r.vec+1)*temp6 + r.vec*(trigamma(r.vec) +
                          temp6^2) + 2 - r.vec - 2*temp5) / sigma^2
        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] * dloc.deta^2
        wz[,iam(2,1,M)] = wz[,iam(2,1,M)] * dsigma.deta * dloc.deta
        wz[,iam(2,2,M)] = wz[,iam(2,2,M)] * dsigma.deta^2
        w * wz
    }), list( .lscale=lscale ))))
}






rgpd = function(n, location=0, scale=1, shape=0) {
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    if(!is.Numeric(n, posit=TRUE, allow=1, integ=TRUE)) 
        stop("bad input for argument \"n\"")
    if(!is.Numeric(location)) stop("bad input for argument \"location\"")
    if(!is.Numeric(shape)) stop("bad input for argument \"shape\"")
    ans = numeric(n)
    shape = rep(shape, len=n); location = rep(location, len=n); 
    scale = rep(scale, len=n)
    scase = abs(shape) < sqrt(.Machine$double.eps)
    nscase = sum(scase)
    if(n - nscase)
        ans[!scase] = location[!scase] +
            scale[!scase]*((runif(n-nscase))^(-shape[!scase])-1) / shape[!scase]
    if(nscase)
        ans[scase] = location[scase] - scale[scase] * log(runif(nscase))
    ans
}

dgpd = function(x, location=0, scale=1, shape=0) {
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(location)) stop("bad input for argument \"location\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    if(!is.Numeric(shape)) stop("bad input for argument \"shape\"")
    n = max(length(x), length(location), length(scale), length(shape))
    ans = rep(as.numeric(NA), len=n)
    shape = rep(shape, len=n); location = rep(location, len=n); 
    scale = rep(scale, len=n); x = rep(x-location, len=n)
    scase = abs(shape) < sqrt(.Machine$double.eps)
    nscase = sum(scase)
    if(n - nscase) {
        pos = (1+shape*x/scale > 0 & x > 0)
        ind8 = pos & !scase
        ans[ind8] = (1+shape[ind8]*x[ind8]/scale[ind8])^(-1/shape[ind8]-1) /
                   scale[ind8]
        ind8 = !pos & !scase
        ans[ind8] = 0
    }
    if(nscase) {
        pos = x>0
        ind5 = pos & scase
        ans[ind5] = exp(-x[ind5]/scale[ind5]) / scale[ind5]
        ans[!pos & scase] = 0
    }
    ans
}

pgpd = function(q, location=0, scale=1, shape=0) {
    if(!is.Numeric(q)) stop("bad input for argument \"q\"")
    if(!is.Numeric(location)) stop("bad input for argument \"location\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    if(!is.Numeric(shape)) stop("bad input for argument \"shape\"")
    n = max(length(q), length(location), length(scale), length(shape))
    ans = numeric(n)
    shape = rep(shape, len=n); location = rep(location, len=n); 
    scale = rep(scale, len=n); q = rep(q-location, len=n)
    scase = abs(shape) < sqrt(.Machine$double.eps)
    nscase = sum(scase)
    if(n - nscase) {
        q[q<0] = 0
        ans = 1 - pmax(0,(1+shape*q/scale))^(-1/shape)
    }
    if(nscase) {
        pos = q>=0
        ind9 =  pos & scase
        ans[ind9] =  -expm1(-q[ind9]/scale[ind9])
        ind9 = !pos & scase
        ans[ind9] = 0
    }
    ans
}

qgpd = function(p, location=0, scale=1, shape=0) {
    if(!is.Numeric(p, posit=TRUE) || any(p >= 1)) stop("0 < p < 1 is required")
    if(!is.Numeric(location)) stop("bad input for argument \"location\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    if(!is.Numeric(shape)) stop("bad input for argument \"shape\"")
    n = max(length(p), length(location), length(scale), length(shape))
    ans = numeric(n)
    shape = rep(shape, len=n); location = rep(location, len=n); 
    scale = rep(scale, len=n); p = rep(p, len=n)
    scase = abs(shape) < sqrt(.Machine$double.eps)
    nscase = sum(scase)
    if(n - nscase) {
        ans[!scase] = location[!scase] + scale[!scase] *
            ((1-p[!scase])^(-shape[!scase]) - 1) / shape[!scase]
    }
    if(nscase) {
        ans[scase] = location[scase] - scale[scase] * log1p(-p[scase])
    }
    ans
}





gpd = function(threshold=0,
               lscale="loge",
               lshape="logoff",
               escale = list(),
               eshape = if(lshape=="logoff") list(offset=0.5) else 
               if(lshape=="elogit") list(min=-0.5, max=0.5) else NULL,
               percentiles=c(90,95),
               iscale=NULL,
               ishape=NULL, 
               tshape0=0.001,
               method.init=1,
               zero=2) {
    if(!is.Numeric(threshold)) 
        stop("bad input for argument \"threshold\"")
    if(!is.Numeric(method.init, allow=1, posit=TRUE, integer=TRUE) ||
       method.init > 2.5)
        stop("argument \"method.init\" must be 1 or 2")
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(length(percentiles) && 
      (!is.Numeric(percentiles, posit=TRUE) || max(percentiles) >= 100))
        stop("bad input for argument \"percentiles\"")
    if(!is.Numeric(tshape0, allow=1, posit=TRUE) || tshape0 > 0.1)
        stop("bad input for argument \"tshape0\"")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(escale)) escale = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("Generalized Pareto Distribution\n",
            "Links:    ",
            namesof("scale", link=lscale, earg= escale ), ", ",
            namesof("shape", link=lshape, earg= eshape)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(as.matrix(y)) != 1)
            stop("response must be a vector or one-column matrix")
        Threshold = if(is.Numeric( .threshold)) .threshold else 0
        if(is.Numeric(  .threshold)) {
            orig.y = y
        }
        ystar = y - Threshold  # Operate on ystar
        extra$threshold = Threshold
        predictors.names=
            c(namesof("scale", .lscale, earg= .escale, short=TRUE),
              namesof("shape", .lshape, earg= .eshape, short=TRUE ))
        if(!length(etastart)) {
            meany = mean(ystar)
            vary = var(ystar)
            xiinit = if(length( .ishape)) .ishape else {
                if( .method.init == 1) -0.5*(meany^2/vary - 1) else
                    0.5 * (1 - median(ystar)^2 / vary)
            }
            siginit = if(length( .iscale)) .iscale else {
                if(.method.init==1) 0.5*meany*(meany^2/vary + 1) else
                    abs(1-xiinit) * median(ystar)
            }
            siginit[siginit <= 0] = 0.01    # sigma > 0
            xiinit[xiinit <= -0.5] = -0.40  # Fisher scoring works if xi > -0.5
            xiinit[xiinit >=  1.0] =  0.90  # Mean/var exists if xi < 1 / 0.5
            if( .lshape == "loge") xiinit[xiinit <=  0.0] =  0.05
            siginit = rep(siginit, leng=length(y))
            xiinit = rep(xiinit, leng=length(y))

            etastart = cbind(theta2eta(siginit, .lscale, earg= .escale ),
                             theta2eta(xiinit,  .lshape, earg= .eshape ))
        }
    }), list( .lscale=lscale, .lshape=lshape, .threshold=threshold,
              .iscale=iscale, .ishape=ishape,
              .escale=escale, .eshape=eshape,
              .method.init=method.init ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        sigma = eta2theta(eta[,1], .lscale, earg= .escale )
        xi = eta2theta(eta[,2], .lshape, earg= .eshape )
        cent = .percentiles
        lp = length(cent)  # NULL means lp==0 and the mean is returned
        Threshold = if(is.Numeric( .threshold)) .threshold else 0
        if(lp) {
            fv = matrix(as.numeric(NA), nrow(eta), lp)
            iszero = (abs(xi) < .tshape0)
            for(i in 1:lp) {
                temp = 1-cent[i]/100
                fv[!iszero,i] = Threshold + (temp^(-xi[!iszero]) -1) *
                                sigma[!iszero] / xi[!iszero]
                fv[iszero,i] = Threshold - sigma[iszero] * log(temp)
            }
            dimnames(fv) = list(dimnames(eta)[[1]],
                                paste(as.character(.percentiles), "%", sep=""))
        } else {
            fv = Threshold + sigma / (1 - xi)   # This is the mean, E(Y)
            fv[xi >= 1] = NA  # Mean exists only if xi < 1.
        }
        fv
    }, list( .lscale=lscale, .lshape=lshape, .threshold=threshold,
             .escale=escale, .eshape=eshape,
             .tshape0=tshape0, .percentiles=percentiles ))),
    last=eval(substitute(expression({
        misc$links = c(scale = .lscale, shape = .lshape)
        misc$true.mu = FALSE     # @fitted is not a true mu
        misc$earg= list(scale= .escale , shape= .eshape )
        misc$percentiles = .percentiles
        misc$threshold = if(is.Numeric( .threshold)) .threshold else 0
        misc$expected = TRUE
        misc$tshape0 = .tshape0
        if(any(xi < -0.5))
            warning("some values of the shape parameter are less than -0.5")
    }), list( .lscale=lscale, .lshape=lshape, .threshold=threshold,
              .escale=escale, .eshape=eshape,
              .tshape0=tshape0, .percentiles=percentiles ))),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals=FALSE, eta, extra=NULL) {
            sigma = eta2theta(eta[,1], .lscale, earg= .escale )
            xi = eta2theta(eta[,2], .lshape, earg= .eshape )
            if(any(iszero <- (abs(xi) < .tshape0))) {
            }
            Threshold = extra$threshold
            ystar = y - Threshold  # Operate on ystar
            A = 1 + xi*ystar/sigma
            mytolerance = .Machine$double.eps
            bad <- (A<=mytolerance)   # Range violation
            if(any(sum(w[bad]))) {
                cat("There are some range violations\n")
                if(exists("flush.console")) flush.console()
            }
            igpd = !iszero &  !bad
            iexp =  iszero &  !bad
          if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w[bad]) * (-1.0e10) +
            sum(w[igpd] * (-log(sigma[igpd]) - (1+1/xi[igpd])*log(A[igpd]))) +
            sum(w[iexp] * (-log(sigma[iexp]) - ystar[iexp]/sigma[iexp]))
        }, list( .tshape0=tshape0, .lscale=lscale,
                 .escale=escale, .eshape=eshape,
                 .lshape=lshape ))),
    vfamily=c("gpd", "vextremes"),
    deriv=eval(substitute(expression({
        sigma = eta2theta(eta[,1], .lscale, earg= .escale )
        xi = eta2theta(eta[,2], .lshape, earg= .eshape )
        Threshold = extra$threshold
        ystar = y - Threshold  # Operate on ystar
        A = 1 + xi*ystar/sigma
        mytolerance = .Machine$double.eps
        bad <- (A <= mytolerance)
        if(any(sum(w[bad]))) {
            cat(sum(w[bad],na.rm=TRUE), # "; ignoring them"
                "observations violating boundary constraints\n")
            if(exists("flush.console")) flush.console()
        }
        if(any(iszero <- (abs(xi) < .tshape0))) {
        }
        igpd = !iszero &  !bad
        iexp =  iszero &  !bad
        dl.dxi = dl.dsigma = rep(0, len=length(y))
        dl.dsigma[igpd] = ((1 + xi[igpd]) * ystar[igpd] / (sigma[igpd] + 
                          xi[igpd]*ystar[igpd]) - 1) / sigma[igpd]
        dl.dxi[igpd] = log(A[igpd])/xi[igpd]^2 - (1 + 1/xi[igpd]) *
                       ystar[igpd] / (A[igpd] * sigma[igpd])
        dl.dxi[iexp] = ystar[iexp] *
                       (0.5*ystar[iexp]/sigma[iexp] - 1) / sigma[iexp]
        dsigma.deta = dtheta.deta(sigma, .lscale, earg= .escale )
        dxi.deta = dtheta.deta(xi, .lshape, earg= .eshape )
        w * cbind(dl.dsigma * dsigma.deta, dl.dxi * dxi.deta)
    }), list( .tshape0=tshape0, .lscale=lscale,
              .escale=escale, .eshape=eshape,
              .lshape=lshape ))),
    weight=eval(substitute(expression({
        n <- length(w) # needed! 
        wz = matrix(as.numeric(NA), n, 3)
        wz[,iam(1,1,M)] = 1 / ((1+2*xi) * sigma^2)
        wz[,iam(2,2,M)] = 2 / ((1+2*xi) * (1+xi))
        wz[,iam(1,2,M)] = 1 / ((1+2*xi) * (1+xi) * sigma)  # Positive!!!
        wz[,iam(1,1,M)] = w * wz[,iam(1,1,M)] * dsigma.deta^2
        wz[,iam(2,2,M)] = w * wz[,iam(2,2,M)] * dxi.deta^2
        wz[,iam(1,2,M)] = w * wz[,iam(1,2,M)] * dsigma.deta * dxi.deta
        wz
    }), list( .lscale=lscale ))))
}





meplot.default = function(y, main="Mean Excess Plot",
    xlab="Threshold", ylab="Mean Excess", lty=c(2,1:2), 
    conf=0.95, col=c("blue","black","blue"), type="l", ...) {
    if(!is.Numeric(y)) stop("bad input for argument \"y\"")
    n = length(y)
    sy = sort(y)
    dsy = rev(sy)  # decreasing sequence
    me = rev(cumsum(dsy))/(n:1) - sy
    me2 = rev(cumsum(dsy^2))
    var = (me2 - (n:1)*(me+sy)^2) / (n:1)
    ci = qnorm((1+conf)/2) * sqrt(abs(var)) / sqrt(n:1)
    mymat = cbind(me-ci, me, me+ci)
    sy = sy - sqrt(.Machine$double.eps)
    matplot(sy, mymat, main=main, xlab=xlab, ylab=ylab, 
            lty=lty, col=col, type=type, ...)
    invisible(list(threshold=sy, meanExcess=me))
}

meplot.vlm = function(object, ...) {
    if(!length(y <- object@y)) stop("y slot is empty")
    ans = meplot(as.numeric(y), ...) 
    invisible(ans)
}

if(!isGeneric("meplot"))
    setGeneric("meplot", function(object, ...) standardGeneric("meplot"))

setMethod("meplot", "numeric",
         function(object, ...)
         meplot.default(y=object, ...))

setMethod("meplot", "vlm",
         function(object, ...)
         meplot.vlm(object, ...))



guplot.default = function(y, main="Gumbel Plot",
    xlab="Reduced data", ylab="Observed data", type="p", ...) {
    if(!is.Numeric(y)) stop("bad input for argument \"y\"")
    n = length(y)
    sy = sort(y)
    x = -log(-log(((1:n) - 0.5) / n))
    plot(x, sy, main=main, xlab=xlab, ylab=ylab, type=type, ...)
    invisible(list(x=x, y=sy))
}

guplot.vlm = function(object, ...) {
    if(!length(y <- object@y)) stop("y slot is empty")
    ans = guplot(as.numeric(y), ...) 
    invisible(ans)
}

if(!isGeneric("guplot"))
    setGeneric("guplot", function(object, ...) standardGeneric("guplot"))

setMethod("guplot", "numeric",
         function(object, ...)
         guplot.default(y=object, ...))

setMethod("guplot", "vlm",
         function(object, ...)
         guplot.vlm(object, ...))



if(FALSE)
bvevd.log.qn <- function(lscale="loge",
                         ldep="logit",
                         exchangeable=FALSE, zero=NULL)
{
    if(mode(ldep) != "character" && mode(ldep) != "name")
        ldep <- as.character(substitute(ldep))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale <- as.character(substitute(lscale))

    new("vglmff",
    blurb=c("Bivariate logistic extreme value distribution\n\n",
           "Link:     ", "loc1, ", namesof("scale1", lscale), ", shape1, ",
                         "loc2, ", namesof("scale2", lscale), ", shape2, ", 
                         namesof("dep", ldep), "\n"),
    constraints=eval(substitute(expression({
        constraints <- cm.vgam(cbind(rbind(diag(3),diag(3),matrix(0,1,3)),
                                     ei(7,7)), x, .exchangeable, constraints)
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .exchangeable=exchangeable, .zero=zero ))),
    deviance= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        nz <- y > 0
        devi <-  - (y - mu)
        devi[nz] <- devi[nz] + y[nz] * log(y[nz]/mu[nz])
        if(residuals) sign(y - mu) * sqrt(2 * abs(devi) * w) else
            2 * sum(w * devi)
    },
    initialize=eval(substitute(expression({
        M = if(is.matrix(y)) ncol(y) else 1
        if(M != 1) stop("response must be a vector or a one-column matrix")
        predictors.names =
        c("loc1",
          namesof("scale1", lscale), "shape1",
          "loc2",
          namesof("scale2", lscale), "shape2",
          namesof("dep", ldep))
        if(!length(etastart)) {
            etastart <- theta2eta(mu, link=.link)
        }
    }), list( .link=link, .estimated.dispersion=estimated.dispersion) )),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mu = eta2theta(eta, link=.link)
        mu
    }, list( .link=link) )),
    last=eval(substitute(expression({
        misc$link = c(loc1="identity", scale1= .lscale, shape1="identity",
                      loc2="identity", scale2= .lscale, shape2="identity",
                      dep= .ldep)
        misc$expected = FALSE
        misc$BFGS = TRUE
    }), list( .lscale=lscale, .ldep=ldep) )),
    link=eval(substitute(function(mu, extra=NULL) {
        theta2eta(mu, link=.link)
    }, list( .link=link) )),
    loglikelihood= function(mu, y, w, residuals = FALSE, eta, extra=NULL) {
        pdf = deriv3( ~ w * ifelse(1+xi1*(y1-mu1)/sigma1>0, y1, 0),
                        name=c("y1"), hessian= TRUE)
        pdf = deriv3( ~ w * exp(-((pmax(0,1+xi1*(y1-mu1)/sigma1)^(-1/xi1))^(1/dep) +
                   (pmax(0,1+xi2*(y2-mu2)/sigma2)^(-1/xi2))^(1/dep))^(dep)),
                        name=c("y1","y2"), hessian= TRUE)
          if(residuals) stop("loglikelihood residuals not implemented yet") else
        pdf
 
    },
    vfamily="bvevd.log.qn",
    deriv=eval(substitute(expression({
        if(iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }

        derivnew =
        if(.link == "loge" && (any(mu < .Machine$double.eps))) {
            w * (y - mu)
        } else {
            lambda <- mu
            dl.dlambda <- (y-lambda) / lambda
            dlambda.deta <- dtheta.deta(theta=lambda, link=.link)
            w * dl.dlambda * dlambda.deta
        }
        derivnew
    }), list( .link=link) )),
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
        }
        wznew
    }), list( .lscale=lscale, .ldep=ldep) )))
}




rgev = function(n, location=0, scale=1, shape=0) {
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    if(!is.Numeric(n, posit=TRUE, allow=1, integ=TRUE)) 
        stop("bad input for argument \"n\"")
    if(!is.Numeric(location)) 
        stop("bad input for argument argument \"location\"")
    if(!is.Numeric(shape)) stop("bad input for argument argument \"shape\"")
    ans = numeric(n)
    shape = rep(shape, len=n); location = rep(location, len=n); 
    scale = rep(scale, len=n)
    scase = abs(shape) < sqrt(.Machine$double.eps)
    nscase = sum(scase)
    if(n - nscase)
        ans[!scase] = location[!scase] + scale[!scase] *
            ((-log(runif(n-nscase)))^(-shape[!scase]) -1) / shape[!scase]
    if(nscase)
        ans[scase] = rgumbel(nscase, location[scase], scale[scase])
    ans
}

dgev = function(x, location=0, scale=1, shape=0) {
    if(!is.Numeric(x)) stop("bad input for argument \"x\"")
    if(!is.Numeric(location)) stop("bad input for argument \"location\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    if(!is.Numeric(shape)) stop("bad input for argument \"shape\"")
    n = max(length(x), length(location), length(scale), length(shape))
    ans = rep(as.numeric(NA), len=n)
    shape = rep(shape, len=n); location = rep(location, len=n); 
    scale = rep(scale, len=n); x = rep(x-location, len=n)
    scase = abs(shape) < sqrt(.Machine$double.eps)
    nscase = sum(scase)
    if(n - nscase) {
        zedd = pmax(0,(1+shape*x/scale))
        ans[!scase] = exp(-zedd[!scase]^(-1/shape[!scase])) *
              zedd[!scase]^(-1/shape[!scase] -1) / scale[!scase]
    }
    if(nscase)
        ans[scase] = dgumbel(x[scase], location[scase], scale[scase])
    ans
}

pgev = function(q, location=0, scale=1, shape=0) {
    if(!is.Numeric(q)) stop("bad input for argument \"q\"")
    if(!is.Numeric(location)) stop("bad input for argument \"location\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    if(!is.Numeric(shape)) stop("bad input for argument \"shape\"")
    n = max(length(q), length(location), length(scale), length(shape))
    ans = numeric(n)
    shape = rep(shape, len=n); location = rep(location, len=n); 
    scale = rep(scale, len=n); q = rep(q-location, len=n)
    scase = abs(shape) < sqrt(.Machine$double.eps)
    nscase = sum(scase)
    if(n - nscase) {
        zedd = pmax(0,(1+shape*q/scale))
        ans[!scase] = exp(-zedd[!scase]^(-1/shape[!scase]))
    }
    if(nscase)
        ans[scase] = pgumbel(q[scase], location[scase], scale[scase])
    ans
}

qgev = function(p, location=0, scale=1, shape=0) {
    if(!is.Numeric(p, posit=TRUE) || any(p >= 1)) stop("0 < p < 1 is required")
    if(!is.Numeric(location)) stop("bad input for argument \"location\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("\"scale\" must be positive")
    if(!is.Numeric(shape)) stop("bad input for argument \"shape\"")
    n = max(length(p), length(location), length(scale), length(shape))
    ans = numeric(n)
    shape = rep(shape, len=n); location = rep(location, len=n); 
    scale = rep(scale, len=n); p = rep(p, len=n)
    scase = abs(shape) < sqrt(.Machine$double.eps)
    nscase = sum(scase)
    if(n - nscase) {
        ans[!scase] = location[!scase] + scale[!scase] *
            ((-log(p[!scase]))^(-shape[!scase]) -1) / shape[!scase]
    }
    if(nscase)
        ans[scase] = qgumbel(p[scase], location[scase], scale[scase])
    ans
}



egumbel = function(llocation="identity",
                   lscale="loge",
                   elocation = list(),
                   escale = list(),
                   iscale=NULL,
                   R=NA, percentiles=c(95,99),
                   mpv=FALSE, zero=NULL)
{
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(!is.logical(mpv) || length(mpv) != 1)
        stop("bad input for argument \"mpv\"")
    if(length(percentiles) &&
       (!is.Numeric(percentiles, posit=TRUE) || max(percentiles) >= 100))
        stop("bad input for argument \"percentiles\"")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(length(iscale) && !is.Numeric(iscale, posit=TRUE))
        stop("bad input for argument \"iscale\"")
    if(!is.list(elocation)) elocation = list()
    if(!is.list(escale)) escale = list()

    new("vglmff",
    blurb=c("Gumbel distribution (univariate response)\n\n",
            "Links:    ",
            namesof("location", llocation, earg= elocation, tag= TRUE), ", ", 
            namesof("scale", lscale, earg= escale , tag= TRUE), "\n",
            "Mean:     location + scale*0.5772..\n",
            "Variance: pi^2 * scale^2 / 6"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        y = cbind(y)
        if(ncol(y) > 1)
            stop("Use gumbel() to handle multivariate responses")
        if(min(y) <= 0)
            stop("all response values must be positive")
        predictors.names =
        c(namesof("location", .llocation, earg= .elocation, tag= FALSE),
          namesof("scale", .lscale, earg= .escale , tag= FALSE))

        extra$R = .R
        extra$mpv = .mpv
        extra$percentiles = .percentiles

        if(!length(etastart)) {
            sc.init =  if(is.Numeric( .iscale, posit=TRUE)) 
                           .iscale else 1.5 * (0.01+sqrt(var(y)*6)) / pi
            sc.init = rep(sc.init, len=n)
            EulerM = -digamma(1)
            loc.init = (y - sc.init * EulerM)
            etastart = cbind(theta2eta(loc.init, .llocation, earg= .elocation),
                             theta2eta(sc.init,  .lscale, earg= .escale ))
        }
    }), list( .llocation=llocation, .lscale=lscale, .iscale=iscale, 
              .elocation=elocation, .escale=escale,
              .R=R, .mpv=mpv, .percentiles=percentiles ))),
    inverse=eval(substitute( function(eta, extra=NULL) {
        loc = eta2theta(eta[,1], .llocation, earg= .elocation)
        sigma = eta2theta(eta[,2], .lscale, earg= .escale )
        EulerM = -digamma(1)
        Percentiles = extra$percentiles
        mpv = extra$mpv
        lp = length(Percentiles)  # may be 0
        if(!lp) return(loc + sigma * EulerM)
        mu = matrix(as.numeric(NA), nrow(eta), lp + mpv)
        Rvec = extra$R
        if(1 <= lp)
        for(i in 1:lp) {
            ci = if(is.Numeric(Rvec)) Rvec * (1 - Percentiles[i] / 100) else
                -log( Percentiles[i] / 100)
            mu[,i] = loc - sigma * log(ci)
        }
        if(mpv)
            mu[,ncol(mu)] = loc - sigma * log(log(2))
        dmn2 = if(lp>=1) paste(as.character(Percentiles), "%", sep="") else NULL
        if(mpv)
            dmn2 = c(dmn2, "MPV")
        dimnames(mu) = list(dimnames(eta)[[1]], dmn2)
        mu
    }, list( .llocation=llocation, .lscale=lscale,
             .elocation=elocation, .escale=escale ))),
    last=eval(substitute(expression({
        misc$link = c(location= .llocation, scale = .lscale) 
        misc$earg= list(location= .elocation, scale= .escale)
        misc$true.mu = !length( .percentiles) # @fitted is not a true mu
        misc$R = .R
        misc$mpv = .mpv
        misc$percentiles = .percentiles
    }), list( .llocation=llocation, .lscale=lscale, .mpv=mpv,
              .elocation=elocation, .escale=escale,
              .R=R, .percentiles=percentiles ))),
    loglikelihood=eval(substitute(
            function(mu,y,w,residuals= FALSE,eta,extra=NULL) {
        loc = eta2theta(eta[,1], .llocation, earg= .elocation)
        sc  = eta2theta(eta[,2], .lscale, earg= .escale )
        zedd = (y-loc) / sc
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-log(sc) - zedd - exp(-zedd)))
    }, list( .llocation=llocation, .lscale=lscale,
             .elocation=elocation, .escale=escale ))),
    vfamily="egumbel",
    deriv=eval(substitute(expression({
        loc = eta2theta(eta[,1], .llocation, earg= .elocation)
        sc = eta2theta(eta[,2], .lscale, earg= .escale )
        zedd = (y-loc) / sc
        temp2 = -expm1(-zedd)
        dl.dloc = temp2 / sc
        dl.dsc = -1/sc + temp2 * zedd / sc
        dloc.deta = dtheta.deta(loc, .llocation, earg= .elocation)
        dsc.deta = dtheta.deta(sc, .lscale, earg= .escale )
        w * cbind(dl.dloc * dloc.deta, dl.dsc * dsc.deta)
    }), list( .llocation=llocation, .lscale=lscale,
              .elocation=elocation, .escale=escale ))),
    weight=expression({
        digamma1 = digamma(1)
        ed2l.dsc2 = ((2+digamma1)*digamma1 + trigamma(1) + 1) / sc^2
        ed2l.dloc2 = 1 / sc^2
        ed2l.dscloc = -(1 + digamma1) / sc^2 
        wz = matrix(as.numeric(NA), n, dimm(M=2))
        wz[,iam(1,1,M)] = ed2l.dloc2 * dloc.deta^2
        wz[,iam(2,2,M)] = ed2l.dsc2 * dsc.deta^2
        wz[,iam(1,2,M)] = ed2l.dscloc * dloc.deta * dsc.deta
        w * wz
    }))
}




cgumbel = function(llocation="identity",
                   lscale="loge",
                   elocation = list(),
                   escale = list(), iscale=NULL,
                   mean=TRUE, percentiles=NULL, zero=2)
{
    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(!is.logical(mean) || length(mean) != 1)
        stop("mean must be a single logical value")
    if(!mean && (!is.Numeric(percentiles, posit=TRUE) ||
                 any(percentiles>=100)))
        stop("valid percentiles values must be given when mean=FALSE")
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.list(elocation)) elocation = list()
    if(!is.list(escale)) escale = list()

    new("vglmff",
    blurb=c("Censored Gumbel distribution\n\n",
            "Links:    ",
            namesof("location", llocation, earg= elocation, tag= TRUE), ", ", 
            namesof("scale", lscale, earg= escale, tag= TRUE),
            "\n",
            "Mean:     location + scale*0.5772..\n",
            "Variance: pi^2 * scale^2 / 6"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        y = cbind(y)
        if(ncol(y) > 1)
            stop("Use gumbel.block() to handle multivariate responses")
        if(any(y) <= 0)
            stop("all response values must be positive")

        if(!length(extra$leftcensored)) extra$leftcensored = rep(FALSE, len=n)
        if(!length(extra$rightcensored)) extra$rightcensored = rep(FALSE, len=n)
        if(any(extra$rightcensored & extra$leftcensored))
            stop("some observations are both right and left censored!")

        predictors.names =
        c(namesof("location", .llocation, earg= .elocation, tag= FALSE),
          namesof("scale", .lscale, earg= .escale , tag= FALSE))
        if(!length(etastart)) {
            sc.init =  if(is.Numeric( .iscale, posit=TRUE)) 
                           .iscale else 1.1 * sqrt(var(y) * 6 ) / pi
            sc.init = rep(sc.init, len=n)
            EulerM = -digamma(1)
            loc.init = (y - sc.init * EulerM)
            loc.init[loc.init <= 0] = min(y)
            etastart = cbind(theta2eta(loc.init, .llocation, earg= .elocation ),
                             theta2eta(sc.init,  .lscale, earg= .escale ))
        }
    }), list( .lscale=lscale, .iscale=iscale,
              .llocation = llocation,
              .elocation = elocation, .escale = escale ))), 
    inverse=eval(substitute( function(eta, extra=NULL) {
        loc  = eta2theta(eta[,1], .llocation)
        sc   = eta2theta(eta[,2], .lscale)
        EulerM = -digamma(1)
        if(.mean) loc + sc * EulerM else {
            lp = length(.percentiles)  # 0 if NULL
            mu = matrix(as.numeric(NA), nrow(eta), lp)
            for(i in 1:lp) {
                ci = -log( .percentiles[i] / 100)
                mu[,i] = loc - sc * log(ci)
            }
            dmn2 = paste(as.character(.percentiles), "%", sep="")
            dimnames(mu) <- list(dimnames(eta)[[1]], dmn2)
            mu
        }
    }, list( .lscale=lscale, .percentiles=percentiles,
             .llocation = llocation,
             .elocation = elocation, .escale = escale ,
             .mean=mean ))), 
    last=eval(substitute(expression({
        misc$link = c(location= .llocation, scale = .lscale) 
        misc$earg= list(location= .elocation, scale= .escale )
        misc$true.mu = .mean    # if FALSE then @fitted is not a true mu 
        misc$percentiles = .percentiles
    }), list( .lscale=lscale, .mean=mean,
              .llocation = llocation,
              .elocation = elocation, .escale = escale ,
              .percentiles=percentiles ))),
    loglikelihood=eval(substitute(
            function(mu,y,w,residuals= FALSE,eta,extra=NULL) {
        loc = eta2theta(eta[,1], .llocation, earg= .elocation )
        sc  = eta2theta(eta[,2], .lscale, earg= .escale )
        zedd = (y-loc) / sc

        cenL = extra$leftcensored
        cenU = extra$rightcensored
        cen0 = !cenL & !cenU   # uncensored obsns
        Fy = exp(-exp(-zedd))
        ell1 = -log(sc[cen0]) - zedd[cen0] - exp(-zedd[cen0])
        ell2 = log(Fy[cenL])
        ell3 = log1p(-Fy[cenU])
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w[cen0] * ell1) + sum(w[cenL] * ell2) + sum(w[cenU] * ell3)
    }, list( .lscale=lscale,
             .llocation = llocation,
             .elocation = elocation, .escale = escale ))),
    vfamily="cgumbel",
    deriv=eval(substitute(expression({
        cenL = extra$leftcensored
        cenU = extra$rightcensored
        cen0 = !cenL & !cenU   # uncensored obsns

        loc = eta2theta(eta[,1], .llocation, earg= .elocation )
        sc  = eta2theta(eta[,2], .lscale, earg= .escale )
        zedd = (y-loc) / sc
        temp2 = -expm1(-zedd)
        dl.dloc = temp2 / sc
        dl.dsc = -1/sc + temp2 * zedd / sc
        dloc.deta = dtheta.deta(loc, .llocation, earg= .elocation )
        dsc.deta = dtheta.deta(sc, .lscale, earg= .escale )

        ezedd = exp(-zedd)
        Fy = exp(-ezedd)
        dFy.dloc = -ezedd * Fy / sc
        dFy.dsc = zedd * dFy.dloc # -zedd * exp(-zedd) * Fy / sc
        if(any(cenL)) {
            dl.dloc[cenL] = -ezedd[cenL] / sc[cenL]
            dl.dsc[cenL] = -zedd[cenL] * ezedd[cenL] / sc[cenL]
        }
        if(any(cenU)) {
            dl.dloc[cenU] = -dFy.dloc[cenU] / (1-Fy[cenU])
            dl.dsc[cenU] = -dFy.dsc[cenU] / (1-Fy[cenU])
        }
        w * cbind(dl.dloc * dloc.deta, dl.dsc * dsc.deta)
    }), list( .lscale=lscale,
              .llocation = llocation,
              .elocation = elocation, .escale = escale ))),
    weight=expression({
        A1 = ifelse(cenL, Fy, 0)
        A3 = ifelse(cenU, 1-Fy, 0)
        A2 = 1 - A1 - A3   # Middle; uncensored
        digamma1 = digamma(1)
        ed2l.dsc2 = ((2+digamma1)*digamma1 + trigamma(1) + 1) / sc^2
        ed2l.dloc2 = 1 / sc^2
        ed2l.dlocsc = -(1 + digamma1) / sc^2 
        wz = matrix(as.numeric(NA), n, dimm(M=2))
        wz[,iam(1,1,M)] = A2 * ed2l.dloc2 * dloc.deta^2
        wz[,iam(2,2,M)] = A2 * ed2l.dsc2 * dsc.deta^2
        wz[,iam(1,2,M)] = A2 * ed2l.dlocsc * dloc.deta * dsc.deta
        d2l.dloc2 = -ezedd / sc^2
        d2l.dsc2 = (2 - zedd) * zedd * ezedd / sc^2
        d2l.dlocsc = (1 - zedd) * ezedd / sc^2
        wz[,iam(1,1,M)]=wz[,iam(1,1,M)]-A1^2 * d2l.dloc2 * dloc.deta^2
        wz[,iam(2,2,M)]=wz[,iam(2,2,M)]-A1^2 * d2l.dsc2 * dsc.deta^2
        wz[,iam(1,2,M)]=wz[,iam(1,2,M)]-A1^2 * d2l.dlocsc * dloc.deta * dsc.deta
        d2Fy.dloc2 = dFy.dloc * dl.dloc + Fy * d2l.dloc2
        d2Fy.dsc2 = dFy.dsc * dl.dsc + Fy * d2l.dsc2
        d2Fy.dlocsc = dFy.dsc * dl.dloc + Fy * d2l.dlocsc
        d2l.dloc2 = -((1-Fy) * d2Fy.dloc2 - dFy.dloc^2) / (1-Fy)^2
        d2l.dsc2 = -((1-Fy) * d2Fy.dsc2 - dFy.dsc^2) / (1-Fy)^2
        d2l.dlocsc =-((1-Fy) * d2Fy.dlocsc - dFy.dloc * dFy.dsc) / (1-Fy)^2
        wz[,iam(1,1,M)]=wz[,iam(1,1,M)]-A3^2 * d2l.dloc2 * dloc.deta^2
        wz[,iam(2,2,M)]=wz[,iam(2,2,M)]-A3^2 * d2l.dsc2 * dsc.deta^2
        wz[,iam(1,2,M)]=wz[,iam(1,2,M)]-A3^2 * d2l.dlocsc * dloc.deta * dsc.deta
        w * wz
    }))
}



dfrechet = function(x, location=0, scale=1, shape) {
    if(!is.Numeric(scale, posit=TRUE)) stop("scale must be positive")
    if(!is.Numeric(shape, posit=TRUE)) stop("shape must be positive")
    rzedd = scale / (x - location)
    ans = rzedd^2 * shape * exp(-(rzedd^shape)) * rzedd^(shape-1) / scale
    ans[x <= location] = 0
    ans
}

pfrechet = function(q, location=0, scale=1, shape) {
    if(!is.Numeric(scale, posit=TRUE)) stop("scale must be positive")
    if(!is.Numeric(shape, posit=TRUE)) stop("shape must be positive")
    rzedd = scale / (q - location)
    ans = exp(-(rzedd^shape))
    ans[q <= location] = 0
    ans
}

qfrechet = function(p, location=0, scale=1, shape) {
    if(!is.Numeric(p, posit=TRUE) || any(p >= 1)) stop("0 < p < 1 is required")
    if(!is.Numeric(scale, posit=TRUE)) stop("scale must be positive")
    if(!is.Numeric(shape, posit=TRUE)) stop("shape must be positive")
    location + scale * (-log(p))^(-1/shape)
}

rfrechet = function(n, location=0, scale=1, shape) {
    if(!is.Numeric(n, posit=TRUE, allow=1, integ=TRUE)) 
        stop("bad input for argument \"n\"")
    if(!is.Numeric(scale, posit=TRUE)) stop("scale must be positive")
    if(!is.Numeric(shape, posit=TRUE)) stop("shape must be positive")
    location + scale * (-log(runif(n)))^(-1/shape)
}

frechet2.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}

frechet2 = function(location=0,
                    lscale="loge",
                    lshape="loglog",
                escale = list(),
                eshape = list(),
                    iscale=NULL, ishape=3,
                    zero=NULL)
{
    if(!is.Numeric(location))
        stop("bad input for argument \"location\"")
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale <- as.character(substitute(lscale))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape <- as.character(substitute(lshape))
    if(!is.list(escale)) escale = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("2-parameter Frechet Distribution\n",
            "Links:    ",
            namesof("scale", link=lscale, earg=escale ), ", ",
            namesof("shape", link=lshape, earg=eshape )),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names =
        c(namesof("scale", .lscale, earg=.escale, short=TRUE),
          namesof("shape", .lshape, earg=.eshape, short=TRUE))
        extra$location = rep( .location, len=n) # stored here
        if(!length(etastart)) {
            # Initial values for limiting case as xi --> 0, r_i==1
            locinit = extra$location
            if(any(y <= locinit))
                stop("initial values for location are out of range")
            shape.init = if(length( .ishape)) rep( .ishape, len=n) else {
                rep(3.0, len=n)   # variance exists if shape>2
            }
            Scale.init = if(length( .iscale)) rep( .iscale, len=n) else {
                if(all(shape.init > 1))
                abs( (y-locinit+0.001) / (gamma(1-1/shape.init)) ) else
                     rep( 1.0, len=n)
            }
            etastart = cbind(theta2eta(Scale.init, .lscale, earg=.escale ), 
                             theta2eta(shape.init, .lshape, earg=.escale ))
        }
    }), list( .lscale=lscale, .lshape=lshape,
              .escale = escale, .eshape= eshape,
              .location=location, .iscale=iscale, .ishape=ishape ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        loc = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale )
        shape = eta2theta(eta[,2], .lshape, earg= .eshape )
        ans = rep(as.numeric(NA), len=length(shape))
        ok = shape > 1
        ans[ok] = loc[ok] + Scale[ok] * gamma(1 - 1/shape[ok])
        ans
    }, list( .lscale=lscale, .lshape=lshape,
             .escale=escale, .eshape=eshape ))),
    last=eval(substitute(expression({
        misc$links <- c("scale"= .lscale, "shape"= .lshape)
        misc$expected = FALSE
        misc$BFGS = TRUE
    }), list( .lscale=lscale, .lshape=lshape,
              .escale=escale, .eshape=eshape ))),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals=FALSE, eta, extra=NULL) {
        loc = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale )
        shape = eta2theta(eta[,2], .lshape, earg= .eshape )
        rzedd = Scale / (y-loc)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (log(shape) + log(Scale) - 2 * log(y-loc) -
                     rzedd^shape + (shape-1) * log(rzedd)))
    }, list( .lscale=lscale, .lshape=lshape,
             .escale=escale, .eshape=eshape ))),
    vfamily=c("frechet2", "vextremes"),
    deriv=eval(substitute(expression({
        loc = extra$location
        Scale = eta2theta(eta[,1], .lscale, earg= .escale )
        shape = eta2theta(eta[,2], .lshape, earg= .eshape )
        rzedd = Scale / (y-loc)   # reciprocial of zedd
        dl.dloc = (shape+1)/(y-loc) - (shape / (y-loc)) * (rzedd)^shape
        dl.dScale = shape * (1-rzedd^shape) / Scale
        dl.dshape = 1/shape + log(rzedd) * (1 -  rzedd^shape)
        if(iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }
        derivnew = w *
        cbind(dl.dScale * dtheta.deta(Scale, .lscale, earg= .escale ),
              dl.dshape * dtheta.deta(shape, .lshape, earg= .eshape ))
        derivnew
    }), list( .lscale=lscale, .lshape=lshape,
              .escale=escale, .eshape=eshape ))),
    weight=eval(substitute(expression({
        if(iter == 1) {
            wznew = cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
        } else {
            wzold = wznew
            wznew = qnupdate(w=w, wzold=wzold, dderiv=(derivold - derivnew),
                             deta=etanew-etaold, M=M,
                             trace=trace)  # weights incorporated in args
        }
        wznew
    }), list( .lscale=lscale, .lshape=lshape ))))
}



frechet3.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}


frechet3= function(anchor=NULL,
                   ldifference="loge",
                   lscale="loge",
                   lshape="loglog",
                   edifference = list(),
                   escale = list(),
                   eshape = list(),
                   ilocation=NULL, iscale=NULL, ishape=3, zero=NULL,
                   effpos = .Machine$double.eps^0.75)
{
    if(mode(ldifference) != "character" && mode(ldifference) != "name")
        ldifference <- as.character(substitute(ldifference))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale <- as.character(substitute(lscale))
    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape <- as.character(substitute(lshape))
    if(!is.Numeric(ishape, allo=1, posi=TRUE)) stop("bad input for argument \"ishape\"")
    if(!is.Numeric(effpos, allo=1)|| effpos<0) stop("bad input for argument \"effpos\"")
    if(!is.list(edifference)) edifference = list()
    if(!is.list(escale)) escale = list()
    if(!is.list(eshape)) eshape = list()

    new("vglmff",
    blurb=c("3-parameter Frechet Distribution\n",
            "Links:    ",
            namesof("difference", link=ldifference, earg=edifference), ", ", 
            namesof("scale", link=lscale, earg=escale), ", ",
            namesof("shape", link=lshape, earg=eshape)),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names =
        c(namesof("difference",  .ldifference, earg= .edifference, short=TRUE),
          namesof("scale",  .lscale, earg= .escale,  short=TRUE),
          namesof("shape",  .lshape, earg= .eshape,  short=TRUE))
        anchorpt = if(is.Numeric( .anchor, allow=1)) .anchor else min(y)
        if(min(y) < anchorpt) stop("anchor point is too large")
        extra$LHSanchor = anchorpt
        if(!length(etastart)) {
            # Initial values for limiting case as xi --> 0, r_i==1
            locinit = if(length( .ilocation)) rep( .ilocation, len=n) else
                      rep(anchorpt - 0.01*diff(range(y)), len=n)
            if(any(y <= locinit))
                stop("initial values for location are out of range")
            if(any(anchorpt <= locinit))
                stop("require anchor point > initial location parameter value")
            shape.init = if(length( .ishape)) rep( .ishape, len=n) else {
                rep(3.0, len=n)   # variance exists if shape>2
            }
            Scale.init = if(length( .iscale)) rep( .iscale, len=n) else {
                if(all(shape.init > 1))
                abs( (y-locinit+0.001) / (gamma(1-1/shape.init)) ) else
                     rep( 1.0, len=n)
            }
            etastart = cbind(theta2eta(anchorpt - locinit, .ldifference),
                             theta2eta(Scale.init, .lscale), 
                             theta2eta(shape.init, .lshape))
        }
    }), list( .ldifference=ldifference, .lscale=lscale, .lshape=lshape, 
              .edifference=edifference, .escale=escale, .eshape=eshape, 
              .anchor=anchor,
              .ilocation=ilocation, .iscale=iscale, .ishape=ishape ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        loc = extra$LHSanchor - eta2theta(eta[,1], .ldifference, earg= .edifference)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale )
        shape = eta2theta(eta[,3], .lshape, earg= .eshape )
        ans = rep(as.numeric(NA), len=length(shape))
        ok = shape > 1
        ans[ok] = loc[ok] + Scale[ok] * gamma(1 - 1/shape[ok])
        ans
    }, list( .ldifference=ldifference, .lscale=lscale, .lshape=lshape,
             .edifference=edifference, .escale=escale, .eshape=eshape ))), 
    last=eval(substitute(expression({
        misc$links <- c("difference"= .ldifference, "scale"= .lscale,
                        "shape"= .lshape)
        misc$expected = FALSE
        misc$BFGS = TRUE
    }), list( .ldifference=ldifference, .lscale=lscale, .lshape=lshape,
              .edifference=edifference, .escale=escale, .eshape=eshape ))),  
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals=FALSE, eta, extra=NULL) {
        loc = extra$LHSanchor - eta2theta(eta[,1], .ldifference, earg= .edifference)
        Scale = eta2theta(eta[,2], .lscale, earg= .escale )
        shape = eta2theta(eta[,3], .lshape, earg= .eshape )
        rzedd = Scale / (y-loc)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (log(shape) + log(Scale) - 2 * log(y-loc) -
                     rzedd^shape + (shape-1) * log(rzedd)))
        }, list( .ldifference=ldifference, .lscale=lscale, .lshape=lshape,
                 .edifference=edifference, .escale=escale, .eshape=eshape ))),
    vfamily=c("frechet3", "vextremes"),
    deriv=eval(substitute(expression({
        difference = eta2theta(eta[,1], .ldifference, earg= .edifference )
        Scale = eta2theta(eta[,2], .lscale, earg= .escale )
        shape = eta2theta(eta[,3], .lshape, earg= .eshape )
        loc = extra$LHSanchor - difference
        extra$location = loc   # Store the location parameter estimate here
        rzedd = Scale / (y-loc)   # reciprocial of zedd
        dl.dloc = (shape+1)/(y-loc) - (shape / (y-loc)) * (rzedd)^shape
        dl.ddiff = -dl.dloc
        dl.dScale = shape * (1-rzedd^shape) / Scale
        dl.dshape = 1/shape + log(rzedd) * (1 -  rzedd^shape)
        if(iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }
        derivnew = w *
        cbind(dl.ddiff  * dtheta.deta(difference, .ldifference, earg= .edifference ),
              dl.dScale * dtheta.deta(Scale, .lscale, earg= .escale ),
              dl.dshape * dtheta.deta(shape, .lshape, earg= .eshape ))
        derivnew
    }), list( .ldifference=ldifference, .lscale=lscale, .lshape=lshape,
              .edifference=edifference, .escale=escale, .eshape=eshape ))),
    weight=eval(substitute(expression({
        if(iter == 1) {
            wznew = cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
        } else {
            wzold = wznew
            wznew = qnupdate(w=w, wzold=wzold, dderiv=(derivold - derivnew),
                             deta=etanew-etaold, M=M, effpos = .effpos,
                             trace=trace)  # weights incorporated in args
        }
        wznew
    }), list( .effpos=effpos ))))
}


recnormal1.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}

recnormal1 = function(lmean="identity", lsd="loge",
                      imean=NULL, isd=NULL, method.init=1, zero=NULL)
{

    if(mode(lmean) != "character" && mode(lmean) != "name")
        lmean = as.character(substitute(lmean))
    if(mode(lsd) != "character" && mode(lsd) != "name")
        lsd = as.character(substitute(lsd))
    if(!is.Numeric(method.init, allow=1, integ=TRUE, positi=TRUE) ||
       method.init > 3.5)
        stop("argument \"method.init\" must be 1 or 2 or 3")

    new("vglmff",
    blurb=c("Upper record values from a univariate normal distribution\n\n",
            "Links:    ",
            namesof("mean", lmean, tag= TRUE), "; ",
            namesof("sd", lsd, tag= TRUE),
            "\n",
            "Variance: sd^2"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names = c(namesof("mean", .lmean, tag= FALSE),
                             namesof("sd",   .lsd, tag= FALSE))
        if(ncol(y <- cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(any(diff(y) <= 0))
            stop("response must have increasingly larger and larger values")
        if(any(w != 1))
            warning("weights should have unit values only")
        if(!length(etastart)) {
            mean.init = if(length( .imean)) rep( .imean, len=n) else {
                if(.lmean == "loge") pmax(1/1024, min(y)) else min(y)}
            sd.init = if(length( .isd)) rep( .isd, len=n) else {
                if(.method.init == 1)  1*(sd(y)) else
                if(.method.init == 2)  5*(sd(y)) else
                                      .5*(sd(y))
                }
            etastart = cbind(theta2eta(rep(mean.init, len=n), .lmean),
                             theta2eta(rep(sd.init, len=n), .lsd))
        }
    }), list( .lmean=lmean, .lsd=lsd, .imean=imean, .isd=isd,
             .method.init=method.init ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta[,1], .lmean)
    }, list( .lmean=lmean ))),
    last=eval(substitute(expression({
        misc$link = c("mu"= .lmean, "sd"= .lsd)
        misc$expected = FALSE
    }), list( .lmean=lmean, .lsd=lsd ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        sd = eta2theta(eta[,2], .lsd)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            zedd = (y - mu) / sd
            NN = nrow(eta)
            sum(w * (-log(sd) - 0.5 * zedd^2)) -
            sum(w[-NN] * (log1p(-pnorm(zedd[-NN]))))
        }
    }, list( .lsd=lsd ))),
    vfamily=c("recnormal1"),
    deriv=eval(substitute(expression({
        NN = nrow(eta)
        mymu = eta2theta(eta[,1], .lmean)
        sd = eta2theta(eta[,2], .lsd)
        zedd = (y - mymu) / sd
        temp200 = dnorm(zedd) / (1-pnorm(zedd))
        dl.dmu = (zedd - temp200) / sd
        dl.dmu[NN] = zedd[NN] / sd[NN]
        dl.dsd = (-1 + zedd^2 - zedd * temp200)  / sd
        dl.dsd[NN] = (-1 + zedd[NN]^2)  / sd[NN]
        dmu.deta = dtheta.deta(mymu, .lmean) 
        dsd.deta = dtheta.deta(sd, .lsd) 
        if(iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }
        derivnew = w * cbind(dl.dmu * dmu.deta, dl.dsd * dsd.deta)
        derivnew
    }), list( .lmean=lmean, .lsd=lsd ))),
    weight=expression({
        if(iter == 1) {
            wznew = cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
        } else {
            wzold = wznew
            wznew = qnupdate(w=w, wzold=wzold, dderiv=(derivold - derivnew),
                             deta=etanew-etaold, M=M,
                             trace=trace)  # weights incorporated in args
        }
        wznew
    }))
}



recexp1.control <- function(save.weight=TRUE, ...)
{
    list(save.weight=save.weight)
}

recexp1 = function(lrate="loge", irate=NULL, method.init=1)
{

    if(mode(lrate) != "character" && mode(lrate) != "name")
        lrate = as.character(substitute(lrate))
    if(!is.Numeric(method.init, allow=1, integ=TRUE, positi=TRUE) ||
       method.init > 3.5)
        stop("argument \"method.init\" must be 1 or 2 or 3")

    new("vglmff",
    blurb=c("Upper record values from a 1-parameter exponential distribution\n\n",
            "Links:    ",
            namesof("rate", lrate, tag= TRUE),
            "\n",
            "Variance: 1/rate^2"),
    initialize=eval(substitute(expression({
        predictors.names = c(namesof("rate", .lrate, tag= FALSE))
        if(ncol(y <- cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(any(diff(y) <= 0))
            stop("response must have increasingly larger and larger values")
        if(any(w != 1))
            warning("weights should have unit values only")
        if(!length(etastart)) {
            rate.init = if(length( .irate)) rep( .irate, len=n) else {
                init.rate =
                    if(.method.init == 1) length(y) / y[length(y),1] else
                    if(.method.init == 2) 1/mean(y) else 1/median(y)
                if(.lrate == "loge") pmax(1/1024, init.rate) else init.rate}
            etastart = cbind(theta2eta(rep(rate.init, len=n), .lrate))
        }
    }), list( .lrate=lrate, .irate=irate, .method.init=method.init ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        eta2theta(eta, .lrate)
    }, list( .lrate=lrate ))),
    last=eval(substitute(expression({
        misc$link = c("rate"= .lrate)
        misc$expected = TRUE
    }), list( .lrate=lrate ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        rate = eta2theta(eta, .lrate)
        if(residuals) stop("loglikelihood residuals not implemented yet") else {
            NN = length(eta)
            y = cbind(y)
            sum(w * log(rate)) - w[NN] * rate[NN] * y[NN,1]
        }
    }, list( .lrate=lrate ))),
    vfamily=c("recexp1"),
    deriv=eval(substitute(expression({
        NN = length(eta)
        rate = c(eta2theta(eta, .lrate))
        dl.drate = 1 / rate 
        dl.drate[NN] = 1/ rate[NN] - y[NN,1]
        drate.deta = dtheta.deta(rate, .lrate)
        w * cbind(dl.drate * drate.deta)
    }), list( .lrate=lrate ))),
    weight=expression({
        ed2l.drate2 = -1 / rate^2
        wz = -w * drate.deta^2 * ed2l.drate2
        wz
    }))
}



if(FALSE)
exbilogi <- function(zero=c(3,6,7),
                   llocation="identity",
                   lscale="loge",
                   lshape="identity",
                   ldependency="logit",
                   percentiles=c(95,99),
                   ishape1=-0.1, ishape2=-0.1, idependency=1)
{
 stop("not working")

    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))

    if(mode(llocation) != "character" && mode(llocation) != "name")
        llocation = as.character(substitute(llocation))

    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))

    if(mode(ldependency) != "character" && mode(ldependency) != "name")
        ldependency = as.character(substitute(ldependency))
    if(!is.Numeric(percentiles, posit=TRUE) || max(percentiles) >= 100)
        stop("bad input for argument \"percentiles\"")

    new("vglmff",
    blurb=c("Extreme Bivariate Logistic Distribution\n",
            "Links:    ",
            namesof("loc1", link=llocation), ", ", 
            namesof("scale1", link=lscale), ", ",
            namesof("shape1", link=lshape), ", ",
            namesof("loc2", link=llocation), ", ", 
            namesof("scale2", link=lscale), ", ",
            namesof("shape2", link=lshape), ", ",
            namesof("dependency", link=ldependency)),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        predictors.names = c(
        namesof("loc1",  .llocation,  short=TRUE),
        namesof("scale1",  .lscale,  short=TRUE),
        namesof("shape1",  .lshape,  short=TRUE),
        namesof("loc2",  .llocation,  short=TRUE),
        namesof("scale2",  .lscale,  short=TRUE),
        namesof("shape2",  .lshape,  short=TRUE),
        namesof("dependency",  .ldependency,  short=TRUE))
        if(ncol(as.matrix(y)) != 2)
            stop("the response must be a two-column matrix")

        if(!length(etastart)) {
            EulerM = -digamma(1)
            # Initial values for limiting case as xi --> 0, r_i==1
            sig1.init = sqrt(6 * var(y[,1]))/pi
            mu1.init = median(y[,1]) - EulerM * sig1.init
            xi1.init = if(length(.ishape1)) .ishape1 else 0.1
            sig2.init = sqrt(6 * var(y[,2]))/pi
            mu2.init = median(y[,2]) - EulerM * sig2.init
            xi2.init = if(length(.ishape2)) .ishape2 else 0.1
            alpha.init = if(length(.idependency)) .idependency else 1
            etastart = cbind(rep(theta2eta(mu1.init, .llocation), nrow(y)),
                                 theta2eta(sig1.init, .lscale), 
                                 theta2eta(xi1.init,  .lshape),
                             rep(theta2eta(mu2.init, .llocation), nrow(y)),
                                 theta2eta(sig2.init, .lscale), 
                                 theta2eta(xi2.init,  .lshape),
                                 theta2eta(alpha.init,  .ldependency))
            bad1 = 1 + xi1.init*(y[,1]-mu1.init)/sig1.init <= 0
            bad2 = 1 + xi2.init*(y[,2]-mu2.init)/sig2.init <= 0
            if(any(bad1, na.rm = TRUE) || any(bad2, na.rm = TRUE))
                stop(paste(sum(bad1, na.rm = TRUE) + sum(bad2, na.rm = TRUE),
                           "observations violating",
                           "boundary constraints while initializing"))
        }
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
             .ldependency=ldependency,
             .ishape1=ishape1, .ishape2=ishape2, .idependency=idependency ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        mmus = eta2theta(eta[,c(1,4)], .llocation)
        sigmas = eta2theta(eta[,c(2,5)], .lscale)
        xis = eta2theta(eta[,c(3,6)], .lshape)
        alpha = eta2theta(eta[,7], .ldependency)
        cent = .percentiles
        lp = length(cent)
        fv = matrix(as.numeric(NA), nrow(eta), lp)
        dimnames(fv) = list(dimnames(eta)[[1]],
                             paste(as.character(.percentiles), "%", sep=""))
        fv
    }, list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
            .ldependency=ldependency, .percentiles=percentiles ))),
    last=eval(substitute(expression({
        misc$links = c(loc1 = .llocation, lscale1 = .lscale, shape1 =.lshape,
                       loc2 = .llocation, lscale2 = .lscale, shape2 =.lshape,
                       dependency= .ldependency)
        misc$true.mu = FALSE    # @fitted.value is not a true mu
        misc$percentiles = .percentiles
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
             .ldependency=ldependency, .percentiles=percentiles ))),
    loglikelihood=eval(substitute(
        function(mu,y,w, residuals=FALSE, eta, extra=NULL) {
        mmus = eta2theta(eta[,c(1,4)], .llocation)
        sigmas = eta2theta(eta[,c(2,5)], .lscale)
        xis = eta2theta(eta[,c(3,6)], .lshape)
        alpha = eta2theta(eta[,7], .ldependency)
        zedds = (y - mmus) / sigmas
        V = (zedds[,1]^(1/alpha) + zedds[,2]^(1/alpha))^alpha
        if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * (log1p(-alpha) - log(alpha) -
                     log(sigmas[,1]) - log(sigmas[,2]) +
                     (xis[,1] + 1/alpha) * log(zedds[,1]) +
                     (xis[,2] + 1/alpha) * log(zedds[,2]) +
                     (1 - 2/alpha) * log(V) - V))
        }, list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
                .ldependency=ldependency ))),
    vfamily=c("exbilogi", "vextremes"),
    deriv=eval(substitute(expression({
        mmus = eta2theta(eta[,c(1,4)], .llocation)
        sigmas = eta2theta(eta[,c(2,5)], .lscale)
        xis = eta2theta(eta[,c(3,6)], .lshape)
        alpha = eta2theta(eta[,7], .ldependency)
        zedds = (y - mmus) / sigmas
        Zee = (1 + zedds)^(-1/xis)
        dZee.dmus = 
        dZee.dsigmas = 
        dZee.dxis = 
        A = (zedds[,1]^(1/alpha) + zedds[,2]^(1/alpha))
        V = A^alpha
        dV.dalpha = V * (log(A) - (zedds[,1]^(1/alpha) *log(zedds[,1]) +
                         zedds[,2]^(1/alpha) *log(zedds[,2])) / (alpha*A))
        temp400 = (1/alpha + xis)
        temp300 = ((1-2/alpha)/V - 1)

        dzedds.dmus =  Zee^(1+xis) / sigmas
        dzedds.dsigmas = (y- zedds) * Zee^(1+xis) / sigmas
        dzedds.dxis = Zee * (log(Zee)/xis - zedds/(1+xis*(y-mmmu)/sigmas)) / xis
        dV.dmus = alpha * V^(1-1/alpha) * dZee.dmus
        dV.dsigmas =  alpha * V^(1-1/alpha) * dZee.dsigmas
        dV.dxis =  alpha * V^(1-1/alpha) * dZee.dxis

        dl.dalpha = -1/(1-alpha) - 1/alpha - (log(zedds[,1]) + log(zedds[,2])) /
                    alpha^2 + 2*log(V)/alpha^2 + temp300 * dV.dalpha
        dl.dmus = temp400 * dzedds.dmus / zedds + temp300 * dV.dmus
        dl.dsigmas = -1/sigmas + temp400 * dzedds.dsigmas / zedds +
                     temp300 * dV.dsigmas
        dl.dxis = log(zedds) + temp400 * dzedds.dxis / zedds + temp300 * dV.dxis
        dmus.dtheta = dtheta.deta(mmus, .llocation)
        dsigmas.dtheta = dtheta.deta(sigmas, .lscale)
        dxis.dtheta = dtheta.deta(xis, .lshape)

        if(iter == 1) {
            etanew = eta
        } else {
            derivold = derivnew
            etaold = etanew
            etanew = eta
        }
        derivnew = w * cbind(dl.dmus[,1]    * dmus.dtheta[,1],
                             dl.dsigmas[,1] * dsigmas.dtheta[,1],
                             dl.dxis[,1]    * dxi.dtheta[,1],
                             dl.dmus[,2]    * dmus.dtheta[,2],
                             dl.dsigmas[,2] * dsigmas.dtheta[,2],
                             dl.dxis[,2]    * dxi.dtheta[,2],
                             dl.dalpha      * dtheta.deta(alpha, .ldependency))
        derivnew
    }), list( .llocation=llocation, .lscale=lscale, .lshape=lshape,
             .ldependency=ldependency ))),
    weight=eval(substitute(expression({
        if(iter == 1) {
            wznew = cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
        } else {
            wzold = wznew
            wznew = qnupdate(w=w, wzold=wzold, dderiv=(derivold - derivnew),
                             deta=etanew-etaold, M=M,
                             trace=trace)  # weights incorporated in args
        }
        wznew
    }), list( .llocation=llocation ))))
}




poissonp = function(ostatistic, dimension=2, link="loge", earg=list(),
                    idensity=NULL, method.init=1) {
    if(!is.Numeric(ostatistic, posit=TRUE, allow=1, integ=TRUE))
        stop("argument 'ostatistic' must be a single positive integer")
    if(!is.Numeric(dimension, posit=TRUE, allow=1, integ=TRUE) ||
       dimension > 3)
        stop("argument 'dimension' must be 2 or 3")
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()
    if(!is.Numeric(method.init, allow=1, posit=TRUE, integer=TRUE) ||
       method.init > 2.5)
        stop("argument \"method.init\" must be 1 or 2")
    if(length(idensity) && !is.Numeric(idensity, posit=TRUE))
        stop("bad input for argument \"idensity\"")

    new("vglmff",
    blurb=c(if(dimension==2)
            "Poisson-points-on-a-plane distances distribution\n" else
            "Poisson-points-on-a-volume distances distribution\n",
            "Link:    ",
            namesof("density", link, earg=earg), "\n\n",
            if(dimension==2)
            "Mean:    gamma(s+0.5) / (gamma(s) * sqrt(density * pi))" else
            "Mean:    gamma(s+1/3) / (gamma(s) * (4*density*pi/3)^(1/3))"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(any(y <= 0))
            stop("response must contain positive values only")
        predictors.names = namesof("density", .link, earg=.earg, tag=FALSE) 
        if(!length(etastart)) {
            use.this = if( .method.init == 1) median(y) + 1/8 else
                       weighted.mean(y,w)
            if( .dimension == 2) {
                myratio = exp(lgamma( .ostatistic +0.5) - lgamma( .ostatistic))
                density.init = if(is.Numeric( .idensity))
                    rep( .idensity, len=n) else
                    rep(myratio^2 / (pi * use.this^2), len=n)
                etastart = theta2eta(density.init, .link, earg= .earg)
            } else {
                myratio = exp(lgamma( .ostatistic +1/3) - lgamma( .ostatistic))
                density.init = if(is.Numeric( .idensity))
                    rep( .idensity, len=n) else
                    rep(3 * myratio^3 / (4 * pi * use.this^3), len=n)
                etastart = theta2eta(density.init, .link, earg= .earg)
            }
        }
    }), list( .link=link, .earg=earg, .ostatistic=ostatistic,
              .dimension=dimension, .method.init=method.init,
              .idensity=idensity ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        density = eta2theta(eta, .link, earg= .earg)
        if( .dimension == 2) {
            myratio = exp(lgamma( .ostatistic +0.5) - lgamma( .ostatistic ))
            myratio / sqrt(density * pi)
        } else {
            myratio = exp(lgamma( .ostatistic +1/3) - lgamma( .ostatistic))
            myratio / (4*density * pi/3)^(1/3)
        }
    }, list( .link=link, .earg=earg, .ostatistic=ostatistic,
             .dimension=dimension ))),
    last=eval(substitute(expression({
        misc$link = c("density"= .link)
        misc$earg = list("density"= .earg)
        misc$expected = TRUE
        misc$ostatistic = .ostatistic
        misc$dimension = .dimension
    }), list( .link=link, .earg=earg, .ostatistic=ostatistic,
              .dimension=dimension ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        density = eta2theta(eta, .link, earg= .earg)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
            if( .dimension == 2)
                sum(w * (log(2) + .ostatistic * log(pi * density) -
                     lgamma( .ostatistic) + (2* .ostatistic-1) * log(y) -
                     density * pi * y^2)) else
                sum(w * (log(3) + .ostatistic * log(4*pi * density/3) -
                     lgamma( .ostatistic) + (3* .ostatistic-1) * log(y) -
                     (4/3) * density * pi * y^3))
    }, list( .link=link, .earg=earg, .ostatistic=ostatistic,
             .dimension=dimension ))),
    vfamily=c("poissonp"),
    deriv=eval(substitute(expression({
        density = eta2theta(eta, .link, earg= .earg)
        if( .dimension == 2) {
            dl.ddensity = .ostatistic / density - pi * y^2
        } else {
            dl.ddensity = .ostatistic / density - (4/3) * pi * y^3
        }
        ddensity.deta = dtheta.deta(density, .link, earg= .earg)
        w * dl.ddensity * ddensity.deta
    }), list( .link=link, .earg=earg, .ostatistic=ostatistic,
              .dimension=dimension ))),
    weight=eval(substitute(expression({
        ed2l.ddensity2 = .ostatistic / density^2
        wz = ddensity.deta^2 * ed2l.ddensity2
        w * wz
    }), list( .link=link, .earg=earg, .ostatistic=ostatistic,
              .dimension=dimension ))))
}


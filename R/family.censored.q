# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.










if(FALSE)
cexpon = 
ecexpon = function(link="loge", location=0)
{
    if(!is.Numeric(location, allow=1))
        stop("bad input for \"location\"")
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))

    new("vglmff",
    blurb=c("Censored exponential distribution\n\n",
            "Link:     ", namesof("rate", link, tag= TRUE), "\n",
            "Mean:     ", "mu = ", location, " + 1 / ",
            namesof("rate", link, tag= FALSE), "\n",
            "Variance: ",
            if(location==0) "Exponential: mu^2" else
            paste("(mu-", location, ")^2", sep="")),
    initialize=eval(substitute(expression({
        extra$location = .location # This is passed into, e.g., link, deriv etc.
        if(any(y[,1] <= extra$location))
            stop(paste("all responses must be greater than", extra$location))
        predictors.names = namesof("rate", .link, tag= FALSE)
        type <- attr(y, "type")
        if (type=="right" || type=="left"){
          mu = y[,1] + (abs(y[,1] - extra$location) < 0.001) / 8
        }else
        if (type=="interval"){
          temp <- y[,3]
          mu = ifelse(temp == 3, y[,2] + (abs(y[,2] - extra$location) < 0.001)
          / 8,y[,1] + (abs(y[,1] - extra$location) < 0.001) / 8)
        }
        if(!length(etastart))
            etastart = theta2eta(1/(mu-extra$location), .link)

        if (type=="right") {
          temp <- y[, 2]
          extra$uncensored = ifelse(temp == 1, TRUE, FALSE)
          extra$rightcensored = ifelse(temp == 0, TRUE, FALSE)
          extra$leftcensored = rep(FALSE, len=n)
          extra$interval = rep(FALSE, len=n)
        } else
        if (type=="left") {
          temp <- y[, 2]
          extra$uncensored = ifelse(temp == 1, TRUE, FALSE)
          extra$rightcensored = rep(FALSE, len=n)
          extra$leftcensored = ifelse(temp == 0, TRUE, FALSE)
          extra$interval = rep(FALSE, len=n)
        } else
        if (type=="counting") {
          stop("type=='counting' not recognized")
          extra$uncensored = rep(temp == 1, TRUE, FALSE)
          extra$interval = rep(FALSE, len=n)
          extra$leftcensored = rep(FALSE, len=n)
          extra$rightcensored = rep(FALSE, len=n)
          extra$counting = ifelse(temp == 0, TRUE, FALSE)
        } else
        if (type=="interval") {
          temp <- y[, 3]
          extra$uncensored = ifelse(temp == 1, TRUE, FALSE)
          extra$rightcensored = ifelse(temp == 0, TRUE, FALSE)
          extra$leftcensored = ifelse(temp == 2, TRUE, FALSE)
          extra$interval = ifelse(temp == 3, TRUE, FALSE)
        } else
          stop("'type' not recognized")
        #if(!length(extra$leftcensored)) extra$leftcensored = rep(FALSE, len=n)
        #if(!length(extra$rightcensored)) extra$rightcensored = rep(FALSE, len=n)
        #if(any(extra$rightcensored & extra$leftcensored))
        #    stop("some observations are both right and left censored!")
    }), list( .location=location, .link=link ))),
    inverse=eval(substitute(function(eta, extra=NULL)
        extra$location + 1 / eta2theta(eta, .link),
    list( .link=link ) )),
    last=eval(substitute(expression({
        misc$location = extra$location
        misc$link = c("rate" = .link)
    }), list( .link=link ))),
    link=eval(substitute(function(mu, extra=NULL)
        theta2eta(1/(mu-extra$location), .link),
    list( .link=link ) )),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        rate = 1 / (mu - extra$location)
        cen0 = extra$uncensored
        cenL = extra$leftcensored
        cenU = extra$rightcensored
        cenI = extra$interval
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w[cenL] * log1p(-exp(-rate[cenL]*(y[cenL,1]-extra$location)))) +
        sum(w[cenU] * (-rate[cenU]*(y[cenU,1]-extra$location))) +
        sum(w[cen0] * (log(rate[cen0]) - rate[cen0]*(y[cen0,1]-extra$location)))+
        sum(w[cenI] * log(-exp(-rate[cenI]*(y[cenI,2]-extra$location))+
        exp(-rate[cenI]*(y[cenI,1]-extra$location))))
    }, list( .link=link ))),
    vfamily=c("ecexpon"),
    deriv=eval(substitute(expression({
        rate = 1 / (mu - extra$location)
        cen0 = extra$uncensored
        cenL = extra$leftcensored
        cenU = extra$rightcensored
        cenI = extra$interval
        dl.drate = 1/rate - (y[,1]-extra$location)  # uncensored
        tmp200 = exp(-rate*(y[,1]-extra$location))
        tmp200b = exp(-rate*(y[,2]-extra$location)) # for interval censored
        if(any(cenL))
            dl.drate[cenL] = (y[cenL,1]-extra$location) *
                             tmp200[cenL] / (1 - tmp200[cenL])
        if(any(cenU))
            dl.drate[cenU] = -(y[cenU,1]-extra$location)
        if(any(cenI))
            dl.drate[cenI] = ((y[cenI,2]-extra$location)*tmp200b[cenI]-
            (y[cenI,1]-extra$location)*tmp200[cenI])/
            (-tmp200b[cenI]+tmp200[cenI])
        drate.deta = dtheta.deta(rate, .link)
        w * dl.drate * drate.deta
    }), list( .link=link ) )),
    weight=eval(substitute(expression({
        A123 = ((mu-extra$location)^2) # uncensored d2l.drate2
        Lowpt = ifelse(cenL, y[,1], extra$location)
        Lowpt = ifelse(cenI, y[,1], Lowpt) #interval censored
        Upppt = ifelse(cenU, y[,1], Inf)
        Upppt = ifelse(cenI, y[,2], Upppt) #interval censored
        tmp300 = exp(-rate*(Lowpt - extra$location))
        d2l.drate2 = 0 * y[,1]
        ind50 = Lowpt > extra$location
        d2l.drate2[ind50] = (Lowpt[ind50]-extra$location)^2 *
                            tmp300[ind50] / (1-tmp300[ind50])
        d2l.drate2 = d2l.drate2 + (exp(-rate*(Lowpt-extra$location)) -
                                   exp(-rate*(Upppt-extra$location))) * A123
        wz = w * (drate.deta^2) * d2l.drate2
        wz
    }), list( .link=link ))))
}



cnormal1 = function(lmu="identity", lsd="loge", imethod=1, zero=2)
{
    if(mode(lmu) != "character" && mode(lmu) != "name")
        lmu = as.character(substitute(lmu))
    if(mode(lsd) != "character" && mode(lsd) != "name")
        lsd = as.character(substitute(lsd))
    if(!is.Numeric(imethod, allow=1, integer=TRUE, positi=TRUE) || imethod > 2)
        stop("imethod must be 1 or 2")

    new("vglmff",
    blurb=c("Censored univariate normal\n\n",
            "Links:    ", namesof("mu", lmu, tag= TRUE), "; ",
                          namesof("sd", lsd, tag= TRUE), "\n",
            "Conditional variance: sd^2"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        y = cbind(y)
        if(ncol(y)>1) stop("the response must be a vector or a 1-column matrix")

        if(!length(extra$leftcensored)) extra$leftcensored = rep(FALSE, len=n)
        if(!length(extra$rightcensored)) extra$rightcensored = rep(FALSE, len=n)
        if(any(extra$rightcensored & extra$leftcensored))
            stop("some observations are both right and left censored!")

        predictors.names = 
        c(namesof("mu", .lmu, tag= FALSE),
          namesof("sd", .lsd, tag= FALSE))
        if(!length(etastart)) {
            anyc = extra$leftcensored | extra$rightcensored
            i11 = if( .imethod == 1) anyc else FALSE  # can be all data
            junk=if(is.R()) lm.wfit(x=cbind(x[!i11,]),y=y[!i11],w=w[!i11]) else
                   lm.wfit(x=cbind(x[!i11,]), y=y[!i11], w=w[!i11],method="qr")
            sd.y.est = sqrt( sum(w[!i11] * junk$resid^2) / junk$df.residual )
            etastart = cbind(mu=y, rep(theta2eta(sd.y.est, .lsd), length=n))
            if(any(anyc)) etastart[anyc,1] = x[anyc,,drop=FALSE] %*% junk$coeff
        }
   }), list( .lmu=lmu, .lsd=lsd, .imethod=imethod ))),
    inverse=eval(substitute( function(eta, extra=NULL) {
        eta2theta(eta[,1], .lmu)
    }, list( .lmu=lmu ))),
    last=eval(substitute(expression({
        misc$link = c("mu"= .lmu, "sd"= .lsd)
        misc$expected = TRUE
    }), list( .lmu=lmu, .lsd=lsd ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        cenL = extra$leftcensored
        cenU = extra$rightcensored
        cen0 = !cenL & !cenU   # uncensored obsns
        mum = eta2theta(eta[,1], .lmu)
        sd = eta2theta(eta[,2], .lsd)
        Lower = ifelse(cenL, y, -Inf)
        Upper = ifelse(cenU, y,  Inf)
        ell1 = -log(sd[cen0]) - 0.5 * ((y[cen0] - mum[cen0])/sd[cen0])^2
        ell2 = log1p(-pnorm((mum[cenL] - Lower[cenL])/sd[cenL]))
        ell3 = log1p(-pnorm(( Upper[cenU] -  mum[cenU])/sd[cenU]))
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w[cen0] * ell1) + sum(w[cenL] * ell2) + sum(w[cenU] * ell3)
    }, list( .lmu=lmu, .lsd=lsd ))),
    vfamily=c("tobit"),
    deriv=eval(substitute(expression({
        cenL = extra$leftcensored
        cenU = extra$rightcensored
        cen0 = !cenL & !cenU   # uncensored obsns
        Lower = ifelse(cenL, y, -Inf)
        Upper = ifelse(cenU, y,  Inf)
        mum = eta2theta(eta[,1], .lmu)
        sd = eta2theta(eta[,2], .lsd)
        dl.dmu = (y-mum) / sd^2
        dl.dsd = (((y-mum)/sd)^2 - 1) / sd
        dmu.deta = dtheta.deta(mum, .lmu) 
        dsd.deta = dtheta.deta(sd, .lsd) 
        if(any(cenL)) {
            mumL = mum - Lower
            temp21L = mumL[cenL] / sd[cenL]
            PhiL = pnorm(temp21L)
            phiL = dnorm(temp21L)
            fred21 = phiL / (1 - PhiL)
            dl.dmu[cenL] = -fred21 / sd[cenL]
            dl.dsd[cenL] = mumL[cenL] * fred21 / sd[cenL]^2
            rm(fred21)
        }
        if(any(cenU)) {
            mumU = Upper - mum
            temp21U = mumU[cenU] / sd[cenU]
            PhiU = pnorm(temp21U)
            phiU = dnorm(temp21U)
            fred21 = phiU / (1 - PhiU)
            dl.dmu[cenU] = fred21 / sd[cenU]   # Negated
            dl.dsd[cenU] = mumU[cenU] * fred21 / sd[cenU]^2
            rm(fred21)
        }
        w * cbind(dl.dmu * dmu.deta, dl.dsd * dsd.deta)
    }), list( .lmu=lmu, .lsd=lsd ))),
    weight=eval(substitute(expression({
        A1 = 1 - pnorm((mum - Lower) / sd)   # Lower
        A3 = 1 - pnorm(( Upper - mum) / sd)  # Upper
        A2 = 1 - A1 - A3                      # Middle; uncensored
        wz = matrix(0, n, 3)
        wz[,iam(1,1,M)] = A2 * 1 / sd^2  # ed2l.dmu2
        wz[,iam(2,2,M)] = A2 * 2 / sd^2  # ed2l.dsd2
        mumL = mum - Lower
        temp21L = mumL / sd
        PhiL = pnorm(temp21L)
        phiL = dnorm(temp21L)
        temp31L = ((1-PhiL) * sd)^2 
        wz.cenL11 = phiL * (phiL - (1-PhiL)*temp21L) / temp31L
        wz.cenL22 = mumL * phiL * ((1-PhiL) * (2 - temp21L^2) +
                    mumL * phiL / sd) / (sd * temp31L)
        wz.cenL12 = phiL * ((1-PhiL)*(temp21L^2 - 1) - temp21L*phiL) / temp31L
        wz.cenL11[!is.finite(wz.cenL11)] = 0
        wz.cenL22[!is.finite(wz.cenL22)] = 0
        wz.cenL12[!is.finite(wz.cenL12)] = 0
        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] + A1 * wz.cenL11
        wz[,iam(2,2,M)] = wz[,iam(2,2,M)] + A1 * wz.cenL22
        wz[,iam(1,2,M)] = A1 * wz.cenL12
        mumU = Upper - mum    # often Inf
        temp21U = mumU / sd    # often Inf
        PhiU = pnorm(temp21U)  # often 1
        phiU = dnorm(temp21U)  # often 0
        temp31U = ((1-PhiU) * sd)^2  # often 0
        tmp8 = (1-PhiU)*temp21U
        wzcenU11 = phiU * (phiU - tmp8) / temp31U
        tmp9 = (1-PhiU) * (2 - temp21U^2)
        wzcenU22 = mumU * phiU * (tmp9 + mumU * phiU / sd) / (sd * temp31U)
        wzcenU12 = -phiU * ((1-PhiU)*(temp21U^2 - 1) - temp21U*phiU) / temp31U
        wzcenU11[!is.finite(wzcenU11)] = 0  # Needed when Upper==Inf
        wzcenU22[!is.finite(wzcenU22)] = 0  # Needed when Upper==Inf
        wzcenU12[!is.finite(wzcenU12)] = 0  # Needed when Upper==Inf
        wz[,iam(1,1,M)] = wz[,iam(1,1,M)] + A3 * wzcenU11
        wz[,iam(2,2,M)] = wz[,iam(2,2,M)] + A3 * wzcenU22
        wz[,iam(1,2,M)] = wz[,iam(1,2,M)] + A3 * wzcenU12
        wz[,iam(1,1,M)] = w * wz[,iam(1,1,M)] * dmu.deta^2
        wz[,iam(2,2,M)] = w * wz[,iam(2,2,M)] * dsd.deta^2
        wz[,iam(1,2,M)] = w * wz[,iam(1,2,M)] * dmu.deta * dsd.deta
        wz
    }), list( .lmu=lmu, .lsd=lsd ))))
}



crayleigh = function(link="loge", earg = list(), expected=FALSE) {
    if(mode(link) != "character" && mode(link) != "name")
        link = as.character(substitute(link))
    if(!is.logical(expected) || length(expected) != 1)
        stop("bad input for argument \"expected\"")
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Censored Rayleigh distribution",
            "f(y) = y*exp(-0.5*(y/a)^2)/a^2, y>0, a>0\n",
            "Link:    ",
            namesof("a", link, earg= earg ), "\n", "\n",
            "Mean:    a * sqrt(pi / 2)"),
    initialize=eval(substitute(expression({
        if(ncol(cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        if(length(extra$leftcensored)) stop("cannot handle left-censored data")
        if(!length(extra$rightcensored)) extra$rightcensored = rep(FALSE, len=n)
        predictors.names = namesof("a", .link, earg= .earg, tag= FALSE) 
        if(!length(etastart)) {
            a.init = (y+1/8) / sqrt(pi/2)
            etastart = theta2eta(a.init, .link, earg= .earg )
        }
    }), list( .link=link, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        a = eta2theta(eta, .link, earg= .earg )
        a * sqrt(pi/2)
    }, list( .link=link, .earg=earg ))),
    last=eval(substitute(expression({
        misc$link = c("a"= .link)
        misc$earg = list(a= .earg)
        misc$expected = .expected
    }), list( .link=link, .earg=earg, .expected=expected ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        a = eta2theta(eta, .link, earg= .earg )
        cen0 = !extra$rightcensored   # uncensored obsns
        cenU = extra$rightcensored
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w[cen0]*(log(y[cen0]) - 2*log(a[cen0]) - 0.5*(y[cen0]/a[cen0])^2)) -
        0.5 * sum(w[cenU] * (y[cenU]/a[cenU])^2)
    }, list( .link=link, .earg=earg ))),
    vfamily=c("crayleigh"),
    deriv=eval(substitute(expression({
        cen0 = !extra$rightcensored   # uncensored obsns
        cenU = extra$rightcensored
        a = eta2theta(eta, .link, earg= .earg )
        dl.da = ((y/a)^2 - 2) / a
        da.deta = dtheta.deta(a, .link, earg= .earg )
        dl.da[cenU] = y[cenU]^2 / a[cenU]^3
        w * dl.da * da.deta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        ed2l.da2 = 4 / a^2
        wz = da.deta^2 * ed2l.da2
        if( .expected) {
            ed2l.da2[cenU] = 6 / (a[cenU])^2
            wz[cenU] = (da.deta[cenU])^2 * ed2l.da2[cenU]
        } else {
            d2l.da2 = 3 * (y[cenU])^2 / (a[cenU])^4
            d2a.deta2 = d2theta.deta2(a[cenU], .link, earg= .earg )
            wz[cenU] = (da.deta[cenU])^2 * d2l.da2 - dl.da[cenU] * d2a.deta2
        }
        w * wz
    }), list( .link=link, .earg=earg, .expected=expected ))))
}


weibull = 
weibull.sev = function(lshape="loge", lscale="loge",
                       eshape=list(), escale=list(),
                       ishape=NULL, iscale=NULL,
                       nrfs = 1,
                       imethod=1, zero=2)
{

    if(mode(lshape) != "character" && mode(lshape) != "name")
        lshape = as.character(substitute(lshape))
    if(mode(lscale) != "character" && mode(lscale) != "name")
        lscale = as.character(substitute(lscale))
    if(length(zero) && !is.Numeric(zero, integer=TRUE, posit=TRUE))
        stop("bad input for argument \"zero\"")
    if(!is.Numeric(imethod, allow=1, integer=TRUE, positi=TRUE) || imethod > 2)
        stop("argument \"imethod\" must be 1 or 2")
    if(!is.list(eshape)) eshape = list()
    if(!is.list(escale)) escale = list()
    if(!is.Numeric(nrfs, allow=1) || nrfs<0 || nrfs > 1)
        stop("bad input for 'nrfs'")

    new("vglmff",
    blurb=c("Weibull distribution\n\n",
            "Links:    ",
            namesof("shape", lshape, earg= eshape), ", ", 
            namesof("scale", lscale, earg= escale), "\n", 
            "Mean:     scale * gamma(1 + 1/shape)\n",
            "Variance: scale^2 * (gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)"),
    constraints=eval(substitute(expression({
        constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        y = cbind(y)
        if(ncol(y)>1) stop("the response must be a vector or a 1-column matrix")

        if(is.SurvS4(y))
            stop("only uncensored observations are allowed; don't use Surv()")

        predictors.names =
        c(namesof("shape", .lshape, earg= .eshape, tag=FALSE),
          namesof("scale", .lscale, earg= .escale, tag=FALSE))
        if(!length(.ishape) || !length(.iscale)) {
            anyc = FALSE  # extra$leftcensored | extra$rightcensored
            i11 = if( .imethod == 1) anyc else FALSE  # can be all data
            qvec = c(.25, .5, .75)   # Arbitrary; could be made an argument
            init.shape = if(length( .ishape)) .ishape else 1
            xvec = log(-log1p(-qvec))
            fit0 = lsfit(x=xvec, y=log(quantile(y[!i11], qvec)))
        }

        if(!length(etastart)) {
            shape = rep(if(length(.ishape)) .ishape else 1/fit0$coef["X"],len=n)
            scale = rep(if(length(.iscale)) .iscale else
                        exp(fit0$coef["Intercept"]), len=n)
            etastart =
            cbind(theta2eta(shape, .lshape, earg= .eshape ),
                  theta2eta(scale, .lscale, earg= .escale ))
        }
    }), list( .lscale=lscale, .lshape=lshape,
              .escale=escale, .eshape=eshape,
              .iscale=iscale, .ishape=ishape, .imethod=imethod ) )),
    inverse=eval(substitute(function(eta, extra=NULL) {
        shape = eta2theta(eta[,1], .lshape, earg= .eshape )
        scale = eta2theta(eta[,2], .lscale, earg= .escale )
        scale * gamma(1+1/shape)
    }, list( .lscale=lscale, .lshape=lshape,
             .escale=escale, .eshape=eshape ) )),
    last=eval(substitute(expression({
        if(regnotok <- any(shape <= 2))
            warning(paste("MLE regularity conditions are violated",
                          "(shape <= 2) at the final iteration"))
        misc$link = c(shape= .lshape, scale= .lscale)
        misc$earg= list(shape= .eshape, scale= .escale)
        misc$nrfs = .nrfs
        misc$RegCondOK = !regnotok   # Save this for later
    }), list( .lscale=lscale, .lshape=lshape,
              .escale=escale, .eshape=eshape, .nrfs=nrfs ) )),
    loglikelihood=eval(substitute(
            function(mu, y, w, residuals= FALSE,eta, extra=NULL) {
        shape = eta2theta(eta[,1], .lshape, earg= .eshape )
        scale = eta2theta(eta[,2], .lscale, earg= .escale )
        ell1 = (log(shape) - log(scale) + (shape-1) *
               log(y/scale) - (y / scale)^shape)
        if(residuals) stop("loglikelihood residuals not implemented yet") else
            sum(w * ell1)
    }, list( .lscale=lscale, .lshape=lshape,
             .escale=escale, .eshape=eshape ) )),
    vfamily=c("weibull.sev"),
    deriv=eval(substitute(expression({
        shape = eta2theta(eta[,1], .lshape, earg= .eshape )
        scale = eta2theta(eta[,2], .lscale, earg= .escale )
        dl.dshape = 1/shape + log(y/scale) - log(y/scale) * (y/scale)^shape
        dl.dscale = (shape/scale) * (-1 + (y/scale)^shape)
        dshape.deta = dtheta.deta(shape, .lshape, earg= .eshape )
        dscale.deta = dtheta.deta(scale, .lscale, earg= .escale )
        w * cbind( dl.dshape * dshape.deta, dl.dscale * dscale.deta )
    }), list( .lscale=lscale, .lshape=lshape,
              .escale=escale, .eshape=eshape ) )),
    weight=eval(substitute(expression({
        EulerM = -digamma(1.0)
        wz = matrix(as.numeric(NA), n, dimm(M))  #3=dimm(M)
        ed2l.dshape = (6*(EulerM-1)^2 +pi^2)/(6*shape^2) # Kleiber & Kotz (2003)
        ed2l.dscale = (shape/scale)^2
        ed2l.dshapescale = (EulerM-1)/scale
        wz[,iam(1,1,M)] = ed2l.dshape * dshape.deta^2
        wz[,iam(2,2,M)] = ed2l.dscale * dscale.deta^2
        wz[,iam(1,2,M)] = ed2l.dshapescale * dscale.deta * dshape.deta
        wz = w * wz
        wz
    }), list( .eshape=eshape, .nrfs=nrfs ))))
}






setOldClass(c("SurvS4","Surv"))


Surv <-
function (time, time2, event, type = c("right", "left", "interval",
    "counting", "interval2"), origin = 0)
{
    nn <- length(time)
    ng <- nargs()
    if (missing(type)) {
        if (ng == 1 || ng == 2)
            type <- "right" else if (ng == 3)
            type <- "counting" else stop("Invalid number of arguments")
    } else {
        type <- match.arg(type)
        ng <- ng - 1
        if (ng != 3 && (type == "interval" || type == "counting"))
            stop("Wrong number of args for this type of survival data")
        if (ng != 2 && (type == "right" || type == "left" ||
            type == "interval2"))
            stop("Wrong number of args for this type of survival data")
    }
    who <- !is.na(time)
    if (ng == 1) {
        if (!is.numeric(time))
            stop("Time variable is not numeric")
        ss <- cbind(time, 1)
        dimnames(ss) <- list(NULL, c("time", "status"))
    } else if (type == "right" || type == "left") {
        if (!is.numeric(time))
            stop("Time variable is not numeric")
        if (length(time2) != nn)
            stop("Time and status are different lengths")
        if (is.logical(time2))
            status <- 1 * time2 else if (is.numeric(time2)) {
            who2 <- !is.na(time2)
            if (max(time2[who2]) == 2)
                status <- time2 - 1 else status <- time2
            if (any(status[who2] != 0 & status[who2] != 1))
                stop("Invalid status value")
        } else stop("Invalid status value")
        ss <- cbind(time, status)
        dimnames(ss) <- list(NULL, c("time", "status"))
    } else if (type == "counting") {
        if (length(time2) != nn)
            stop("Start and stop are different lengths")
        if (length(event) != nn)
            stop("Start and event are different lengths")
        if (!is.numeric(time))
            stop("Start time is not numeric")
        if (!is.numeric(time2))
            stop("Stop time is not numeric")
        who3 <- who & !is.na(time2)
        if (any(time[who3] >= time2[who3]))
            stop("Stop time must be > start time")
        if (is.logical(event))
            status <- 1 * event else if (is.numeric(event)) {
            who2 <- !is.na(event)
            if (max(event[who2]) == 2)
                status <- event - 1 else status <- event
            if (any(status[who2] != 0 & status[who2] != 1))
                stop("Invalid status value")
        } else stop("Invalid status value")
        ss <- cbind(time - origin, time2 - origin, status)
    } else {
        if (type == "interval2") {
            event <- ifelse(is.na(time), 2, ifelse(is.na(time2),
                0, ifelse(time == time2, 1, 3)))
            if (any(time[event == 3] > time2[event == 3]))
                stop("Invalid interval: start > stop")
            time <- ifelse(event != 2, time, time2)
            type <- "interval"
        } else {
            temp <- event[!is.na(event)]
            if (!is.numeric(temp))
                stop("Status indicator must be numeric")
            if (length(temp) > 0 && any(temp != floor(temp) |
                temp < 0 | temp > 3))
                stop("Status indicator must be 0, 1, 2 or 3")
        }
        status <- event
        ss <- cbind(time, ifelse(!is.na(event) & event == 3,
            time2, 1), status)
    }
    attr(ss, "type") <- type
    class(ss) <- "SurvS4"
    ss
}



is.SurvS4 <- function(x) inherits(x, "SurvS4")




as.character.SurvS4 <-
function (x, ...)
{
    class(x) <- NULL
    type <- attr(x, "type")

    if (type == "right") {
        temp <- x[, 2]
        temp <- ifelse(is.na(temp), "?", ifelse(temp == 0, "+", " "))
        paste(format(x[, 1]), temp, sep = "")
    } else if (type == "counting") {
        temp <- x[, 3]
        temp <- ifelse(is.na(temp), "?", ifelse(temp == 0, "+", " "))
        paste("(", format(x[, 1]), ",", format(x[, 2]), temp, "]", sep = "")
    } else if (type == "left") {
        temp <- x[, 2]
        temp <- ifelse(is.na(temp), "?", ifelse(temp == 0, "<", " "))
        paste(temp, format(x[, 1]), sep = "")
    } else {
        stat <- x[, 3]
        temp <- c("+", "", "-", "]")[stat + 1]
        temp2 <- ifelse(stat == 3, paste("[", format(x[, 1]),
            ", ", format(x[, 2]), sep = ""), format(x[, 1]))
        ifelse(is.na(stat), as.character(NA), paste(temp2, temp, sep = ""))
    }
}



"[.SurvS4" <- function(x, i,j, drop=FALSE) {
    if (missing(j)) {
        temp <- class(x)
        type <- attr(x, "type")
        class(x) <- NULL
        x <- x[i, , drop=FALSE]
        class(x) <- temp
        attr(x, "type") <- type
        x
    } else {

        class(x) <- NULL
        NextMethod("[")
    }
}

is.na.SurvS4 <- function(x) {
    as.vector( (1* is.na(unclass(x)))%*% rep(1, ncol(x)) >0)
}





print.SurvS4 <-
function (x, quote = FALSE, ...)
invisible(print(as.character.SurvS4(x), quote = quote, ...))


setMethod("print", "SurvS4",
         function(x, ...)
         invisible(print.SurvS4(x, ...)))

setMethod("show", "SurvS4",
         function(object)
         invisible(print.SurvS4(object)))












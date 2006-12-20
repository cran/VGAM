# These functions are
# Copyright (C) 1998-2006 T.W. Yee, University of Auckland. All rights reserved.







dzipois = function(x, lambda, phi=0) {
    L = max(length(x), length(lambda), length(phi))
    x = rep(x, len=L); lambda = rep(lambda, len=L); phi = rep(phi, len=L);
    ans = dpois(x, lambda)
    if(!is.Numeric(phi) || any(phi < 0) || any(phi > 1))
        stop("phi must be between 0 and 1 inclusive")
    ifelse(x==0, phi + (1-phi) * ans, (1-phi) * ans)
}

pzipois = function(q, lambda, phi=0) {
    ans = ppois(q, lambda)
    phi = rep(phi, length=length(ans))
    if(!is.Numeric(phi) || any(phi < 0) || any(phi > 1))
        stop("phi must be between 0 and 1 inclusive")
    phi + (1-phi) * ans
}

qzipois = function(p, lambda, phi=0) {
    nn = max(length(p), length(lambda), length(phi))
    p = rep(p, len=nn)
    lambda = rep(lambda, len=nn)
    phi = rep(phi, len=nn)
    if(!is.Numeric(phi) || any(phi < 0) || any(phi > 1))
        stop("phi must be between 0 and 1 inclusive")
    ans = p 
    ans[p<=phi] = 0 
    ans[p>phi] = qpois((p[p>phi]-phi[p>phi])/(1-phi[p>phi]), lam=lambda[p>phi])
    ans
}

rzipois = function(n, lambda, phi=0) {
    if(!is.Numeric(n, positive=TRUE, integer=TRUE, allow=1))
        stop("n must be a single positive integer")
    ans = rpois(n, lambda)
    phi = rep(phi, len=length(ans))
    if(!is.Numeric(phi) || any(phi < 0) || any(phi > 1))
        stop("phi must be between 0 and 1 inclusive")
    ifelse(runif(n) < phi, 0, ans)
}


yip88 = function(link.lambda="loge", n.arg=NULL)
{
    if(mode(link.lambda) != "character" && mode(link.lambda) != "name")
        link.lambda = as.character(substitute(link.lambda))

    new("vglmff",
    blurb=c("Zero-inflated Poisson (based on Yip (1988))\n\n",
           "Link:     ", namesof("lambda", link.lambda), "\n",
           "Variance: (1-phi)*lambda"),
    first=eval(substitute(expression({
        zero <- y==0
        if(any(zero)) {
            if(length(extra)) extra$sumw = sum(w) else
                extra = list(sumw=sum(w))
            if(is.numeric(.n.arg) && extra$sumw != .n.arg) 
                stop(paste("value of n.arg conflicts with data",
                           "(it need not be specified anyway)"))
            warning("trimming out the zero observations")

            axa.save =  attr(x, "assign")
            x = x[!zero,,drop=FALSE]
            attr(x, "assign") = axa.save    # Don't lose these!!
            w = w[!zero]
            y = y[!zero]
        } else 
            if(!is.numeric(.n.arg)) 
                stop("n.arg must be supplied")
        
    }), list( .n.arg=n.arg ))),
    initialize=eval(substitute(expression({
        narg = if(is.numeric(.n.arg)) .n.arg else extra$sumw
        if(sum(w) > narg)
            stop("sum(w) > narg")

        predictors.names = namesof("lambda", .link.lambda, tag=FALSE)
        if(!length(etastart)) {
            lambda.init = rep(median(y), length=length(y))
            etastart = theta2eta(lambda.init, .link.lambda)
        }
        if(length(extra)) {
            extra$sumw = sum(w)
            extra$narg = narg   # For @inverse
        } else 
            extra = list(sumw=sum(w), narg = narg)
    }), list( .link.lambda=link.lambda, .n.arg=n.arg ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        lambda = eta2theta(eta, .link.lambda)
        temp5 = exp(-lambda)
        phi = (1 - temp5 - extra$sumw/extra$narg) / (1 - temp5)
        if(any(phi) <= 0)
            stop("non-positive value(s) of phi")
        (1-phi) * lambda
    }, list( .link.lambda=link.lambda ))),
    last=eval(substitute(expression({
        misc$link = c(lambda = .link.lambda)

        if(ncol(x)==1 && dimnames(x)[[2]]=="(Intercept)") {
            suma = extra$sumw
            phi = (1 - temp5[1] - suma/narg) / (1 - temp5[1])
            phi = if(phi < 0 || phi>1) NA else phi  # phi is a probability
            misc$phi = phi    # zz call it $p0 = phi ?? 
        }
    }), list( .link.lambda=link.lambda ))),
    loglikelihood=eval(substitute( 
        function(mu,y,w,residuals=FALSE, eta, extra=NULL) {
        lambda = eta2theta(eta, .link.lambda)
        lstar = -lambda + y * log(lambda) - log(1-exp(-lambda))
        sum(w * lstar)
    }, list( .link.lambda=link.lambda ))),
    vfamily=c("yip88"),
    deriv=eval(substitute(expression({
        lambda = eta2theta(eta, .link.lambda)
        temp5 = exp(-lambda)
        dl.dlambda = -1 + y/lambda - temp5/(1-temp5)
        dlambda.deta = dtheta.deta(lambda, .link.lambda)
        w * dl.dlambda * dlambda.deta
    }), list( .link.lambda=link.lambda ))),
    weight=eval(substitute(expression({
        d2lambda.deta2 = d2theta.deta2(lambda, .link.lambda)
        d2l.dlambda2 = -y / lambda^2 + temp5 / (1-temp5)^2
        -w * (d2l.dlambda2*dlambda.deta^2 + dl.dlambda*d2lambda.deta2)
    }), list( .link.lambda=link.lambda ))))
}




zapoisson = function(lp0="logit", llambda="loge",
                     ep0=list(), elambda=list())
{
    if(mode(lp0) != "character" && mode(lp0) != "name")
        lp0 = as.character(substitute(lp0))
    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(!is.list(ep0)) ep0 = list()
    if(!is.list(elambda)) elambda = list()

    new("vglmff",
    blurb=c(
  "Zero-altered Poisson (binomial and positive-Poisson conditional model)\n\n",
           "Links:    ",
           namesof("p0", lp0, earg=ep0, tag=FALSE), ", ",
           namesof("lambda", llambda, earg= .elambda, tag=FALSE),
           "\n"),
    initialize=eval(substitute(expression({
        y = as.matrix(y)
        extra$y0 = y0 = ifelse(y==0, 1, 0)
        extra$ymat = ymat = cbind(y0=y0, y=y)
        extra$NOS = NOS = ncoly = ncol(y)  # Number of species
        extra$skip.these = skip.these = matrix(as.logical(y0), n, NOS)

        mynames1 = if(ncoly==1) "p0" else paste("p0", 1:ncoly, sep="")
        mynames2 = if(ncoly==1) "lambda" else paste("lambda", 1:ncoly, sep="")
        predictors.names = 
            c(namesof(mynames1, .lp0, earg= .ep0, tag=FALSE),
              namesof(mynames2, .llambda, earg= .elambda, tag=FALSE))
        if(!length(etastart)) {
            etastart = cbind(theta2eta((0.5+w*y0)/(1+w), .lp0, earg= .ep0 ),
                             matrix(1, n, NOS))  # 1 here is any old value
            for(spp. in 1:NOS)
                etastart[!skip.these[,spp.],NOS+spp.] =
                    theta2eta(y[!skip.these[,spp.],spp.] /
                              (1-exp(-y[!skip.these[,spp.],spp.])), .llambda,
                              earg= .elambda )
        }
    }), list( .lp0=lp0, .llambda=llambda, .ep0= ep0, .elambda= elambda ))), 
    inverse=eval(substitute(function(eta, extra=NULL) {
        NOS = extra$NOS
        p0 = eta2theta(eta[,1:NOS], .lp0, earg= .ep0)
        lambda = eta2theta(eta[,NOS+(1:NOS)], .llambda, earg= .elambda)
        (1-p0) * (lambda / (1-exp(-lambda)))
    }, list( .lp0=lp0, .llambda=llambda, .ep0= ep0, .elambda= elambda ))),
    last=eval(substitute(expression({
        misc$link = c(rep( .lp0, len=NOS), rep( .llambda, len=NOS))
        names(misc$link) = c(mynames1, mynames2)
        misc$earg = vector("list", 2*NOS)
        names(misc$earg) = c(mynames1, mynames2)
        for(ii in 1:NOS) {
            misc$earg[[      ii]] = .ep0
            misc$earg[[NOS + ii]] = .elambda
        }
    }), list( .lp0=lp0, .llambda=llambda, .ep0= ep0, .elambda= elambda ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals=FALSE, eta,extra=NULL) {
        NOS = extra$NOS
        p0 = cbind(eta2theta(eta[,1:NOS], .lp0, earg= .ep0))
        skip = extra$skip.these
        lambda = cbind(eta2theta(eta[,NOS+(1:NOS)], .llambda, earg= .elambda ))
        ans = 0
        for(spp. in 1:NOS) {
            ans = ans + sum(w[skip[,spp.]] * log(p0[skip[,spp.],spp.])) +
                  sum(w[!skip[,spp.]] * (log(1-p0[!skip[,spp.],spp.]) -
                      log(1-exp(-lambda[!skip[,spp.],spp.])) -
                      lambda[!skip[,spp.],spp.] +
                      y[!skip[,spp.],spp.]*log(lambda[!skip[,spp.],spp.])))
        }
        ans
    }, list( .lp0=lp0, .llambda=llambda, .ep0= ep0, .elambda= elambda ))),
    vfamily=c("zapoisson"),
    deriv=eval(substitute(expression({
        NOS = extra$NOS
        y0 = extra$y0
        skip = extra$skip.these
        p0 = cbind(eta2theta(eta[,1:NOS], .lp0, earg= .ep0))
        lambda = cbind(eta2theta(eta[,NOS+(1:NOS)], .llambda, earg= .ep0))
        dl.dlambda = y/lambda - 1 - 1/(exp(lambda)-1)
        for(spp. in 1:NOS)
            dl.dlambda[skip[,spp.],spp.] = 0
        dlambda.deta = dtheta.deta(lambda, .llambda, earg= .ep0)
        mup0 = p0
        temp3 = if(.lp0 == "logit") {
            w * (y0 - mup0)
        } else
            w * dtheta.deta(mup0, link=.lp0, earg= .ep0) * (y0/mup0 - 1) / (1-mup0)
        ans = cbind(temp3, w * dl.dlambda * dlambda.deta)
        ans
    }), list( .lp0=lp0, .llambda=llambda, .ep0= ep0, .elambda= elambda ))),
    weight=eval(substitute(expression({
        wz = matrix( .Machine$double.eps^0.8, n, 2*NOS)
        for(spp. in 1:NOS) {
            temp4 = exp(lambda[!skip[,spp.], spp.])
            ed2l.dlambda2 = -temp4 * (1/lambda[!skip[,spp.],spp.] -
                            1/(temp4-1)) / (temp4-1)
            wz[!skip[,spp.],NOS+spp.] = -w[!skip[,spp.]] *
                                      (dlambda.deta[!skip[,spp.],spp.]^2) *
                                      ed2l.dlambda2
        }

        tmp100 = mup0*(1-mup0)
        tmp200 = if(.lp0 == "logit") {
            cbind(w * tmp100)
        } else {
            cbind(w * dtheta.deta(mup0, link= .lp0, earg= .ep0)^2 / tmp100)
        }
        for(ii in 1:NOS) {
            index200 = abs(tmp200[,ii]) < .Machine$double.eps
            if(any(index200)) {
                tmp200[index200,ii] = .Machine$double.eps # Diagonal 0's are bad 
            }
        }
        wz[,1:NOS] =  tmp200
        wz
    }), list( .lp0=lp0, .llambda=llambda, .ep0= ep0, .elambda= elambda ))))
}



zanegbinomial = function(lp0="logit", lmunb = "loge", lk = "loge",
                         ep0=list(), emunb =list(), ek = list(),
                         ik = 1, zero = -3, cutoff = 0.995, method.init=3)
{

    if(!is.Numeric(cutoff, positiv=TRUE, allow=1) || cutoff<0.8 || cutoff>=1)
        stop("range error in the argument cutoff")
    if(!is.Numeric(ik, positiv=TRUE))
        stop("ik must contain positive values only")
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 3) stop("argument \"method.init\" must be 1, 2 or 3")

    if(mode(lmunb) != "character" && mode(lmunb) != "name")
        lmunb = as.character(substitute(lmunb))
    if(mode(lk) != "character" && mode(lk) != "name")
        lk = as.character(substitute(lk))
    if(mode(lp0) != "character" && mode(lp0) != "name")
        lp0 = as.character(substitute(lp0))
    if(!is.list(ep0)) ep0 = list()
    if(!is.list(emunb)) emunb = list()
    if(!is.list(ek)) ek = list()

    new("vglmff",
    blurb=c("Zero-altered negative binomial (binomial and\n",
            "positive-negative binomial conditional model)\n\n",
           "Links:    ",
           namesof("p0", lp0, earg= ep0, tag=FALSE), ", ",
           namesof("munb", lmunb, earg= emunb, tag=FALSE), ", ",
           namesof("k", lk, earg= ek, tag=FALSE), "\n",
           "Mean:     (1-p0) * munb / [1 - (k/(k+munb))^k]"),
    constraints=eval(substitute(expression({
        temp752 = .zero
        if(length(temp752) && all(temp752 == -3))
            temp752 = 3*(1:ncol(y))
        constraints = cm.zero.vgam(constraints, x, temp752, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        y = as.matrix(y)
        extra$NOS = NOS = ncoly = ncol(y)  # Number of species
        M = 3 * ncoly # 

        mynames1 = if(NOS==1) "p0" else paste("p0", 1:NOS, sep="")
        mynames2 = if(NOS==1) "munb" else paste("munb", 1:NOS, sep="")
        mynames3 = if(NOS==1) "k" else paste("k", 1:NOS, sep="")
        predictors.names =
            c(namesof(mynames1, .lp0, earg= .ep0, tag= FALSE),
              namesof(mynames2, .lmunb, earg= .emunb, tag= FALSE),
              namesof(mynames3, .lk, earg= .ek, tag= FALSE))
        predictors.names = predictors.names[interleave.VGAM(3*NOS, M=3)]
        extra$y0 = y0 = ifelse(y==0, 1, 0)
        extra$ymat = ymat = cbind(y0=y0, y=y)
        extra$skip.these = skip.these = matrix(as.logical(y0), n, NOS)

        if(!length(etastart)) {
            if( .method.init == 3) {
                mu.init = y + 1/16
            } else {
                mu.init = y
                for(iii in 1:ncol(y))
                    mu.init[,iii] = if( .method.init == 2)
                        weighted.mean(y[,iii], w=w) else
                        median(rep(y[,iii], w)) + 1/8
            }
            kmat0 = matrix( .ik, nrow(y), ncoly, byrow=TRUE) # Initial kmat
            pnb0 = (kmat0 / (kmat0 + mu.init))^kmat0
            etastart = cbind(theta2eta((0.5 + w * y0) / (1 + w), .lp0, earg= .ep0 ),
                             theta2eta(mu.init*(1-pnb0), .lmunb, earg= .emunb ),
                             theta2eta(kmat0, .lk, earg= .ek ))
            etastart = etastart[,interleave.VGAM(ncol(etastart),M=3)]
        }
    }), list( .lp0=lp0, .lmunb=lmunb, .lk=lk, .ik=ik,
              .ep0=ep0, .emunb=emunb, .ek=ek,
              .method.init=method.init ))), 
    inverse=eval(substitute(function(eta, extra=NULL) {
        NOS = extra$NOS
        p0 = eta2theta(eta[,3*(1:NOS)-2], .lp0, earg= .ep0 )
        munb = eta2theta(eta[,3*(1:NOS)-1,drop=FALSE], .lmunb, earg= .emunb )
        kmat = eta2theta(eta[,3*(1:NOS),drop=FALSE], .lk, earg= .ek )
        pnb0 = (kmat / (kmat + munb))^kmat # p(0) from negative binomial
        (1 - p0) * munb / (1 - pnb0)
    }, list( .lp0=lp0, .lk=lk, .lmunb=lmunb,
             .ep0=ep0, .emunb=emunb, .ek=ek ))),
    last=eval(substitute(expression({
        misc$link = c(rep( .lp0, length=NOS), rep( .lmunb, length=NOS),
                      rep( .lk, length=NOS))
        temp.names = c(mynames1, mynames2, mynames3)
        temp.names = temp.names[interleave.VGAM(3*NOS, M=3)]
        names(misc$link) = temp.names
        misc$earg = vector("list", 3*NOS)
        names(misc$earg) = temp.names
        for(ii in 1:NOS) {
            misc$earg[[3*ii-2]] = .ep0
            misc$earg[[3*ii-1]] = .emunb
            misc$earg[[3*ii  ]] = .ek
        }
        misc$cutoff = .cutoff
        misc$method.init = .method.init
    }), list( .lp0=lp0, .lmunb=lmunb, .lk=lk, .cutoff=cutoff,
              .ep0=ep0, .emunb=emunb, .ek=ek,
              .method.init=method.init ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals=FALSE, eta,extra=NULL) {
        NOS = extra$NOS
        p0 = eta2theta(eta[,3*(1:NOS)-2,drop=FALSE], .lp0, earg= .ep0 )
        munb = eta2theta(eta[,3*(1:NOS)-1,drop=FALSE], .lmunb, earg= .emunb )
        kmat = eta2theta(eta[,3*(1:NOS),drop=FALSE], .lk, earg= .ek )
        skip = extra$skip.these
        pnb0 = (kmat / (kmat + munb))^kmat
        ans = 0.0
        for(spp. in 1:NOS) {
            i8 = skip[,spp.]
            ans = ans + sum(w[i8] * log(p0[i8,spp.])) +
            sum(w[!i8] * (log(1-p0[!i8,spp.]) + y[!i8,spp.] * 
                log(munb[!i8,spp.]/(munb[!i8,spp.]+
                kmat[!i8,spp.])) + kmat[!i8,spp.]*log(kmat[!i8,spp.] /
                (munb[!i8,spp.]+kmat[!i8,spp.])) +
                lgamma(y[!i8,spp.]+kmat[!i8,spp.]) - 
                lgamma(kmat[!i8,spp.]) - lgamma(y[!i8,spp.]+1) -
                (if(is.R())
                log1p(-pnb0[!i8,spp.]) else log(1 - pnb0[!i8,spp.]))))
        }
        ans
    }, list( .lp0=lp0, .lmunb=lmunb, .lk=lk,
             .ep0=ep0, .emunb=emunb, .ek=ek ))),
    vfamily=c("zanegbinomial"),
    deriv=eval(substitute(expression({
        NOS = extra$NOS
        y0 = extra$y0
        p0 = eta2theta(eta[,3*(1:NOS)-2], .lp0, earg= .ep0 )
        munb = eta2theta(eta[,3*(1:NOS)-1,drop=FALSE], .lmunb, earg= .emunb )
        kmat = eta2theta(eta[,3*(1:NOS),drop=FALSE], .lk, earg= .ek )
        skip = extra$skip.these

        d3 = deriv3(~ -log(1 - (kmat. /(kmat. + munb. ))^kmat. ),
                    c("munb.", "kmat."), hessian= TRUE) # Extra term
        dl0.dthetas =  array(NA, c(n, NOS, 2))
        d2l0.dthetas2 =  array(NA, c(n, NOS, 3))  # matrix-band format
        for(spp. in 1:NOS) {
            kmat. = kmat[,spp.]
            munb. = munb[,spp.]
            eval.d3 = eval(d3)  # Evaluated for one species
            dl0.dthetas[,spp.,1] =  attr(eval.d3, "gradient")[,1]
            dl0.dthetas[,spp.,2] =  attr(eval.d3, "gradient")[,2]
            d2l0.dthetas2[,spp.,1] =  attr(eval.d3, "hessian")[,1,1]
            d2l0.dthetas2[,spp.,2] =  attr(eval.d3, "hessian")[,2,2]
            d2l0.dthetas2[,spp.,3] =  attr(eval.d3, "hessian")[,1,2]
        }
        dl.dmunb = y/munb - (y+kmat)/(kmat+munb) + dl0.dthetas[,,1]  
        dl.dk = digamma(y+kmat) - digamma(kmat) - (y+kmat)/(munb+kmat) + 1 +
                log(kmat/(kmat+munb)) + dl0.dthetas[,,2]  
        for(spp. in 1:NOS)
            dl.dk[skip[,spp.],spp.] = dl.dmunb[skip[,spp.],spp.] = 0

        dmunb.deta = dtheta.deta(munb, .lmunb, earg= .emunb )
        dk.deta = dtheta.deta(kmat, .lk, earg= .ek )
        myderiv = w * cbind(dl.dmunb * dmunb.deta, dl.dk * dk.deta)

        mup0 = p0
        temp3 = if(.lp0 == "logit") {
            w * (y0 - mup0)
        } else
            w * dtheta.deta(mup0, link=.lp0, earg= .ep0 ) * (y0/mup0 - 1) / (1-mup0)

        ans = cbind(temp3, myderiv)
        ans = ans[,interleave.VGAM(ncol(ans), M=3)]
        ans
    }), list( .lp0=lp0, .lmunb=lmunb, .lk=lk,
              .ep0=ep0, .emunb=emunb, .ek=ek ))),
    weight=eval(substitute(expression({
        wz = matrix(0, n, 6*NOS-1)  # wz is not 'diagonal' 
        pnb0 = (kmat / (kmat + munb))^kmat
        ed2l.dmunb2 = (1/munb - (munb + kmat*(1-pnb0))/(munb+kmat)^2) /
                      (1-pnb0) - d2l0.dthetas2[,,1]
        fred = dotFortran(name="enbin8",
                      ans=double(n*NOS), as.double(kmat),
                      as.double(kmat/(munb+kmat)), as.double(.cutoff),
                      as.integer(n), ok=as.integer(1), as.integer(NOS),
                      sumpdf=double(1), macheps=as.double(.Machine$double.eps))
        if(fred$ok != 1)
            stop("error in Fortran subroutine exnbin")
        dim(fred$ans) = c(n, NOS)
        ed2l.dk2 = -fred$ans/(1-pnb0) - 1/kmat + 1/(kmat+munb) -
                   munb * pnb0 / ((1-pnb0)*(munb+kmat)^2) - d2l0.dthetas2[,,2]
        wz[,3*(1:NOS)-1] = w * dmunb.deta^2 * ed2l.dmunb2
        wz[,3*(1:NOS)] = w * dk.deta^2 * ed2l.dk2

        wz[,3*NOS+3*(1:NOS)-1] = -w * d2l0.dthetas2[,,3] * dmunb.deta * dk.deta

        tmp100 = mup0*(1-mup0)
        tmp200 = if(.lp0 == "logit") {
            cbind(w * tmp100)
        } else {
            cbind(w * dtheta.deta(mup0, link= .lp0, earg= .ep0 )^2 / tmp100)
        }
        for(ii in 1:NOS) {
            index200 = abs(tmp200[,ii]) < .Machine$double.eps
            if(any(index200)) {
                tmp200[index200,ii] = .Machine$double.eps # Diagonal 0's are bad 
            }
        }
        wz[,3*(1:NOS)-2] =  tmp200

        for(spp. in 1:NOS) {
            wz[skip[,spp.],3*spp. - 1] = 
            wz[skip[,spp.],3*spp.] = .Machine$double.eps^0.5
            wz[skip[,spp.],3*NOS+3*(spp.)-1] = 0
        }
        wz
    }), list( .lp0=lp0, .ep0=ep0, .cutoff=cutoff ))))
}


rposnegbin = function(n, munb, k) {
    if(!is.Numeric(k, posit=TRUE))
        stop("argument \"k\" must be positive")
    if(!is.Numeric(munb, posit=TRUE))
        stop("argument \"munb\" must be positive")
    if(!is.Numeric(n, posit=TRUE, integ=TRUE, allow=1))
        stop("argument \"n\" must be a positive integer")
    ans = rnbinom(n=n, mu=munb, size=k)
    munb = rep(munb, len=n)
    k = rep(k, len=n)
    index = ans == 0
    while(any(index)) {
        more = rnbinom(n=sum(index), mu=munb[index], size=k[index])
        ans[index] = more
        index = ans == 0
    }
    ans
} 



zipoisson = function(lphi="logit", llambda="loge",
                     ephi=list(), elambda =list(),
                     iphi=NULL, zero=NULL)
{
    if(mode(lphi) != "character" && mode(lphi) != "name")
        lphi = as.character(substitute(lphi))
    if(mode(llambda) != "character" && mode(llambda) != "name")
        llambda = as.character(substitute(llambda))
    if(is.Numeric(iphi))
        if(!is.Numeric(iphi, allow=1, posit=TRUE) || iphi >= 1)
            stop("iphi must be a single number inside the interval (0,1)")
    if(!is.list(ephi)) ephi = list()
    if(!is.list(elambda)) elambda = list()

    new("vglmff",
    blurb=c("Zero-inflated Poisson\n\n",
            "Links:    ", namesof("phi", lphi, earg= ephi), ", ",
            namesof("lambda", llambda, earg= elambda), "\n",
            "Mean:     (1-phi)*lambda"),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(ncol(as.matrix(y)) != 1) stop("multivariate responses not allowed")
        predictors.names = c(
            namesof("phi", .lphi, earg= .ephi, tag=FALSE),
            namesof("lambda", .llambda, earg= .ephi, tag=FALSE))
        if(!length(etastart)) {
            phi.init = if(length( .iphi)) .iphi else {
                sum(w[y==0]) / sum(w)
            }
            if(phi.init <= 0 || phi.init >=1) phi.init = 0.1  # Last resort
            lambda.init = y + 1/8
            etastart = cbind(theta2eta(rep(phi.init, len=n), .lphi, earg= .ephi ),
                             theta2eta(lambda.init, .llambda, earg= .ephi ))
        }
    }), list( .lphi=lphi, .llambda=llambda,
              .ephi=ephi, .elambda=elambda,
              .iphi=iphi ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        phi = eta2theta(eta[,1], .lphi, earg= .ephi )
        lambda = eta2theta(eta[,2], .llambda, earg= .elambda )
        (1-phi) * lambda
    }, list( .lphi=lphi, .llambda=llambda,
             .ephi=ephi, .elambda=elambda ))),
    last=eval(substitute(expression({
        misc$link <- c("phi" = .lphi, "lambda" = .llambda)
        misc$earg <- list("phi" = .ephi, "lambda" = .elambda)
        if(intercept.only) {
            phi = eta2theta(eta[1,1], .lphi, earg= .ephi )
            lambda = eta2theta(eta[1,2], .llambda, earg= .elambda )
            misc$prob0 = phi + (1-phi) * exp(-lambda) # P(Y=0)
        }
    }), list( .lphi=lphi, .llambda=llambda,
              .ephi=ephi, .elambda=elambda ))),
    loglikelihood=eval(substitute( 
        function(mu,y,w,residuals=FALSE, eta, extra=NULL) {
        phi = eta2theta(eta[,1], .lphi, earg= .ephi )
        lambda = eta2theta(eta[,2], .llambda, earg= .elambda )
        index = (y==0)
        tmp8 = phi + (1-phi)*exp(-lambda)
        ell0 = log(tmp8[index])
        ell1 = log((1-phi[!index]) * dpois(y[!index], lambda= lambda[!index]))
        sum(w[index] * ell0) + sum(w[!index] * ell1)
    }, list( .lphi=lphi, .llambda=llambda,
             .ephi=ephi, .elambda=elambda ))),
    vfamily=c("zipoisson"),
    deriv=eval(substitute(expression({
        phi = eta2theta(eta[,1], .lphi, earg= .ephi )
        lambda = eta2theta(eta[,2], .llambda, earg= .elambda )
        tmp8 = phi + (1-phi)*exp(-lambda)
        index = (y==0)
        dl.dphi = (1-exp(-lambda)) / tmp8
        dl.dphi[!index] = -1 / (1-phi[!index])
        dl.dlambda = -(1-phi) * exp(-lambda) / tmp8
        dl.dlambda[!index] = (y[!index] - lambda[!index]) / lambda[!index] 
        dphi.deta = dtheta.deta(phi, .lphi, earg= .ephi)
        dlambda.deta = dtheta.deta(lambda, .llambda, earg= .elambda )
        ans = w * cbind(dl.dphi * dphi.deta, dl.dlambda * dlambda.deta)
        if(.llambda == "loge" && (any(lambda[!index] < .Machine$double.eps))) {
            ans[!index,2] = w[!index] * (y[!index] - lambda[!index])
        }
        ans
    }), list( .lphi=lphi, .llambda=llambda,
              .ephi=ephi, .elambda=elambda ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), nrow=n, ncol=dimm(M))
        d2l.dphi2 = (1-exp(-lambda)) / ((1-phi)*tmp8)
        d2l.dlambda2 = (1-phi)/lambda - phi*(1-phi)*exp(-lambda) / tmp8
        d2l.dphilambda = -exp(-lambda) / tmp8
        wz[,iam(1,1,M)] = d2l.dphi2 * dphi.deta^2
        wz[,iam(2,2,M)] = d2l.dlambda2 * dlambda.deta^2
        wz[,iam(1,2,M)] = d2l.dphilambda * dphi.deta * dlambda.deta
        if(.llambda == "loge" && (any(lambda[!index] < .Machine$double.eps))) {
            ind5 = !index & (lambda < .Machine$double.eps)
            if(any(ind5))
                wz[ind5,iam(2,2,M)] = (1-phi[ind5]) * .Machine$double.eps
        }
        w * wz
    }), list( .lphi=lphi, .llambda=llambda,
              .ephi=ephi, .elambda=elambda ))))
}




zibinomial = function(lphi="logit", link.mu="logit",
                      ephi=list(), emu=list(),
                      iphi=NULL, zero=1, mv=FALSE)
{
    if(as.logical(mv)) stop("argument \"mv\" must be FALSE")
    if(mode(lphi) != "character" && mode(lphi) != "name")
        lphi = as.character(substitute(lphi))
    if(mode(link.mu) != "character" && mode(link.mu) != "name")
        link.mu = as.character(substitute(link.mu))
    if(is.Numeric(iphi))
        if(!is.Numeric(iphi, allow=1, posit=TRUE) || iphi >= 1)
            stop("iphi must be a single number inside the interval (0,1)")
    if(!is.list(ephi)) ephi = list()
    if(!is.list(emu)) emu = list()

    new("vglmff",
    blurb=c("Zero-inflated binomial\n\n",
            "Links:    ", namesof("phi", lphi, earg= ephi ), ", ",
            namesof("mu", link.mu, earg= emu ), "\n",
            "Mean:     (1-phi) * mu / (1 - (1-mu)^w)"),
    constraints=eval(substitute(expression({
        constraints <- cm.zero.vgam(constraints, x, .zero, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        {
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
                stop("Response not of the right form (1 or 2 columns required)")
        }

        predictors.names = c( namesof("phi", .lphi, earg= .ephi, tag=FALSE),
                              namesof("mu",  .link.mu, earg= .emu, tag=FALSE))
        phi.init = if(length( .iphi)) .iphi else {
            sum(w[y==0]) / sum(w)
        }
        if(phi.init <= 0 || phi.init >=1) phi.init = 0.1  # Last resort
        mustart = cbind(rep(phi.init, len=n), mustart) # 1st coln not a real mu
    }), list( .lphi=lphi, .link.mu=link.mu,
              .ephi=ephi, .emu=emu,
              .iphi=iphi ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        phi = eta2theta(eta[,1], .lphi, earg= .ephi )
        mubin = eta2theta(eta[,2], .link.mu, earg= .emu )
        (1-phi) * mubin
    }, list( .lphi=lphi, .link.mu=link.mu,
             .ephi=ephi, .emu=emu ))),
    last=eval(substitute(expression({
        misc$link <- c("phi" = .lphi, "mu" = .link.mu)
        misc$earg <- list("phi" = .ephi, "mu" = .emu )
        if(intercept.only && all(w==w[1])) {
            phi = eta2theta(eta[1,1], .lphi, earg= .ephi )
            mubin = eta2theta(eta[1,2], .link.mu, earg= .emu )
            misc$p0 = phi + (1-phi) * (1-mubin)^w[1] # P(Y=0)
        }
    }), list( .lphi=lphi, .link.mu=link.mu,
              .ephi=ephi, .emu=emu ))),
    link=eval(substitute(function(mu, extra=NULL)
        cbind(theta2eta(mu[,1], .lphi, earg= .ephi ),
              theta2eta(mu[,2], .link.mu, earg= .emu ))
    , list( .lphi=lphi, .link.mu=link.mu,
            .ephi=ephi, .emu=emu ))),
    loglikelihood=eval(substitute( 
        function(mu,y,w,residuals=FALSE, eta, extra=NULL) {
        phi = eta2theta(eta[,1], .lphi, earg= .ephi )
        mubin = eta2theta(eta[,2], .link.mu, earg= .emu )
        index = (y==0)
        tmp8 = phi + (1-phi)*(1-mubin)^w
        ell0 = log(tmp8[index])
        ell1 = log(1-phi[!index]) + dbinom(x=round(w[!index]*y[!index]), 
               size=w[!index], prob=mubin[!index], log=TRUE)
        sum(ell0) + sum(ell1)
    }, list( .lphi=lphi, .link.mu=link.mu,
             .ephi=ephi, .emu=emu ))),
    vfamily=c("zibinomial"),
    deriv=eval(substitute(expression({
        phi = eta2theta(eta[,1], .lphi, earg= .ephi )
        mubin = eta2theta(eta[,2], .link.mu, earg= .emu )
        prob0 = (1-mubin)^w    # Actually q^w
        tmp8 = phi + (1-phi)*prob0
        index = (y==0)
        dl.dphi = (1-prob0) / tmp8
        dl.dphi[!index] = -1 / (1-phi[!index])
        dl.dmubin = -w * (1-phi) * (1-mubin)^(w-1) / tmp8
        dl.dmubin[!index] = w[!index] * (y[!index]/mubin[!index] - 
            (1-y[!index]) / (1-mubin[!index]))
        dphi.deta = dtheta.deta(phi, .lphi, earg= .ephi )
        dmubin.deta = dtheta.deta(mubin, .link.mu, earg= .emu )
        ans = cbind(dl.dphi * dphi.deta, dl.dmubin * dmubin.deta)
        if(.link.mu == "logit") {
            ans[!index,2] = w[!index] * (y[!index] - mubin[!index])
        }
        ans
    }), list( .lphi=lphi, .link.mu=link.mu,
              .ephi=ephi, .emu=emu ))),
    weight=eval(substitute(expression({
        wz = matrix(as.numeric(NA), nrow=n, ncol=dimm(M))
        d2l.dphi2 = (1-mubin^w) / ((1-phi) * tmp8)
        d2l.dmubin2 = w * (1-phi) * ((1 - mubin * (1-mubin)^(w-1)) / 
                      (mubin*(1-mubin)) - mubin^(w-2) * (w*phi-tmp8) / tmp8)
        d2l.dphimubin = -w * (1-mubin)^(w-1) / tmp8
        wz[,iam(1,1,M)] = d2l.dphi2 * dphi.deta^2
        wz[,iam(2,2,M)] = d2l.dmubin2 * dmubin.deta^2
        wz[,iam(1,2,M)] = d2l.dphimubin * dphi.deta * dmubin.deta
        if(TRUE) {
            ind6 = wz[,iam(2,2,M)] < .Machine$double.eps
            if(any(ind6))
                wz[ind6,iam(2,2,M)] = .Machine$double.eps
        }
        wz
    }), list( .lphi=lphi, .link.mu=link.mu,
              .ephi=ephi, .emu=emu ))))
}



dzibinom = function(x, size, prob, log = FALSE, phi=0) {
    L = max(length(x), length(size), length(prob), length(phi))
    x = rep(x, len=L); size = rep(size, len=L);
    prob = rep(prob, len=L); phi = rep(phi, len=L);
    ans = dbinom(x, size, prob, log=log)
    if(!is.Numeric(phi) || any(phi < 0) || any(phi > 1))
        stop("phi must be between 0 and 1 inclusive")
    ifelse(x==0, phi + (1-phi) * ans, (1-phi) * ans)
}

pzibinom = function(q, size, prob, lower.tail = TRUE, log.p = FALSE, phi=0) {
    ans = pbinom(q, size, prob, lower.tail = lower.tail, log.p = log.p)
    phi = rep(phi, length=length(ans))
    if(!is.Numeric(phi) || any(phi < 0) || any(phi > 1))
        stop("phi must be between 0 and 1 inclusive")
    phi + (1-phi) * ans
}

qzibinom = function(p, size, prob, lower.tail = TRUE, log.p = FALSE, phi=0) {
    nn = max(length(p), length(size), length(prob), length(phi))
    p = rep(p, len=nn)
    size = rep(size, len=nn)
    prob = rep(prob, len=nn)
    phi = rep(phi, len=nn)
    if(!is.Numeric(phi) || any(phi < 0) || any(phi > 1))
        stop("phi must be between 0 and 1 inclusive")
    ans = p 
    ans[p<=phi] = 0 
    ans[p>phi] = qbinom((p[p>phi]-phi[p>phi])/(1-phi[p>phi]), size[p>phi],
                        prob[p>phi], lower.tail = lower.tail, log.p = log.p)
    ans
}

rzibinom = function(n, size, prob, phi=0) {
    if(!is.Numeric(n, positive=TRUE, integer=TRUE, allow=1))
        stop("n must be a single positive integer")
    ans = rbinom(n, size, prob)
    phi = rep(phi, len=length(ans))
    if(!is.Numeric(phi) || any(phi < 0) || any(phi > 1))
        stop("phi must be between 0 and 1 inclusive")
    ifelse(runif(n) < phi, 0, ans)
}







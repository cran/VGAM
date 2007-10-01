# These functions are
# Copyright (C) 1998-2007 T.W. Yee, University of Auckland. All rights reserved.





posnegbinomial = function(lmunb = "loge", lk = "loge",
                          emunb =list(), ek = list(),
                          ik = NULL, zero = -2, cutoff = 0.995,
                          method.init=1)
{
    if(!is.Numeric(cutoff, allow=1) || cutoff<0.8 || cutoff>=1)
        stop("range error in the argument \"cutoff\"")
    if(!is.Numeric(method.init, allow=1, integ=TRUE, posit=TRUE) ||
       method.init > 3) stop("argument \"method.init\" must be 1, 2 or 3")
    if(length(ik) && !is.Numeric(ik, posit=TRUE))
        stop("bad input for argument \"ik\"")

    if(mode(lmunb) != "character" && mode(lmunb) != "name")
        lmunb = as.character(substitute(lmunb))
    if(mode(lk) != "character" && mode(lk) != "name")
        lk = as.character(substitute(lk))
    if(!is.list(emunb)) emunb = list()
    if(!is.list(ek)) ek = list()

    new("vglmff",
    blurb=c("Positive-negative binomial distribution\n\n",
           "Links:    ",
           namesof("munb", lmunb, earg= emunb ), ", ",
           namesof("k", lk, earg= ek ), "\n",
           "Mean:     munb / (1 - (k/(k+munb))^k)"),
    constraints=eval(substitute(expression({
        temp752 = .zero
        if(length(temp752) && all(temp752 == -2))
            temp752 = 2*(1:ncol(y))
        constraints = cm.zero.vgam(constraints, x, temp752, M)
    }), list( .zero=zero ))),
    initialize=eval(substitute(expression({
        if(any(y==0)) stop("there are zero values in the response")
        y = as.matrix(y) 
        M = 2 * ncol(y) 
        extra$NOS = NOS = ncoly = ncol(y)  # Number of species
        predictors.names = c(namesof(if(NOS==1) "munb" else
            paste("munb", 1:NOS, sep=""), .lmunb, earg= .emunb, tag= FALSE),
            namesof(if(NOS==1) "k" else paste("k", 1:NOS, sep=""),
            .lk, earg= .ek, tag= FALSE))
        predictors.names = predictors.names[interleave.VGAM(M, M=2)]
        if(!length(etastart)) {
            if( .method.init == 3) {
                mu.init = y + 1/8
            } else {
                mu.init = y
                for(iii in 1:ncol(y))
                    mu.init[,iii] = if( .method.init == 2)
                        weighted.mean(y[,iii], w=w) else
                        median(rep(y[,iii], w)) + 1/8
            }
            if( is.Numeric( .ik )) {
                kmat0 = matrix( .ik, nr=n, nc=NOS, byrow=TRUE)
            } else {
                kmat0 = matrix(0, nr=n, nc=NOS)
                Loglikfun = function(y, munb, kmat, w) {
                    p0 = (kmat / (kmat + munb))^kmat
        sum(w * (y * log(munb/(munb+kmat)) + kmat*log(kmat/(munb+kmat)) +
                 lgamma(y+kmat) - lgamma(kmat) - lgamma(y+1) - 
                 (if(is.R()) log1p(-p0) else log1p(-p0)))) }

                k.grid = rvar = 2^((-3):6)
                for(spp. in 1:NOS) {
                    for(ii in 1:length(k.grid))
                        rvar[ii] = Loglikfun(y=y[,spp.], mu=mu.init[,spp.],
                                             kmat=k.grid[ii], w=w)
                    try.this = k.grid[rvar == max(rvar)]
                    kmat0[,spp.] = try.this
                }
            }
            p00 = (kmat0 / (kmat0 + mu.init))^kmat0
            etastart = cbind(theta2eta(mu.init*(1-p00), .lmunb, earg= .emunb ),
                             theta2eta(kmat0, .lk, earg= .ek ))
            etastart = etastart[,interleave.VGAM(M, M=2),drop=FALSE]
        }
    }), list( .lmunb=lmunb, .lk=lk, .ik=ik,
              .emunb=emunb, .ek=ek,
              .method.init=method.init ))),
    inverse=eval(substitute(function(eta, extra=NULL) {
        NOS = ncol(eta) / 2
        munb = eta2theta(eta[,2*(1:NOS)-1,drop=FALSE], .lmunb, earg= .emunb )
        kmat = eta2theta(eta[,2*(1:NOS),drop=FALSE], .lk, earg= .ek )
        p0 = (kmat / (kmat + munb))^kmat
        munb / (1 - p0)
    }, list( .lk=lk, .lmunb=lmunb,
             .emunb=emunb, .ek=ek ))),
    last=eval(substitute(expression({
        temp0303 = c(rep( .lmunb, length=NOS), rep( .lk, length=NOS))
        names(temp0303) = c(if(NOS==1) "munb" else paste("munb", 1:NOS, sep=""), 
                            if(NOS==1) "k" else paste("k", 1:NOS, sep=""))
        temp0303 = temp0303[interleave.VGAM(M, M=2)]
        misc$link = temp0303  # Already named
        misc$earg = vector("list", 2*NOS)
        names(misc$earg) = names(misc$link)
        for(ii in 1:NOS) {
            misc$earg[[2*ii-1]] = .emunb
            misc$earg[[2*ii  ]] = .ek
        }
        misc$cutoff = .cutoff 
        misc$method.init = .method.init
    }), list( .lmunb=lmunb, .lk=lk, .cutoff=cutoff,
              .emunb=emunb, .ek=ek,
              .method.init=method.init ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals= FALSE,eta, extra=NULL) {
        NOS = ncol(eta) / 2
        munb = eta2theta(eta[,2*(1:NOS)-1,drop=FALSE], .lmunb, earg= .emunb )
        kmat = eta2theta(eta[,2*(1:NOS),drop=FALSE], .lk, earg= .ek )
        p0 = (kmat / (kmat + munb))^kmat
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (y * log(munb/(munb+kmat)) + kmat*log(kmat/(munb+kmat)) +
                 lgamma(y+kmat) - lgamma(kmat) - lgamma(y+1) -
                 (if(is.R()) log1p(-p0) else log1p(-p0))))
    }, list( .lmunb=lmunb, .lk=lk,
             .emunb=emunb, .ek=ek ))),
    vfamily=c("posnegbinomial"),
    deriv=eval(substitute(expression({
        NOS= extra$NOS
        munb = eta2theta(eta[,2*(1:NOS)-1,drop=FALSE], .lmunb, earg= .emunb )
        kmat = eta2theta(eta[,2*(1:NOS),drop=FALSE], .lk, earg= .ek )
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

        NOS = ncol(eta) / 2
        dl.dmunb = y/munb - (y+kmat)/(kmat+munb) + dl0.dthetas[,,1]
        dl.dk = digamma(y+kmat) - digamma(kmat) - (y+kmat)/(munb+kmat) + 1 +
                log(kmat/(kmat+munb)) + dl0.dthetas[,,2]
        dmunb.deta = dtheta.deta(munb, .lmunb, earg= .emunb )
        dk.deta = dtheta.deta(kmat, .lk, earg= .ek )
        myderiv = w * cbind(dl.dmunb * dmunb.deta, dl.dk * dk.deta)
        myderiv[,interleave.VGAM(M, M=2)]
    }), list( .lmunb=lmunb, .lk=lk,
              .emunb=emunb, .ek=ek ))),
    weight=eval(substitute(expression({
        wz = matrix(0, n, 4*NOS-1)  # wz is no longer 'diagonal' 
        p0 = (kmat / (kmat + munb))^kmat
        ed2l.dmunb2 = (1/munb - (munb + kmat*(1-p0))/(munb+kmat)^2) / (1-p0) -
                      d2l0.dthetas2[,,1]
        fred = dotFortran(name="enbin8",
                      ans=double(n*NOS),
                      as.double(kmat),
                      as.double(kmat/(munb+kmat)),
                      as.double(.cutoff),
                      as.integer(n), ok=as.integer(1), as.integer(NOS),
                      sumpdf=double(1), macheps=as.double(.Machine$double.eps))
        if(fred$ok != 1)
            stop("error in Fortran subroutine exnbin")
        dim(fred$ans) = c(n, NOS)
        ed2l.dk2 = -fred$ans/(1-p0) - 1/kmat + 1/(kmat+munb) -
                   munb * p0 / ((1-p0)*(munb+kmat)^2) - d2l0.dthetas2[,,2]
        wz[,2*(1:NOS)-1] = dmunb.deta^2 * ed2l.dmunb2
        wz[,2*(1:NOS)] = dk.deta^2 * ed2l.dk2
        wz[,2*NOS+2*(1:NOS)-1] = -d2l0.dthetas2[,,3] * dmunb.deta * dk.deta
        w * wz
    }), list( .cutoff=cutoff ))))
}



dpospois = function(x, lambda) {
    if(!is.Numeric(lambda, posit=TRUE))
        stop("bad input for argument \"lambda\"")
    L = max(length(x), length(lambda))
    x = rep(x, len=L); lambda = rep(lambda, len=L); 
    ans = ifelse(x==0, 0, dpois(x, lambda) / (1 - exp(-lambda)))
    ans
}


ppospois = function(q, lambda) {
    if(!is.Numeric(lambda, posit=TRUE))
        stop("bad input for argument \"lambda\"")
    L = max(length(q), length(lambda))
    q = rep(q, len=L); lambda = rep(lambda, len=L); 
    ifelse(q<1, 0, (ppois(q, lambda) - exp(-lambda)) / (1 - exp(-lambda)))
}

qpospois = function(p, lambda) {
    if(!is.Numeric(lambda, posit=TRUE))
        stop("bad input for argument \"lambda\"")
    if(!is.Numeric(p, posit=TRUE) || any(p >= 1))
        stop("bad input for argument \"p\"")
    qpois(p * (1 - exp(-lambda)) + exp(-lambda), lambda)
}


rpospois = function(n, lambda) {
    if(!is.Numeric(n, posit=TRUE, allow=1, integ=TRUE))
        stop("bad input for argument \"n\"")
    if(!is.Numeric(lambda, posit=TRUE))
        stop("bad input for argument \"lambda\"")
    ans = rpois(n, lambda)
    lambda = rep(lambda, len=n)
    index = ans == 0
    while(any(index)) {
        more = rpois(sum(index), lambda[index])
        ans[index] = more
        index = ans == 0
    }
    ans
}




pospoisson = function(link="loge", earg=list())
{
    if(!missing(link))
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()

    new("vglmff",
    blurb=c("Positive-Poisson distribution\n\n",
           "Links:    ",
           namesof("lambda", link, earg= earg, tag=FALSE),
           "\n"),
    initialize=eval(substitute(expression({
        y = as.matrix(y)
        predictors.names = namesof(if(ncol(y)==1) "lambda"
            else paste("lambda", 1:ncol(y), sep=""), .link,
            earg= .earg, tag=FALSE)
        if(!length(etastart))
            etastart = theta2eta(y / (1-exp(-y)), .link, earg= .earg )
    }), list( .link=link, .earg= earg ))), 
    inverse=eval(substitute(function(eta, extra=NULL) {
        lambda = eta2theta(eta, .link, earg= .earg )
        lambda / (1-exp(-lambda))
    }, list( .link=link, .earg= earg ))),
    last=eval(substitute(expression({
        misc$link = rep( .link, len=M)
        names(misc$link) = if(M==1) "lambda" else paste("lambda", 1:M, sep="")
        misc$earg = vector("list", M)
        names(misc$earg) = names(misc$link)
        for(ii in 1:M)
            misc$earg[[ii]] = .earg
    }), list( .link=link, .earg= earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals=FALSE, eta,extra=NULL) {
        lambda = eta2theta(eta, .link, earg= .earg ) 
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(w * (-log1p(-exp(-lambda)) - lambda + y*log(lambda)))
    }, list( .link=link, .earg= earg ))),
    vfamily=c("pospoisson"),
    deriv=eval(substitute(expression({
         lambda = eta2theta(eta, .link, earg= .earg ) 
         dl.dlambda = y/lambda - 1 - 1/(exp(lambda)-1)
         dlambda.deta = dtheta.deta(lambda, .link, earg= .earg )
         w * dl.dlambda * dlambda.deta
    }), list( .link=link, .earg= earg ))),
    weight=eval(substitute(expression({
         temp = exp(lambda)
         ed2l.dlambda2 = -temp * (1/lambda - 1/(temp-1)) / (temp-1)
         wz = -w * (dlambda.deta^2) * ed2l.dlambda2
         wz
    }), list( .link=link, .earg= earg ))))
}



posbinomial = function(link="logit", earg=list())
{

    if(!missing(link))
        link = as.character(substitute(link))
    if(!is.list(earg)) earg = list()
       
    new("vglmff",
    blurb=c("Positive-Binomial distribution\n\n",
           "Links:    ",
           namesof("p", link, earg= earg, tag=FALSE), "\n"),
    initialize=eval(substitute(expression({
      	eval(binomialff(link= .link)@initialize)
        predictors.names = namesof("p", .link, earg= .earg , tag=FALSE)
	if(length(extra)) extra$w = w else extra = list(w=w)
        if(!length(etastart))
	    etastart = cbind(theta2eta(mustart, .link, earg= .earg ))
    }), list( .link = link, .earg=earg ))),
    inverse=eval(substitute(function(eta, extra=NULL){
        theta = eta2theta(eta, .link, earg= .earg )
        theta/(1-(1-theta)^(extra$w))},
    list(.link=link, .earg=earg ))),
    last=eval(substitute(expression({
        extra$w = NULL   # Kill it off 
        misc$link = c(p = .link)
        misc$earg = list(p = .earg )
    }), list( .link=link, .earg=earg ))),
    loglikelihood=eval(substitute(
        function(mu,y,w,residuals=FALSE,eta,extra=NULL) {
        yi = round(y*w)
        theta = eta2theta(eta, .link, earg= .earg )
        if(residuals) stop("loglikelihood residuals not implemented yet") else
        sum(yi*log(theta)+(w-yi)*log1p(-theta)-log1p(-(1-theta)^w))
    }, list( .link=link, .earg=earg ))),
    vfamily=c("posbinomial"),
    deriv=eval(substitute(expression({
        yi = round(y*w)     
        theta = eta2theta(eta, .link, earg= .earg )
        dldtheta = yi/theta-(w-yi)/(1-theta)-w*(1-theta)^(w-1) /
                    (1-(1-theta)^w)
        dthetadeta = dtheta.deta(theta, .link, earg= .earg )
        dldtheta * dthetadeta
    }), list( .link=link, .earg=earg ))),
    weight=eval(substitute(expression({
        temp = 1 - (1-theta)^w
        temp2 = (1-theta)^2
        ed2ldtheta2 = -w/(theta*temp) - w/temp2 + w*theta/(temp2*temp) +
            w*(w-1)* (1-theta)^(w-2) /temp +
            w^2 * temp2^(w-1) / temp^2
        wz = -(dthetadeta^2) * ed2ldtheta2
        wz
    }), list( .link=link, .earg=earg ))))
}



dposbinom = function(x, size, prob, log = FALSE) {
    if(!is.Numeric(prob, positive=TRUE)) 
        stop("no zero or non-numeric values allowed for argument \"prob\"")
    L = max(length(x), length(size), length(prob))
    x = rep(x, len=L); size = rep(size, len=L); prob = rep(prob, len=L);
    ifelse(x==0, 0, dbinom(x=x, size=size, prob=prob, log=log) /
          (1 - (1-prob)^size))
}


pposbinom = function(q, size, prob, lower.tail = TRUE, log.p = FALSE) {
    if(!is.Numeric(prob, positive=TRUE)) 
        stop("no zero or non-numeric values allowed for argument \"prob\"")
    L = max(length(q), length(size), length(prob))
    q = rep(q, len=L); size = rep(size, len=L); prob = rep(prob, len=L);
    ifelse(q<1, 0, (pbinom(q=q, size=size, prob=prob, lower.tail=lower.tail,
         log.p=log.p) - (1-prob)^size) / (1 - (1-prob)^size))
}

qposbinom = function(p, size, prob, lower.tail = TRUE, log.p = FALSE) {
    if(!is.Numeric(prob, positive=TRUE)) 
        stop("no zero or non-numeric values allowed for argument \"prob\"")
    if(!is.Numeric(p, posit=TRUE) || any(p >= 1))
        stop("bad input for argument \"p\"")
    qbinom(p=p * (1 - (1-prob)^size) + (1-prob)^size, size=size, prob=prob,
           lower.tail=lower.tail, log.p=log.p)
}


rposbinom = function(n, size, prob) {
    if(!is.Numeric(n, posit=TRUE, allow=1, integ=TRUE))
        stop("bad input for argument \"n\"")
    if(!is.Numeric(prob, positive=TRUE)) 
        stop("no zero or non-numeric values allowed for argument \"prob\"")
    ans = rbinom(n=n, size=size, prob=prob)
    index = ans == 0
    size = rep(size, len=length(ans))
    prob = rep(prob, len=length(ans))
    while(any(index)) {
        more = rbinom(n=sum(index), size[index], prob=prob[index])
        ans[index] = more
        index = ans == 0
    }
    ans
}





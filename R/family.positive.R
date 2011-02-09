# These functions are
# Copyright (C) 1998-2011 T.W. Yee, University of Auckland.
# All rights reserved.







dposnegbin = function(x, size, prob = NULL, munb = NULL, log = FALSE) {
    if (length(munb)) {
        if (length(prob))
            stop("'prob' and 'munb' both specified")
        prob <- size/(size + munb)
    }
    if (!is.logical(log.arg <- log)) stop("bad input for 'log'")
    rm(log)

    L = max(length(x), length(prob), length(size))
    x = rep(x, len=L); prob = rep(prob, len=L); size = rep(size, len=L);

    ans = dnbinom(x=x, size=size, prob=prob, log=log.arg)
    index0 = x==0

    if (log.arg) {
        ans[ index0] = log(0.0)
        ans[!index0] = ans[!index0] - log1p(-dnbinom(x=0 * x[!index0],
                       size=size[!index0], prob=prob[!index0]))
    } else {
        ans[ index0] = 0.0
        ans[!index0] = ans[!index0] / pnbinom(q=0 * x[!index0],
                       size=size[!index0], prob=prob[!index0], lower.tail=FALSE)
    }
    ans
}


pposnegbin = function(q, size, prob = NULL, munb = NULL) {
    if (length(munb)) {
        if (length(prob))
            stop("'prob' and 'munb' both specified")
        prob <- size/(size + munb)
    }
    L = max(length(q), length(prob), length(size))
    q = rep(q, len=L); prob = rep(prob, len=L); size = rep(size, len=L);

    ifelse(q < 1, 0, (pnbinom(q,   size=size, prob=prob) -
                      dnbinom(q*0, size=size, prob=prob)) / pnbinom(q*0,
                            size=size, prob=prob, lower.tail = FALSE))
}


qposnegbin = function(p, size, prob = NULL, munb = NULL) {
    if (length(munb)) {
        if (length(prob))
            stop("'prob' and 'munb' both specified")
        prob <- size/(size + munb)
    }
    if (!is.Numeric(p, posit = TRUE) || any(p >= 1))
        stop("bad input for argument 'p'")
    qnbinom(p * pnbinom(q=p*0, size=size, prob=prob, lower.tail = FALSE) +
            dnbinom(x=p*0, size=size, prob=prob), size=size, prob=prob)
}


rposnegbin = function(n, size, prob = NULL, munb = NULL) {
    use.n = if ((length.n <- length(n)) > 1) length.n else
            if (!is.Numeric(n, integ = TRUE, allow = 1, posit = TRUE))
                stop("bad input for argument 'n'") else n

    if (length(munb)) {
        if (length(prob))
            stop("'prob' and 'munb' both specified")
        prob <- size/(size + munb)
    }
    ans = rnbinom(use.n, size=size, prob=prob)
    index = (ans == 0)
    size = rep(size, len=length(ans))
    prob = rep(prob, len=length(ans))
    while(any(index)) {
        more = rnbinom(n=sum(index), size=size[index], prob=prob[index])
        ans[index] = more
        index = (ans == 0)
    }
    ans
}







 posnegbinomial = function(lmunb = "loge", lk = "loge",
                           emunb =list(), ek = list(),
                           ik = NULL, zero = -2, cutoff = 0.995,
                           shrinkage.init = 0.95, method.init = 1)
{


    lkayy <- lk
    ekayy <- ek
    ikayy <- ik

    if (!is.Numeric(cutoff, allow = 1) || cutoff<0.8 || cutoff>=1)
        stop("range error in the argument 'cutoff'")
    if (!is.Numeric(method.init, allow = 1, integ = TRUE, posit = TRUE) ||
       method.init > 2) stop("argument 'method.init' must be 1 or 2")
    if (length(ikayy) && !is.Numeric(ikayy, posit = TRUE))
        stop("bad input for argument 'ik'")
    if (!is.Numeric(shrinkage.init, allow = 1) || shrinkage.init < 0 ||
       shrinkage.init > 1) stop("bad input for argument 'shrinkage.init'")

    if (mode(lmunb) != "character" && mode(lmunb) != "name")
        lmunb = as.character(substitute(lmunb))
    if (mode(lkayy) != "character" && mode(lkayy) != "name")
        lkayy = as.character(substitute(lkayy))
    if (!is.list(emunb)) emunb = list()
    if (!is.list(ek)) ek = list()

    new("vglmff",
    blurb = c("Positive-negative binomial distribution\n\n",
              "Links:    ",
              namesof("munb", lmunb, earg = emunb ), ", ",
              namesof("k",    lkayy, earg = ek ), "\n",
              "Mean:     munb / (1 - (k/(k+munb))^k)"),
    constraints = eval(substitute(expression({

        dotzero = .zero
        Musual = 2
        eval(negzero.expression)
    }), list( .zero = zero ))),
    initialize = eval(substitute(expression({
        if (any(y==0)) stop("there are zero values in the response")
        y = as.matrix(y) 
        M = 2 * ncol(y) 
        extra$NOS = NOS = ncoly = ncol(y)  # Number of species
        predictors.names = c(
            namesof(if (NOS == 1) "munb" else paste("munb", 1:NOS, sep = ""),
                    .lmunb, earg = .emunb, tag = FALSE),
            namesof(if (NOS == 1) "k"    else paste("k",    1:NOS, sep = ""),
                    .lkayy, earg = .ekayy, tag = FALSE))
        predictors.names = predictors.names[interleave.VGAM(M, M = 2)]
        if (!length(etastart)) {
            mu.init = y
            for(iii in 1:ncol(y)) {
                use.this = if ( .method.init == 2) {
                    weighted.mean(y[,iii], w)
                } else {
                    median(y[,iii])
                }
                mu.init[,iii] = (1 - .sinit) * y[,iii] + .sinit * use.this
            }

            if ( is.Numeric( .ikayy )) {
                kmat0 = matrix( .ikayy, nr = n, nc = NOS, byrow = TRUE)
            } else {
                posnegbinomial.Loglikfun =
                    function(kmat, y, x, w, extraargs) {
                    munb = extraargs
                    sum(w * dposnegbin(x = y, size = kmat, munb = munb,
                                       log = TRUE))
                    }
                k.grid = 2^((-6):6)
                kmat0 = matrix(0, nr=n, nc=NOS)
                for(spp. in 1:NOS) {
                    kmat0[,spp.] = getMaxMin(k.grid,
                                      objfun=posnegbinomial.Loglikfun,
                                      y=y[,spp.], x=x, w=w,
                                      extraargs= mu.init[,spp.])
                }
            }
            p00 = (kmat0 / (kmat0 + mu.init))^kmat0
            etastart =
              cbind(theta2eta(mu.init*(1-p00), .lmunb, earg = .emunb ),
                    theta2eta(kmat0,           .lkayy, earg = .ekayy ))
            etastart = etastart[,interleave.VGAM(M, M = 2), drop = FALSE]
        }
    }), list( .lmunb = lmunb, .lkayy = lkayy, .ikayy = ikayy,
              .emunb = emunb, .ekayy = ekayy,
              .sinit = shrinkage.init,
              .method.init = method.init ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        NOS = ncol(eta) / 2
        munb = eta2theta(eta[,2*(1:NOS)-1, drop = FALSE],
                         .lmunb, earg = .emunb )
        kmat = eta2theta(eta[,2*(1:NOS),   drop = FALSE],
                         .lkayy, earg = .ekayy )
        p0 = (kmat / (kmat + munb))^kmat
        munb / (1 - p0)
    }, list( .lkayy = lkayy, .lmunb = lmunb,
             .ekayy = ekayy, .emunb = emunb ))),
    last = eval(substitute(expression({
        temp0303 = c(rep( .lmunb, length = NOS),
                     rep( .lkayy, length = NOS))
        names(temp0303) =
           c(if (NOS == 1) "munb" else paste("munb", 1:NOS, sep = ""),
             if (NOS == 1) "k"    else paste("k",    1:NOS, sep = ""))
        temp0303 = temp0303[interleave.VGAM(M, M = 2)]
        misc$link = temp0303  # Already named
        misc$earg = vector("list", 2*NOS)
        names(misc$earg) = names(misc$link)
        for(ii in 1:NOS) {
            misc$earg[[2*ii-1]] = .emunb
            misc$earg[[2*ii  ]] = .ekayy
        }
        misc$cutoff = .cutoff 
        misc$method.init = .method.init
    }), list( .lmunb = lmunb, .lkayy = lkayy,
              .emunb = emunb, .ekayy = ekayy,
              .cutoff = cutoff, .method.init = method.init ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE,eta, extra = NULL) {
        NOS = ncol(eta) / 2
        munb = eta2theta(eta[,2*(1:NOS)-1, drop = FALSE],
                         .lmunb, earg = .emunb )
        kmat = eta2theta(eta[,2*(1:NOS),   drop = FALSE],
                         .lkayy, earg = .ekayy )
        if (residuals)
          stop("loglikelihood residuals not implemented yet") else {
          sum(w * dposnegbin(x = y, size = kmat, munb = munb, log = TRUE))
        }

    }, list( .lmunb = lmunb, .lkayy = lkayy,
             .emunb = emunb, .ekayy = ekayy ))),
    vfamily = c("posnegbinomial"),
    deriv = eval(substitute(expression({
        NOS = extra$NOS
        munb = eta2theta(eta[,2*(1:NOS)-1, drop = FALSE],
                         .lmunb, earg = .emunb )
        kmat = eta2theta(eta[,2*(1:NOS),   drop = FALSE],
                         .lkayy, earg = .ekayy )
        d3 = deriv3(~ -log(1 - (kmat. / (kmat. + munb. ))^kmat. ),
                    c("munb.", "kmat."), hessian = TRUE) # Extra term
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
        dl.dkayy = digamma(y+kmat) - digamma(kmat) -
                   (y+kmat) / (munb+kmat) + 1 +
                   log(kmat /(kmat+munb)) + dl0.dthetas[,,2]
        dmunb.deta = dtheta.deta(munb, .lmunb, earg = .emunb )
        dkayy.deta = dtheta.deta(kmat, .lkayy, earg = .ekayy )
        myderiv = w * cbind(dl.dmunb * dmunb.deta,
                            dl.dkayy * dkayy.deta)
        myderiv[, interleave.VGAM(M, M = 2)]
    }), list( .lmunb = lmunb, .lkayy = lkayy,
              .emunb = emunb, .ekayy = ekayy ))),
    weight = eval(substitute(expression({
        wz = matrix(0, n, 4*NOS-1)  # wz is no longer 'diagonal' 
        p0 = (kmat / (kmat + munb))^kmat
        ed2l.dmunb2 = (1/munb -
                      (munb + kmat*(1-p0))/(munb+kmat)^2) / (1-p0) -
                      d2l0.dthetas2[,,1]
        fred = dotFortran(name="enbin8", ans=double(n*NOS),
                      as.double(kmat), as.double(kmat/(munb+kmat)),
                      as.double(.cutoff), as.integer(n),
                      ok=as.integer(1), as.integer(NOS), sumpdf=double(1),
                      macheps = as.double(.Machine$double.eps))
        if (fred$ok != 1)
            stop("error in Fortran subroutine exnbin")
        dim(fred$ans) = c(n, NOS)
        ed2l.dk2 = -fred$ans/(1-p0) - 1/kmat + 1/(kmat+munb) -
                   munb * p0 / ((1-p0)*(munb+kmat)^2) - d2l0.dthetas2[,,2]
        wz[,2*(1:NOS)-1] = dmunb.deta^2 * ed2l.dmunb2
        wz[,2*(1:NOS)  ] = dkayy.deta^2 * ed2l.dk2
        wz[,2*NOS+2*(1:NOS)-1] = -d2l0.dthetas2[,,3] *
                                 dmunb.deta * dkayy.deta
        w * wz
    }), list( .cutoff = cutoff ))))
}





dpospois = function(x, lambda, log = FALSE) {
    if (!is.logical(log.arg <- log)) stop("bad input for 'log'")
    rm(log)

    if (!is.Numeric(lambda, posit = TRUE))
        stop("bad input for argument 'lambda'")
    L = max(length(x), length(lambda))
    x = rep(x, len=L); lambda = rep(lambda, len=L); 
    ans = if (log.arg) {
        ifelse(x == 0, log(0.0), dpois(x, lambda, log = TRUE) -
               log1p(-exp(-lambda)))
    } else {
        ifelse(x == 0, 0, -dpois(x, lambda) / expm1(-lambda))
    }
    ans
}


ppospois = function(q, lambda) {
    if (!is.Numeric(lambda, posit = TRUE))
        stop("bad input for argument 'lambda'")
    L = max(length(q), length(lambda))
    q = rep(q, len=L); lambda = rep(lambda, len=L); 
    ifelse(q < 1, 0, (ppois(q, lambda) - exp(-lambda)) / (-expm1(-lambda)))
}

qpospois = function(p, lambda) {
    if (!is.Numeric(lambda, posit = TRUE))
        stop("bad input for argument 'lambda'")
    if (!is.Numeric(p, posit = TRUE) || any(p >= 1))
        stop("bad input for argument 'p'")
    qpois(p * (-expm1(-lambda)) + exp(-lambda), lambda)
}


rpospois = function(n, lambda) {
    use.n = if ((length.n <- length(n)) > 1) length.n else
            if (!is.Numeric(n, integ = TRUE, allow = 1, posit = TRUE))
                stop("bad input for argument 'n'") else n

    if (any(lambda == 0))
        stop("no zero values allowed for argument 'lambda'")
    ans = rpois(use.n, lambda)
    lambda = rep(lambda, len=use.n)
    index = (ans == 0)
    while(any(index)) {
        more = rpois(n=sum(index), lambda[index])
        ans[index] = more
        index = (ans == 0)
    }
    ans
}




 pospoisson = function(link = "loge", earg = list(), expected = TRUE,
                       ilambda = NULL, method.init = 1)
{

  if (!missing(link))
    link <- as.character(substitute(link))
  if (!is.list(earg)) earg <- list()

  if (!is.logical(expected) || length(expected) != 1)
    stop("bad input for argument 'expected'")
  if (length( ilambda) && !is.Numeric(ilambda, posit = TRUE))
    stop("bad input for argument 'ilambda'")
  if (!is.Numeric(method.init, allow = 1, integ = TRUE, posit = TRUE) ||
    method.init > 3) stop("argument 'method.init' must be 1 or 2 or 3")

  new("vglmff",
  blurb = c("Positive-Poisson distribution\n\n",
            "Links:    ",
            namesof("lambda", link, earg = earg, tag = FALSE)),
  initialize = eval(substitute(expression({
    y <- as.matrix(y)

    if (any(y < 1))
        stop("all y values must be in 1,2,3,...")
    if (any(y != round(y )))
        stop("the response must be integer-valued")

    predictors.names <- namesof(
      paste("lambda", if (ncol(y) > 1) 1:ncol(y) else "", sep = ""),
      .link, earg = .earg, tag = FALSE)

    if ( .method.init == 1) {
      lambda.init <- apply(y, 2, median) + 1/8
      lambda.init <- matrix(lambda.init, n, ncol(y), byrow = TRUE)
    } else if ( .method.init == 2) {
      lambda.init <- apply(y, 2, weighted.mean, w=w) + 1/8
      lambda.init <- matrix(lambda.init, n, ncol(y), byrow = TRUE)
    } else {
      lambda.init <- -y / expm1(-y)
    }
    if (length( .ilambda))
      lambda.init <- lambda.init * 0 + .ilambda
    if (!length(etastart))
        etastart <- theta2eta(lambda.init, .link, earg = .earg)
  }), list( .link = link, .earg = earg,
            .ilambda = ilambda, .method.init = method.init ))),
  inverse = eval(substitute(function(eta, extra = NULL) {
    lambda <- eta2theta(eta, .link, earg = .earg )
    -lambda / expm1(-lambda)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$expected <- .expected
    misc$link <- rep( .link, len = M)
    names(misc$link) <- if (M == 1) "lambda" else
                       paste("lambda", 1:M, sep = "")
    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for(ii in 1:M)
        misc$earg[[ii]] <- .earg
  }), list( .link = link, .earg = earg, .expected = expected ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    lambda <- eta2theta(eta, .link, earg = .earg ) 
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      sum(w * dpospois(x = y, lambda = lambda, log = TRUE))
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("pospoisson"),
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta, .link, earg = .earg ) 
    temp6 <- expm1(lambda)
    dl.dlambda <- y / lambda - 1 - 1 / temp6
    dlambda.deta <- dtheta.deta(lambda, .link, earg = .earg )
    w * dl.dlambda * dlambda.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    if ( .expected ) {
      ed2l.dlambda2 <- (temp6 + 1) * (1/lambda - 1/temp6) / temp6
      wz <- (dlambda.deta^2) * ed2l.dlambda2
    } else {
      d2l.dlambda2 <- y / lambda^2 - (temp6 + 1) / temp6^2
      d2lambda.deta2 <- d2theta.deta2(lambda, .link, earg = .earg)
      wz <- (dlambda.deta^2) * d2l.dlambda2 - dl.dlambda * d2lambda.deta2
    }
    w * wz
  }), list( .link = link, .earg = earg, .expected = expected ))))
}




 posbinomial = function(link = "logit", earg = list()) {

    if (!missing(link))
        link = as.character(substitute(link))
    if (!is.list(earg)) earg = list()
       
    new("vglmff",
    blurb = c("Positive-binomial distribution\n\n",
              "Links:    ",
           namesof("prob", link, earg = earg, tag = FALSE), "\n"),
    initialize = eval(substitute(expression({

        mustart.orig = mustart
      	eval(binomialff(link = .link, earg = .earg )@initialize)
        predictors.names = namesof("prob", .link, earg = .earg , tag = FALSE)
	if (length(extra)) extra$w = w else extra = list(w = w)
        if (!length(etastart)) {
	  mustart.use = if (length(mustart.orig)) mustart.orig else mustart
	  etastart = cbind(theta2eta(mustart.use, .link, earg = .earg ))
        }
        mustart = NULL
    }), list( .link = link, .earg = earg ))),
    inverse = eval(substitute(function(eta, extra = NULL) {
        w = extra$w
        mymu = eta2theta(eta, .link, earg = .earg )
        nvec = if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                  round(w)
        mymu / (1-(1-mymu)^(nvec))
    },
    list( .link = link, .earg = earg ))),
    last = eval(substitute(expression({
        extra$w = NULL   # Kill it off 
        misc$link = c(prob = .link)
        misc$earg = list(prob = .earg )
    }), list( .link = link, .earg = earg ))),
    loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta, extra = NULL) {

          ycounts = if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                    y * w # Convert proportions to counts
          nvec = if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                    round(w)
          use.orig.w = if (is.numeric(extra$orig.w)) extra$orig.w else 1
        mymu = eta2theta(eta, .link, earg = .earg )

        if (residuals) stop("loglikelihood residuals ",
                            "not implemented yet") else {

          sum(use.orig.w * dposbinom(x = ycounts, size = nvec,
                                     prob = mymu, log = TRUE))
        }
    }, list( .link = link, .earg = earg ))),
    vfamily = c("posbinomial"),
    deriv = eval(substitute(expression({
          use.orig.w = if (is.numeric(extra$orig.w)) extra$orig.w else
                       rep(1, n)
          nvec = if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                    round(w)

        mymu = eta2theta(eta, .link, earg = .earg )
        dl.dmymu = y / mymu - (1 - y) / (1 - mymu) -
                   ((1 - mymu)^(nvec-1)) / (1 - (1 - mymu)^nvec)

        dmymu.deta = dtheta.deta(mymu, .link, earg = .earg )
        w * dl.dmymu * dmymu.deta
    }), list( .link = link, .earg = earg ))),
    weight = eval(substitute(expression({
        temp1 = 1 - (1-mymu)^nvec
        temp2 =     (1-mymu)^2
        ed2l.dmymu2 = -1/(mymu*temp1) - 1/temp2 + mymu/(temp2*temp1) +
                       (nvec-1) * (1-mymu)^(nvec-2) / temp1 +
                       nvec * (temp2^(nvec-1)) / temp1^2
        wz = -w * (dmymu.deta^2) * ed2l.dmymu2
        wz
    }), list( .link = link, .earg = earg ))))
}









dposbinom = function(x, size, prob, log = FALSE) {
    log.arg = log
    rm(log)
    L = max(length(x), length(size), length(prob))
    x = rep(x, len=L); size = rep(size, len=L); prob = rep(prob, len=L);

    answer = NaN * x
    is0 <- (x == 0)
    ok2 <- prob > 0 & prob <= 1 & size == round(size) & size > 0
    answer = dbinom(x = x, size = size, prob = prob, log = TRUE) -
             log1p(-dbinom(x=0*x, size=size, prob=prob))
    answer[!ok2] = NaN
    if (log.arg) {
        answer[is0 & ok2]  = log(0.0)
    } else {
        answer = exp(answer)
        answer[is0 & ok2] = 0.0
    }
    answer
}

pposbinom = function(q, size, prob, lower.tail = TRUE, log.p = FALSE) {
    if (!is.Numeric(prob, positive = TRUE)) 
        stop("no zero or non-numeric values allowed for argument 'prob'")
    L = max(length(q), length(size), length(prob))
    q = rep(q, len=L); size = rep(size, len=L); prob = rep(prob, len=L);
    ifelse(q < 1, 0, (pbinom(q=q, size=size, prob=prob, lower.tail=lower.tail,
         log.p=log.p) - (1-prob)^size) / (1 - (1-prob)^size))
}

qposbinom = function(p, size, prob, lower.tail = TRUE, log.p = FALSE) {
    if (!is.Numeric(prob, positive = TRUE)) 
        stop("no zero or non-numeric values allowed for argument 'prob'")
    if (!is.Numeric(p, posit = TRUE) || any(p >= 1))
        stop("bad input for argument 'p'")
    qbinom(p=p * (1 - (1-prob)^size) + (1-prob)^size, size=size, prob=prob,
           lower.tail=lower.tail, log.p=log.p)
}

rposbinom = function(n, size, prob) {
    use.n = if ((length.n <- length(n)) > 1) length.n else
            if (!is.Numeric(n, integ = TRUE, allow = 1, posit = TRUE))
                stop("bad input for argument 'n'") else n

    if (any(prob == 0))
        stop("no zero values allowed for argument 'prob'")
    ans = rbinom(n=use.n, size=size, prob=prob)
    index = (ans == 0)
    size = rep(size, len=length(ans))
    prob = rep(prob, len=length(ans))
    while(any(index)) {
        more = rbinom(n=sum(index), size[index], prob=prob[index])
        ans[index] = more
        index = (ans == 0)
    }
    ans
}





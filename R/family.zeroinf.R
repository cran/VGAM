# These functions are
# Copyright (C) 1998-2012 T.W. Yee, University of Auckland.
# All rights reserved.











dzanegbin = function(x, size, prob = NULL, munb = NULL, pobs0 = 0,
                     log = FALSE) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  LLL = max(length(x), length(pobs0), length(prob), length(size))
  if (length(x)     != LLL) x     = rep(x,    len = LLL)
  if (length(pobs0) != LLL) pobs0 = rep(pobs0,  len = LLL);
  if (length(prob)  != LLL) prob  = rep(prob, len = LLL)
  if (length(size)  != LLL) size  = rep(size, len = LLL);

  ans = rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  if (!is.Numeric(prob, positive = TRUE))
    stop("argument 'prob' must be in [0,Inf)")
  if (!is.Numeric(size, positive = TRUE))
    stop("argument 'size' must be in [0,Inf)")
  index0 = x == 0

  if (log.arg) {
    ans[ index0] = log(pobs0[index0])
    ans[!index0] = log1p(-pobs0[!index0]) +
                   dposnegbin(x[!index0], prob = prob[!index0],
                              size = size[!index0], log = TRUE)
  } else {
    ans[ index0] = pobs0[index0]
    ans[!index0] = (1-pobs0[!index0]) * dposnegbin(x[!index0],
                     prob = prob[!index0], size = size[!index0])
  }
  ans
}



pzanegbin = function(q, size, prob = NULL, munb = NULL, pobs0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  LLL = max(length(q), length(pobs0), length(prob), length(size))
  if (length(q)     != LLL) q     = rep(q,     len = LLL);
  if (length(pobs0) != LLL) pobs0 = rep(pobs0, len = LLL);
  if (length(prob)  != LLL) prob  = rep(prob,  len = LLL);
  if (length(size)  != LLL) size  = rep(size,  len = LLL);
  ans = rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  qindex = (q >  0)
  ans[ qindex] = pobs0[qindex] + (1 - pobs0[qindex]) *
                 pposnegbin(q[qindex], size = size[qindex],
                                       prob = prob[qindex])
  ans[q <  0] = 0
  ans[q == 0] = pobs0[q == 0]
  ans
}


qzanegbin = function(p, size, prob = NULL, munb = NULL, pobs0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size/(size + munb)
  }

  LLL = max(length(p), length(pobs0), length(prob), length(size))
  if (length(p)     != LLL) p      = rep(p,     len = LLL);
  if (length(pobs0) != LLL) pobs0  = rep(pobs0, len = LLL);
  if (length(prob)  != LLL) prob   = rep(prob,  len = LLL);
  if (length(size)  != LLL) size   = rep(size,  len = LLL);
  ans = rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ans = p
  ans[p <= pobs0] = 0
  pindex = (p > pobs0)
  ans[pindex] = qposnegbin((p[pindex] -
                            pobs0[pindex]) / (1 - pobs0[pindex]),
                            prob = prob[pindex],
                            size = size[pindex])
  ans
}


rzanegbin = function(n, size, prob = NULL, munb = NULL, pobs0 = 0) {
  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE,
                          allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  ans = rposnegbin(n = use.n, prob = prob, size = size)
  if (length(pobs0) != use.n) pobs0 = rep(pobs0, len = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")

  ifelse(runif(use.n) < pobs0, 0, ans)
}





dzapois = function(x, lambda, pobs0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  LLL = max(length(x), length(lambda), length(pobs0))
  if (length(x)      != LLL) x      = rep(x,      len = LLL);
  if (length(lambda) != LLL) lambda = rep(lambda, len = LLL);
  if (length(pobs0)  != LLL) pobs0  = rep(pobs0,  len = LLL);
  ans = rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  index0 = (x == 0)

  if (log.arg) {
    ans[ index0] = log(pobs0[index0])
    ans[!index0] = log1p(-pobs0[!index0]) +
                   dpospois(x[!index0], lambda[!index0], log = TRUE)
  } else {
    ans[ index0] = pobs0[index0]
    ans[!index0] = (1 - pobs0[!index0]) *
                   dpospois(x[!index0], lambda[!index0])
  }
  ans
}



pzapois = function(q, lambda, pobs0 = 0) {
  LLL = max(length(q), length(lambda), length(pobs0))
  if (length(q)      != LLL) q      = rep(q,      len = LLL);
  if (length(lambda) != LLL) lambda = rep(lambda, len = LLL);
  if (length(pobs0)  != LLL) pobs0  = rep(pobs0,  len = LLL);
  ans = rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  ans[q >  0] =    pobs0[q > 0] +
                (1-pobs0[q > 0]) * ppospois(q[q > 0], lambda[q > 0])
  ans[q <  0] = 0
  ans[q == 0] = pobs0[q == 0]
  ans
}


qzapois = function(p, lambda, pobs0 = 0) {
  LLL = max(length(p), length(lambda), length(pobs0))
  if (length(p)      != LLL) p      = rep(p,      len = LLL);
  if (length(lambda) != LLL) lambda = rep(lambda, len = LLL);
  if (length(pobs0)  != LLL) pobs0  = rep(pobs0,  len = LLL);

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ans = p
  ind4 = (p > pobs0)
  ans[!ind4] = 0
  ans[ ind4] = qpospois((p[ind4] - pobs0[ind4]) / (1 - pobs0[ind4]),
                        lambda = lambda[ind4])
  ans
}


rzapois = function(n, lambda, pobs0 = 0) {
  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE,
                          allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  ans = rpospois(use.n, lambda)
  if (length(pobs0) != use.n)
    pobs0 = rep(pobs0, length = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")

  ifelse(runif(use.n) < pobs0, 0, ans)
}





dzipois = function(x, lambda, pstr0 = 0, log = FALSE) {



  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  LLL = max(length(x), length(lambda), length(pstr0))
  if (length(x)      != LLL) x      = rep(x,      len = LLL);
  if (length(lambda) != LLL) lambda = rep(lambda, len = LLL);
  if (length(pstr0)  != LLL) pstr0  = rep(pstr0,  len = LLL);

  ans = x + lambda + pstr0



  index0 = (x == 0)
  if (log.arg) {
    ans[ index0] = log(pstr0[ index0] + (1 - pstr0[ index0]) *
                       dpois(x[ index0], lambda[ index0]))
    ans[!index0] = log1p(-pstr0[!index0]) +
                   dpois(x[!index0], lambda[!index0], log = TRUE)
  } else {
    ans[ index0] =      pstr0[ index0] + (1 - pstr0[ index0]) *
                       dpois(x[ index0], lambda[ index0])
    ans[!index0] = (1 - pstr0[!index0]) * dpois(x[!index0], lambda[!index0])
  }


  deflat_limit = -1 / expm1(lambda)
  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN

  ans
}


pzipois = function(q, lambda, pstr0 = 0) {

  LLL = max(length(pstr0), length(lambda), length(q))
  if (length(pstr0)  != LLL) pstr0  = rep(pstr0,  len = LLL);
  if (length(lambda) != LLL) lambda = rep(lambda, len = LLL);
  if (length(q)      != LLL) q      = rep(q,      len = LLL);

  ans = ppois(q, lambda)
  ans = ifelse(q < 0, 0, pstr0 + (1 - pstr0) * ans)


  deflat_limit = -1 / expm1(lambda)
  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN


  ans
}


qzipois = function(p, lambda, pstr0 = 0) {

  LLL = max(length(p), length(lambda), length(pstr0))
  ans =
  p      = rep(p,      len = LLL)
  lambda = rep(lambda, len = LLL)
  pstr0  = rep(pstr0,  len = LLL)

  ans[p <= pstr0] = 0 
  pindex = (p > pstr0)
  ans[pindex] = qpois((p[pindex] - pstr0[pindex]) / (1 - pstr0[pindex]),
                      lambda = lambda[pindex])


  deflat_limit = -1 / expm1(lambda)
  ind0 = (deflat_limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 = pstr0[ind0] + (1 - pstr0[ind0]) * exp(-lambda[ind0])
    ans[p[ind0] <= pobs0] = 0 
    pindex = (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 = pstr0[pindex] + (1 - pstr0[pindex]) * exp(-lambda[pindex])
    ans[pindex] = qpospois((p[pindex] - Pobs0) / (1 - Pobs0),
                           lambda = lambda[pindex])
  }


  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN


  ans[p < 0] = NaN
  ans[p > 1] = NaN
  ans
}


rzipois = function(n, lambda, pstr0 = 0) {

  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE,
                          allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (length(pstr0)  != use.n) pstr0  = rep(pstr0,  len = use.n);
  if (length(lambda) != use.n) lambda = rep(lambda, len = use.n);
 
  ans = rpois(use.n, lambda)
  ans = ifelse(runif(use.n) < pstr0, 0, ans)



  prob0 = exp(-lambda)
  deflat_limit = -1 / expm1(lambda)
  ind0 = (deflat_limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 = pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] = rpospois(sum(ind0), lambda[ind0]) 
    ans[ind0] = ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN

  ans
}






 yip88 = function(link.lambda = "loge", n.arg = NULL) {




  if (mode(link.lambda) != "character" && mode(link.lambda) != "name")
    link.lambda = as.character(substitute(link.lambda))

  new("vglmff",
  blurb = c("Zero-inflated Poisson (based on Yip (1988))\n\n",
            "Link:     ", namesof("lambda", link.lambda), "\n",
            "Variance: (1 - pstr0) * lambda"),
  first = eval(substitute(expression({
    zero <- y == 0
    if (any(zero)) {
      if (length(extra)) extra$sumw = sum(w) else
        extra = list(sumw=sum(w))
      if (is.numeric(.n.arg) && extra$sumw != .n.arg) 
        stop("value of 'n.arg' conflicts with data ",
             "(it need not be specified anyway)")
      warning("trimming out the zero observations")

      axa.save =  attr(x, "assign")
      x = x[!zero,, drop = FALSE]
      attr(x, "assign") = axa.save    # Don't lose these!!
      w = w[!zero]
      y = y[!zero]
    } else {
      if (!is.numeric(.n.arg)) 
        stop("n.arg must be supplied")
    }
        
  }), list( .n.arg = n.arg ))),

  initialize = eval(substitute(expression({
    narg = if (is.numeric(.n.arg)) .n.arg else extra$sumw
    if (sum(w) > narg)
      stop("sum(w) > narg")

    predictors.names = namesof("lambda", .link.lambda, tag = FALSE)
    if (!length(etastart)) {
      lambda.init = rep(median(y), length = length(y))
      etastart = theta2eta(lambda.init, .link.lambda)
    }
    if (length(extra)) {
      extra$sumw = sum(w)
      extra$narg = narg   # For @linkinv
    } else {
      extra = list(sumw = sum(w), narg = narg)
    }
  }), list( .link.lambda = link.lambda, .n.arg = n.arg ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    lambda = eta2theta(eta, .link.lambda)
    temp5 = exp(-lambda)
    pstr0 = (1 - temp5 - extra$sumw/extra$narg) / (1 - temp5)
    if (any(pstr0 <= 0))
      stop("non-positive value(s) of pstr0")
    (1-pstr0) * lambda
  }, list( .link.lambda = link.lambda ))),

  last = eval(substitute(expression({
    misc$link = c(lambda = .link.lambda )

    if (intercept.only) {
      suma = extra$sumw
      pstr0 = (1 - temp5[1] - suma / narg) / (1 - temp5[1])
      pstr0 = if (pstr0 < 0 || pstr0 > 1) NA else pstr0
      misc$pstr0 = pstr0
    }
  }), list( .link.lambda = link.lambda ))),

  loglikelihood = eval(substitute(function(mu, y, w, residuals = FALSE,
                                           eta, extra = NULL) {
    lambda = eta2theta(eta, .link.lambda)
    temp5 = exp(-lambda)
    pstr0 = (1 - temp5 - extra$sumw / extra$narg) / (1 - temp5)
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(w * dzipois(x = y, pstr0 = pstr0, lambda = lambda, log = TRUE))
    }
  }, list( .link.lambda = link.lambda ))),

  vfamily = c("yip88"),
  deriv = eval(substitute(expression({
    lambda = eta2theta(eta, .link.lambda)
    temp5 = exp(-lambda)
    dl.dlambda = -1 + y/lambda - temp5/(1-temp5)
    dlambda.deta = dtheta.deta(lambda, .link.lambda)
    w * dl.dlambda * dlambda.deta
  }), list( .link.lambda = link.lambda ))),
  weight = eval(substitute(expression({
    d2lambda.deta2 = d2theta.deta2(lambda, .link.lambda)
    d2l.dlambda2 = -y / lambda^2 + temp5 / (1 - temp5)^2
    -w * (d2l.dlambda2*dlambda.deta^2 + dl.dlambda*d2lambda.deta2)
  }), list( .link.lambda = link.lambda ))))
}




 zapoisson = function(lpobs0 = "logit", llambda = "loge",
                      epobs0 = list(),  elambda = list(), zero = NULL) {



  lpobs_0 = lpobs0
  epobs_0 = epobs0

  if (mode(lpobs_0) != "character" && mode(lpobs_0) != "name")
    lpobs_0 = as.character(substitute(lpobs_0))
  if (mode(llambda) != "character" && mode(llambda) != "name")
    llambda = as.character(substitute(llambda))

  if (!is.list(epobs_0)) epobs_0 = list()
  if (!is.list(elambda)) elambda = list()

  new("vglmff",
  blurb = c("Zero-altered Poisson ",
            "(Bernoulli and positive-Poisson conditional model)\n\n",
            "Links:    ",
            namesof("pobs0",  lpobs_0, earg = epobs_0, tag = FALSE), ", ",
            namesof("lambda", llambda, earg = elambda, tag = FALSE), "\n",
            "Mean:     (1 - pobs0) * lambda / (1 - exp(-lambda))"),

  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero)
  }, list( .zero = zero ))),
  initialize = eval(substitute(expression({
    Musual <- 2
    y <- as.matrix(y)
    if (any(y != round(y )))
      stop("the response must be integer-valued")
    if (any(y < 0))
      stop("the response must not have negative values")

    extra$y0 = y0 = ifelse(y == 0, 1, 0)
    extra$NOS = NOS = ncoly = ncol(y)  # Number of species
    extra$skip.these = skip.these = matrix(as.logical(y0), n, NOS)

    mynames1 = if (ncoly == 1) "pobs0"    else
               paste("pobs0",    1:ncoly, sep = "")
    mynames2 = if (ncoly == 1) "lambda" else
               paste("lambda", 1:ncoly, sep = "")
    predictors.names = 
        c(namesof(mynames1, .lpobs_0, earg = .epobs_0, tag = FALSE),
          namesof(mynames2, .llambda, earg = .elambda, tag = FALSE))[
          interleave.VGAM(Musual*NOS, M = Musual)]

    if (!length(etastart)) {
      etastart =
        cbind(theta2eta((0.5 + w*y0) / (1+w), .lpobs_0, earg = .epobs_0 ),
              matrix(1, n, NOS))  # 1 here is any old value
      for(spp. in 1:NOS) {
        sthese = skip.these[, spp.]
        etastart[!sthese, NOS+spp.] =
          theta2eta(y[!sthese, spp.] / (-expm1(-y[!sthese, spp.])),
                    .llambda, earg = .elambda )
      }
      etastart = etastart[, interleave.VGAM(ncol(etastart), M = Musual)]
    }
  }), list( .lpobs_0 = lpobs_0, .llambda = llambda,
            .epobs_0 = epobs_0, .elambda = elambda ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS = extra$NOS
    Musual <- 2


    pobs_0 = cbind(eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                             .lpobs_0, earg = .epobs_0 ))
    lambda = cbind(eta2theta(eta[, Musual*(1:NOS)-0, drop = FALSE],
                             .llambda, earg = .elambda ))

    (1 - pobs_0) * lambda / (-expm1(-lambda))
  }, list( .lpobs_0 = lpobs_0, .llambda = llambda,
           .epobs_0 = epobs_0, .elambda = elambda ))),
  last = eval(substitute(expression({
    temp.names = c(rep( .lpobs_0 , len = NOS),
                   rep( .llambda , len = NOS))
    temp.names = temp.names[interleave.VGAM(Musual*NOS, M = Musual)]
    misc$link  = temp.names
    misc$expected = TRUE
    misc$earg = vector("list", Musual * NOS)

    names(misc$link) <-
    names(misc$earg) <-
        c(mynames1, mynames2)[interleave.VGAM(Musual*NOS, M = Musual)]

    for(ii in 1:NOS) {
      misc$earg[[Musual*ii-1]] = .epobs_0
      misc$earg[[Musual*ii  ]] = .elambda
    }
  }), list( .lpobs_0 = lpobs_0, .llambda = llambda,
            .epobs_0 = epobs_0, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    NOS = extra$NOS
    Musual <- 2

    pobs0    = cbind(eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                             .lpobs_0, earg = .epobs_0))
    lambda = cbind(eta2theta(eta[, Musual*(1:NOS)-0, drop = FALSE],
                             .llambda, earg = .elambda ))

    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(w * dzapois(x = y, pobs0 = pobs0, lambda = lambda, log = TRUE))
    }
  }, list( .lpobs_0 = lpobs_0, .llambda = llambda,
           .epobs_0 = epobs_0, .elambda = elambda ))),
  vfamily = c("zapoisson"),
  deriv = eval(substitute(expression({
    Musual <- 2
    NOS = extra$NOS
    y0 = extra$y0
    skip = extra$skip.these

    phimat = cbind(eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                             .lpobs_0, earg = .epobs_0 ))
    lambda = cbind(eta2theta(eta[, Musual*(1:NOS)-0, drop = FALSE],
                             .llambda, earg = .elambda ))

    dl.dlambda = y / lambda + 1 / expm1(-lambda)
    dl.dphimat = -1 / (1 - phimat) # For y > 0 obsns

    for(spp. in 1:NOS) {
      dl.dphimat[skip[, spp.], spp.] = 1 / phimat[skip[, spp.], spp.]
      dl.dlambda[skip[, spp.], spp.] = 0
    }
    dlambda.deta = dtheta.deta(lambda, .llambda, earg = .elambda)
    mu.phi0 = phimat

    temp3 = if (.lpobs_0 == "logit") {
      c(w) * (y0 - mu.phi0)
    } else {
      c(w) * dtheta.deta(mu.phi0, link = .lpobs_0 , earg = .epobs_0 ) *
            dl.dphimat
    }

    ans <- cbind(temp3,
                 c(w) * dl.dlambda * dlambda.deta)
    ans = ans[, interleave.VGAM(ncol(ans), M = Musual)]
    ans
  }), list( .lpobs_0 = lpobs_0, .llambda = llambda,
            .epobs_0 = epobs_0, .elambda = elambda ))),
  weight = eval(substitute(expression({

    wz = matrix(0.0, n, Musual*NOS)



    temp5 = expm1(lambda)
    ed2l.dlambda2 = (1 - phimat) * (temp5 + 1) *
                    (1 / lambda - 1 / temp5) / temp5
    wz[, NOS+(1:NOS)] = w * ed2l.dlambda2 * dlambda.deta^2


    tmp100 = mu.phi0 * (1.0 - mu.phi0)
    tmp200 = if ( .lpobs_0 == "logit" && is.empty.list( .epobs_0 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi0, link = .lpobs_0, earg = .epobs_0)^2)
    }


  if (FALSE)
    for(ii in 1:NOS) {
      index200 = abs(tmp200[, ii]) < .Machine$double.eps
      if (any(index200)) {
        tmp200[index200, ii] = 10.0 * .Machine$double.eps^(3/4)
      }
    }


    wz[, 1:NOS] =  tmp200

    wz = wz[, interleave.VGAM(ncol(wz), M = Musual)]

    wz
  }), list( .lpobs_0 = lpobs_0,
            .epobs_0 = epobs_0 ))))
} #   End of zapoisson





zanegbinomial.control <- function(save.weight = TRUE, ...)
{
  list(save.weight = save.weight)
}



 zanegbinomial =
  function(lpobs0 = "logit", lmunb = "loge", lsize = "loge",
           epobs0 = list(),  emunb = list(), esize = list(),
           ipobs0 = NULL,                    isize = NULL,
           zero = c(-1, -3),
           imethod = 1,
           nsimEIM = 250,
           shrinkage.init = 0.95)
{






  if (!is.Numeric(nsimEIM, allowable.length = 1,
                  positive = TRUE, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 30)
    warning("argument 'nsimEIM' should be greater than 30, say")


  if (length(ipobs0) && (!is.Numeric(ipobs0, positive = TRUE) ||
     max(ipobs0) >= 1))
    stop("If given, argument 'ipobs0' must contain values in (0,1) only")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("If given, argument 'isize' must contain positive values only")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")
  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
     shrinkage.init < 0 ||
     shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")

  if (mode(lmunb) != "character" && mode(lmunb) != "name")
      lmunb = as.character(substitute(lmunb))
  if (mode(lsize) != "character" && mode(lsize) != "name")
      lsize = as.character(substitute(lsize))
  if (mode(lpobs0) != "character" && mode(lpobs0) != "name")
      lpobs0 = as.character(substitute(lpobs0))

  if (!is.list(epobs0)) epobs0 = list()
  if (!is.list(emunb)) emunb = list()
  if (!is.list(esize)) esize = list()



  new("vglmff",
  blurb = c("Zero-altered negative binomial (Bernoulli and\n",
            "positive-negative binomial conditional model)\n\n",
            "Links:    ",
            namesof("pobs0", lpobs0, earg = epobs0, tag = FALSE), ", ",
            namesof("munb",  lmunb,  earg = emunb,  tag = FALSE), ", ",
            namesof("size",  lsize,  earg = esize,  tag = FALSE), "\n",
            "Mean:     (1 - pobs0) * munb / (1 - (size / (size + ",
                                                  "munb))^size)"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 3
    eval(negzero.expression)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    Musual <- 3
    y <- as.matrix(y)
    extra$NOS = NOS = ncoly = ncol(y)  # Number of species
    M = Musual * ncoly # 

    if (any(y != round(y)))
      stop("the response must be integer-valued")
    if (any(y < 0))
      stop("the response must not have negative values")

    mynames1 = if (NOS == 1) "pobs0" else paste("pobs0", 1:NOS, sep = "")
    mynames2 = if (NOS == 1) "munb"  else paste("munb",  1:NOS, sep = "")
    mynames3 = if (NOS == 1) "size"  else paste("size",  1:NOS, sep = "")
    predictors.names =
        c(namesof(mynames1, .lpobs0 , earg = .epobs0 , tag = FALSE),
          namesof(mynames2, .lmunb  , earg = .emunb  , tag = FALSE),
          namesof(mynames3, .lsize  , earg = .esize  , tag = FALSE))[
          interleave.VGAM(Musual*NOS, M = Musual)]


    extra$y0 = y0 = ifelse(y == 0, 1, 0)
    extra$skip.these = skip.these = matrix(as.logical(y0), n, NOS)


    if (!length(etastart)) {
      mu.init = y
      for(iii in 1:ncol(y)) {
        index.posy = (y[, iii] > 0)
        use.this = if ( .imethod == 2) {
          weighted.mean(y[index.posy, iii], w[index.posy])
        } else {
          median(rep(y[index.posy, iii], w[index.posy])) + 1/2
        }
        mu.init[ index.posy, iii] = (1 - .sinit ) * y[index.posy, iii] +
                                        .sinit   * use.this
        mu.init[!index.posy, iii] = use.this
        max.use.this =  7 * use.this + 10
        vecTF = (mu.init[, iii] > max.use.this)
        if (any(vecTF))
          mu.init[vecTF, iii] = max.use.this
      }

      pnb0 = matrix(if (length( .ipobs0 )) .ipobs0 else -1,
                      nrow = n, ncol = NOS, byrow = TRUE)
      for(spp. in 1:NOS) {
        if (any(pnb0[, spp.] < 0)) {
          index.y0 = y[, spp.] < 0.5
          pnb0[, spp.] = max(min(sum(index.y0)/n, 0.97), 0.03)
        }
      }


      if ( is.Numeric( .isize )) {
        kmat0 = matrix( .isize , nrow = n, ncol = ncoly, byrow = TRUE)
      } else {
        posnegbinomial.Loglikfun = function(kmat, y, x, w, extraargs) {
         munb = extraargs
         sum(w * dposnegbin(x = y, munb = munb, size = kmat,
                            log = TRUE))
        }
        k.grid = 2^((-6):6)
        kmat0 = matrix(0, nrow = n, ncol = NOS) 
        for(spp. in 1:NOS) {
          index.posy = y[, spp.] > 0
          posy = y[index.posy, spp.]
          kmat0[, spp.] = getMaxMin(k.grid,
                                   objfun = posnegbinomial.Loglikfun,
                                   y = posy, x = x[index.posy,],
                                   w = w[index.posy],
                                   extraargs = mu.init[index.posy, spp.])
        }
      }

      etastart = cbind(theta2eta(pnb0,    .lpobs0 , earg = .epobs0 ),
                       theta2eta(mu.init, .lmunb  , earg = .emunb  ),
                       theta2eta(kmat0,   .lsize  , earg = .esize  ))
      etastart = etastart[, interleave.VGAM(ncol(etastart), M = Musual)]
    } # End of if (!length(etastart))


  }), list( .lpobs0 = lpobs0, .lmunb = lmunb, .lsize = lsize,
            .epobs0 = epobs0, .emunb = emunb, .esize = esize,
            .ipobs0 = ipobs0,                 .isize = isize,
            .imethod = imethod, .sinit = shrinkage.init ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Musual <- 3
    NOS <- extra$NOS
    phi0 <- eta2theta(eta[, Musual*(1:NOS)-2], .lpobs0 , earg = .epobs0 )
    munb <- eta2theta(eta[, Musual*(1:NOS)-1], .lmunb  , earg = .emunb  )
    kmat <- eta2theta(eta[, Musual*(1:NOS)  ], .lsize  , earg = .esize  )
    pnb0 <- (kmat / (kmat + munb))^kmat # p(0) from negative binomial
    (1 - phi0) * munb / (1 - pnb0)
  }, list( .lpobs0 = lpobs0, .lsize = lsize, .lmunb = lmunb,
           .epobs0 = epobs0, .emunb = emunb, .esize = esize ))),
  last = eval(substitute(expression({
    misc$link =
      c(rep( .lpobs0 , length = NOS),
        rep( .lmunb  , length = NOS),
        rep( .lsize  , length = NOS))[interleave.VGAM(Musual*NOS,
                                                      M = Musual)]
    temp.names = c(mynames1,
                   mynames2,
                   mynames3)[interleave.VGAM(Musual*NOS, M = Musual)]
    names(misc$link) = temp.names

    misc$earg = vector("list", Musual*NOS)
    names(misc$earg) = temp.names
    for(ii in 1:NOS) {
      misc$earg[[Musual*ii-2]] = .epobs0
      misc$earg[[Musual*ii-1]] = .emunb
      misc$earg[[Musual*ii  ]] = .esize
    }

    misc$nsimEIM = .nsimEIM
    misc$imethod = .imethod
    misc$ipobs0  = .ipobs0
    misc$isize = .isize
  }), list( .lpobs0 = lpobs0, .lmunb = lmunb, .lsize = lsize,
            .epobs0 = epobs0, .emunb = emunb, .esize = esize,
            .ipobs0 = ipobs0, .isize = isize,
            .nsimEIM = nsimEIM,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    NOS = extra$NOS
    Musual <- 3
    phi0 = eta2theta(eta[, Musual*(1:NOS)-2], .lpobs0 , earg = .epobs0 )
    munb = eta2theta(eta[, Musual*(1:NOS)-1], .lmunb  , earg = .emunb  )
    kmat = eta2theta(eta[, Musual*(1:NOS)  ], .lsize  , earg = .esize  )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(w * dzanegbin(x = y, pobs0 = phi0, munb = munb, size = kmat,
                        log = TRUE))
    }
  }, list( .lpobs0 = lpobs0, .lmunb = lmunb, .lsize = lsize,
           .epobs0 = epobs0, .emunb = emunb, .esize = esize ))),
  vfamily = c("zanegbinomial"),
  deriv = eval(substitute(expression({
    Musual <- 3
    NOS = extra$NOS
    y0 = extra$y0

    phi0 = eta2theta(eta[, Musual*(1:NOS)-2, drop = FALSE],
                     .lpobs0 , earg = .epobs0 )
    munb = eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                     .lmunb , earg = .emunb )
    kmat = eta2theta(eta[, Musual*(1:NOS)  , drop = FALSE],
                     .lsize , earg = .esize )
    skip = extra$skip.these


    dphi0.deta = dtheta.deta(phi0, .lpobs0 , earg = .epobs0 )
    dmunb.deta = dtheta.deta(munb, .lmunb  , earg = .emunb  )
    dsize.deta = dtheta.deta(kmat, .lsize  , earg = .esize  )


    tempk = kmat / (kmat + munb)
    tempm = munb / (kmat + munb)
    prob0  = tempk^kmat
    oneminusf0  = 1 - prob0
    df0.dmunb   = -tempk * prob0
    df0.dkmat   = prob0 * (tempm + log(tempk))


    dl.dphi0 = -1 / (1 - phi0)
    dl.dmunb = y / munb - (y + kmat) / (munb + kmat) +
               df0.dmunb / oneminusf0
    dl.dsize = digamma(y + kmat) - digamma(kmat) -
               (y + kmat)/(munb + kmat) + 1 + log(tempk) +
               df0.dkmat / oneminusf0



    dl.dphi0[y == 0] = 1 / phi0[y == 0]  # Do it in one line
    skip = extra$skip.these
    for(spp. in 1:NOS) {
      dl.dsize[skip[, spp.], spp.] =
      dl.dmunb[skip[, spp.], spp.] = 0
    }

    dl.deta23 = c(w) * cbind(dl.dmunb * dmunb.deta,
                             dl.dsize * dsize.deta)


    muphi0 = phi0
    dl.deta1 = if ( .lpobs0 == "logit") {
      c(w) * (y0 - muphi0)
    } else {
      c(w) * dphi0.deta * (y0 / muphi0 - 1) / (1 - muphi0)
    }
    ans = cbind(dl.deta1, dl.deta23)
    ans = ans[, interleave.VGAM(ncol(ans), M = Musual)]
    ans
  }), list( .lpobs0 = lpobs0 , .lmunb = lmunb , .lsize = lsize ,
            .epobs0 = epobs0 , .emunb = emunb , .esize = esize  ))),

  weight = eval(substitute(expression({

    six = dimm(Musual)
    wz =
    run.varcov = matrix(0.0, n, six*NOS-1)
    Musualm1 = Musual - 1







    ind2 = iam(NA, NA, M = Musual - 1, both = TRUE, diag = TRUE)


    for(ii in 1:( .nsimEIM )) {
      ysim = rzanegbin(n = n*NOS, pobs0 = phi0,
                       size = kmat, mu = munb)
      dim(ysim) = c(n, NOS)




      dl.dphi0 = -1 / (1 - phi0)
      dl.dmunb = ysim / munb - (ysim + kmat) / (munb + kmat) +
                 df0.dmunb / oneminusf0
      dl.dsize = digamma(ysim + kmat) - digamma(kmat) -
                 (ysim + kmat)/(munb + kmat) + 1 + log(tempk) +
                 df0.dkmat / oneminusf0




      dl.dphi0[ysim == 0] = 1 / phi0[ysim == 0]  # Do it in one line
      ysim0 = ifelse(ysim == 0, 1, 0)
      skip.sim = matrix(as.logical(ysim0), n, NOS)
      for(spp. in 1:NOS) {
        dl.dsize[skip.sim[, spp.], spp.] =
        dl.dmunb[skip.sim[, spp.], spp.] = 0
      }


      for(kk in 1:NOS) {
        temp2 = cbind(dl.dmunb[, kk] * dmunb.deta[, kk],
                      dl.dsize[, kk] * dsize.deta[, kk])
        small.varcov = temp2[, ind2$row.index] *
                       temp2[, ind2$col.index]




        run.varcov[, ((kk-1)*Musual+2):(kk*Musual)] =
        run.varcov[, ((kk-1)*Musual+2):(kk*Musual)] +
          c(small.varcov[, 1:Musualm1])
        run.varcov[, M + (kk-1)*Musual + 2] =
        run.varcov[, M + (kk-1)*Musual + 2] +
          c(small.varcov[, Musualm1 + 1])
      } # kk; end of NOS
    } # ii; end of nsimEIM


    run.varcov = cbind(run.varcov / .nsimEIM )
    run.varcov = if (intercept.only)
      matrix(colMeans(run.varcov),
             n, ncol(run.varcov), byrow = TRUE) else run.varcov




    wzind1 = sort(c( Musual*(1:NOS) - 1,
                     Musual*(1:NOS) - 0,
                 M + Musual*(1:NOS) - 1))
    wz[, wzind1] = c(w) * run.varcov[, wzind1]




    tmp100 = muphi0 * (1 - muphi0)
    tmp200 = if ( .lpobs0 == "logit") {
      cbind(c(w) * tmp100)
    } else {
      c(w) * cbind(dphi0.deta^2 / tmp100)
    }
    for(ii in 1:NOS) {
      index200 = abs(tmp200[, ii]) < .Machine$double.eps
      if (any(index200)) {
        tmp200[index200, ii] = .Machine$double.eps # Diagonal 0's are bad 
      }
    }
    wz[, Musual*(1:NOS)-2] =  tmp200


    wz
  }), list( .lpobs0 = lpobs0,
            .epobs0 = epobs0,
            .nsimEIM = nsimEIM ))))
} # End of zanegbinomial()










 if (FALSE)
rposnegbin = function(n, munb, size) {
  if (!is.Numeric(size, positive = TRUE))
    stop("argument 'size' must be positive")
  if (!is.Numeric(munb, positive = TRUE))
    stop("argument 'munb' must be positive")
  if (!is.Numeric(n, positive = TRUE, integer.valued = TRUE,
                  allowable.length = 1))
    stop("argument 'n' must be a positive integer")
  ans = rnbinom(n=n, mu = munb, size=size)
  munb = rep(munb, length = n)
  size = rep(size, length = n)
  index = ans == 0
  while(any(index)) {
    more = rnbinom(n=sum(index), mu = munb[index], size=size[index])
    ans[index] = more
    index = ans == 0
  }
  ans
}

 if (FALSE)
dposnegbin = function(x, munb, size, log = FALSE) {
    if (!is.Numeric(size, positive = TRUE))
        stop("argument 'size' must be positive")
    if (!is.Numeric(munb, positive = TRUE))
        stop("argument 'munb' must be positive")
    ans = dnbinom(x = x, mu = munb, size=size, log=log)
    ans0 = dnbinom(x=0, mu = munb, size=size, log = FALSE)
    ans = if (log) ans - log1p(-ans0) else ans/(1-ans0)
    ans[x == 0] = if (log) -Inf else 0
    ans
}











 zipoisson = function(lpstr0 = "logit", llambda = "loge",
                      epstr0 = list(),  elambda = list(),
                      ipstr0 = NULL,    ilambda = NULL,
                      imethod = 1,
                      shrinkage.init = 0.8, zero = NULL)
{

  if (mode(lpstr0) != "character" && mode(lpstr0) != "name")
    lpstr0 = as.character(substitute(lpstr0))
  if (mode(llambda) != "character" && mode(llambda) != "name")
    llambda = as.character(substitute(llambda))

  lpstr00 <- lpstr0
  epstr00 <- epstr0
  ipstr00 <- ipstr0

  if (length(ipstr00))
    if (!is.Numeric(ipstr00, positive = TRUE) ||
        any(ipstr00 >= 1))
      stop("argument 'ipstr0' values must be inside the interval (0,1)")
  if (length(ilambda))
    if (!is.Numeric(ilambda, positive = TRUE))
      stop("argument 'ilambda' values must be positive")

  if (!is.list(epstr00)) epstr00 = list()
  if (!is.list(elambda)) elambda = list()

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
     shrinkage.init < 0 ||
     shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")


  new("vglmff",
  blurb = c("Zero-inflated Poisson\n\n",
            "Links:    ",
            namesof("pstr0",  lpstr00, earg = epstr00 ), ", ",
            namesof("lambda", llambda, earg = elambda ), "\n",
            "Mean:     (1 - pstr0) * lambda"),

  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero)
  }, list( .zero = zero ))),
  initialize = eval(substitute(expression({
    y <- as.matrix(y)
    ncoly <- ncol(y)

    Musual <- 2
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly

    if (any(round(y) != y))
      stop("integer-valued responses only allowed for ",
           "the 'zipoisson' family")

    mynames1 <- paste("pstr0",   if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("lambda",  if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lpstr00 , earg = .epstr00 , tag = FALSE),
          namesof(mynames2, .llambda , earg = .elambda , tag = FALSE))[
          interleave.VGAM(M, M = Musual)]



    if (!length(etastart)) {

      matL <- matrix(if (length( .ilambda )) .ilambda else 0,
                     n, ncoly, byrow = TRUE)
      matP <- matrix(if (length( .ipstr00 )) .ipstr00 else 0,
                     n, ncoly, byrow = TRUE)


      for (spp. in 1:ncoly) {
        yvec <- y[, spp.]

        Phi.init <- 1 - 0.85 * sum(w[yvec > 0]) / sum(w)
        Phi.init[Phi.init <= 0.02] = 0.02  # Last resort
        Phi.init[Phi.init >= 0.98] = 0.98  # Last resort

        if ( length(mustart)) {
          mustart <- matrix(mustart, n, ncoly) # Make sure right size
          Lambda.init <- mustart / (1 - Phi.init)
        } else if ( .imethod == 2) {
          mymean <- weighted.mean(yvec[yvec > 0], w[yvec > 0]) + 1/16
          Lambda.init <- (1 - .sinit) * (yvec + 1/8) + .sinit * mymean
        } else {
          use.this <- median(yvec[yvec > 0]) + 1 / 16
          Lambda.init <- (1 - .sinit) * (yvec + 1/8) + .sinit * use.this
        }

        zipois.Loglikfun <- function(phival, y, x, w, extraargs) {
          sum(w * dzipois(x = y, pstr0 = phival,
                          lambda = extraargs$lambda,
                          log = TRUE))
        }
        phi.grid <- seq(0.02, 0.98, len = 21)
        Phimat.init <- getMaxMin(phi.grid,
                                 objfun = zipois.Loglikfun,
                                 y = y, x = x, w = w,
                                 extraargs = list(lambda = Lambda.init))

        if (length(mustart)) {
          Lambda.init <- Lambda.init / (1 - Phimat.init)
        }

        if (!length( .ipstr00 ))
          matP[, spp.] <- Phimat.init
        if (!length( .ilambda ))
          matL[, spp.] <- Lambda.init
      } # spp.

      etastart <- cbind(theta2eta(matP, .lpstr00, earg = .epstr00 ),
                        theta2eta(matL, .llambda, earg = .elambda ))[,
                        interleave.VGAM(M, M = Musual)]
      mustart <- NULL  # Since etastart has been computed.
    } # End of !length(etastart)
  }), list( .lpstr00 = lpstr00, .llambda = llambda,
            .epstr00 = epstr00, .elambda = elambda,
            .ipstr00 = ipstr00, .ilambda = ilambda,
            .imethod = imethod, .sinit = shrinkage.init ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    phimat = eta2theta(eta[, c(TRUE, FALSE)], .lpstr00 , earg = .epstr00 )
    lambda = eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )
    (1 - phimat) * lambda
  }, list( .lpstr00 = lpstr00, .llambda = llambda,
           .epstr00 = epstr00, .elambda = elambda ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <-
      c(rep( .lpstr00 , length = ncoly),
        rep( .llambda , length = ncoly))[interleave.VGAM(M, M = Musual)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = Musual)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for(ii in 1:ncoly) {
        misc$earg[[Musual*ii-1]] <- .epstr00
        misc$earg[[Musual*ii  ]] <- .elambda
    }

    misc$Musual <- Musual
    misc$imethod <- .imethod
    misc$expected <- TRUE

      misc$pobs0 = phimat + (1 - phimat) * exp(-lambda) # P(Y=0)
      if (length(dimnames(y)[[2]]) > 0)
        dimnames(misc$pobs0) = dimnames(y)
        
  }), list( .lpstr00 = lpstr00, .llambda = llambda,
            .epstr00 = epstr00, .elambda = elambda,
            .imethod = imethod ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    phimat = eta2theta(eta[, c(TRUE, FALSE)], .lpstr00 , earg = .epstr00 )
    lambda = eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(w * dzipois(x = y, pstr0 = phimat, lambda = lambda,
                      log = TRUE))
    }
    }, list( .lpstr00 = lpstr00, .llambda = llambda,
             .epstr00 = epstr00, .elambda = elambda ))),
  vfamily = c("zipoisson"),
  deriv = eval(substitute(expression({
    Musual <- 2
    phimat = eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr00 ,
                       earg = .epstr00 )
    lambda = eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .llambda ,
                       earg = .elambda )

    prob0 = phimat + (1 - phimat) * exp(-lambda)
    index0 = as.matrix(y == 0)

    dl.dphimat = -expm1(-lambda) / prob0
    dl.dphimat[!index0] = -1 / (1 - phimat[!index0])
    dl.dlambda = -(1 - phimat) * exp(-lambda) / prob0
    dl.dlambda[!index0] = (y[!index0] - lambda[!index0]) / lambda[!index0]

    dphimat.deta = dtheta.deta(phimat, .lpstr00 , earg = .epstr00 )
    dlambda.deta = dtheta.deta(lambda, .llambda , earg = .elambda )

    ans = c(w) * cbind(dl.dphimat * dphimat.deta,
                       dl.dlambda * dlambda.deta)
    ans <- ans[, interleave.VGAM(M, M = Musual)]


    if ( .llambda == "loge" && is.empty.list( .elambda ) &&
       any(lambda[!index0] < .Machine$double.eps)) {
      for(spp. in 1:(M / Musual)) {
        ans[!index0[, spp.], Musual * spp.] =
          w[!index0[, spp.]] *
         (y[!index0[, spp.], spp.] - lambda[!index0[, spp.], spp.])
      }
    }

    ans
  }), list( .lpstr00 = lpstr00, .llambda = llambda,
            .epstr00 = epstr00, .elambda = elambda ))),
  weight = eval(substitute(expression({
    wz = matrix(0.0, nrow = n, ncol = M + M-1)

    d2l.dphimat2 = -expm1(-lambda) / ((1 - phimat) * prob0)
    d2l.dlambda2 = (1 - phimat) / lambda -
                   phimat * (1 - phimat) * exp(-lambda) / prob0
    d2l.dphimatlambda = -exp(-lambda) / prob0

    d2l.dphimat2 = as.matrix(d2l.dphimat2)
    d2l.dlambda2 = as.matrix(d2l.dlambda2)
    d2l.dphimatlambda = as.matrix(d2l.dphimatlambda)

    for (ii in 1:(M / Musual)) {
      wz[, iam(Musual * ii - 1, Musual * ii - 1, M)] <-
        d2l.dphimat2[, ii] * dphimat.deta[, ii]^2
      wz[, iam(Musual * ii    , Musual * ii    , M)] <-
        d2l.dlambda2[, ii] * dlambda.deta[, ii]^2
      wz[, iam(Musual * ii - 1, Musual * ii    , M)] <-
        d2l.dphimatlambda[, ii] * dphimat.deta[, ii] * dlambda.deta[, ii]

    }



      c(w) * wz
  }), list( .llambda = llambda, .elambda = elambda ))))
} # zipoisson









 zibinomial = function(lpstr0 = "logit", lprob = "logit",
                       epstr0 = list(),  eprob = list(),
                       ipstr0 = NULL,
                       zero = 1, mv = FALSE, imethod = 1)
{
  if (as.logical(mv))
    stop("argument 'mv' must be FALSE")

  if (mode(lpstr0) != "character" && mode(lpstr0) != "name")
    lpstr0 = as.character(substitute(lpstr0))
  if (mode(lprob) != "character" && mode(lprob) != "name")
    lprob = as.character(substitute(lprob))

  if (is.Numeric(ipstr0))
    if (!is.Numeric(ipstr0, positive = TRUE) || any(ipstr0 >= 1))
      stop("'ipstr0' values must be inside the interval (0,1)")
  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.list(epstr0)) epstr0 = list()
  if (!is.list(eprob ))  eprob = list()


  new("vglmff",
  blurb = c("Zero-inflated binomial\n\n",
            "Links:    ",
            namesof("pstr0", lpstr0, earg = epstr0), ", ",
            namesof("prob" , lprob , earg = eprob ), "\n",
            "Mean:     (1 - pstr0) * prob / (1 - (1 - prob)^w)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
            if (!all(w == 1))
                extra$orig.w = w


    {
        NCOL = function (x)
            if (is.array(x) && length(dim(x)) > 1 ||
            is.data.frame(x)) ncol(x) else as.integer(1)

        if (NCOL(y) == 1) {
            if (is.factor(y)) y <- y != levels(y)[1]
            nn = rep(1, n)
            if (!all(y >= 0 & y <= 1))
                stop("response values must be in [0, 1]")
            if (!length(mustart) && !length(etastart))
                mustart = (0.5 + w * y) / (1.0 + w)


            no.successes = y
            if (min(y) < 0)
                stop("Negative data not allowed!")
            if (any(abs(no.successes - round(no.successes)) > 1.0e-8))
                stop("Number of successes must be integer-valued")

        } else if (NCOL(y) == 2) {
            if (min(y) < 0)
                stop("Negative data not allowed!")
            if (any(abs(y - round(y)) > 1.0e-8))
                stop("Count data must be integer-valued")
            y = round(y)
            nvec = y[, 1] + y[, 2]
            y = ifelse(nvec > 0, y[, 1] / nvec, 0)
            w = w * nvec
            if (!length(mustart) && !length(etastart))
              mustart = (0.5 + nvec * y) / (1 + nvec)
        } else {
            stop("for the binomialff family, response 'y' must be a ",
                 "vector of 0 and 1's\n",
                 "or a factor ",
                 "(first level = fail, other levels = success),\n",
                 "or a 2-column matrix where col 1 is the no. of ",
                 "successes and col 2 is the no. of failures")
        }

    }





    predictors.names =
        c(namesof("pstr0", .lpstr0 , earg = .epstr0 , tag = FALSE),
          namesof("prob" , .lprob  , earg = .eprob  , tag = FALSE))


    phi.init = if (length( .ipstr0 )) .ipstr0 else {
        prob0.est = sum(w[y == 0]) / sum(w)
        if ( .imethod == 1) {
          (prob0.est - (1 - mustart)^w) / (1 - (1 - mustart)^w)
        } else {
          prob0.est
        }
    }

    phi.init[phi.init <= -0.10] = 0.50  # Lots of sample variation
    phi.init[phi.init <=  0.01] = 0.05  # Last resort
    phi.init[phi.init >=  0.99] = 0.95  # Last resort

    if ( length(mustart) && !length(etastart))
      mustart = cbind(rep(phi.init, len = n),
                      mustart) # 1st coln not a real mu
  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob,
            .ipstr0 = ipstr0,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    phi   = eta2theta(eta[, 1], .lpstr0 , earg = .epstr0 )
    mubin = eta2theta(eta[, 2], .lprob  , earg = .eprob  )
    (1 - phi) * mubin
  }, list( .lpstr0 = lpstr0, .lprob = lprob,
           .epstr0 = epstr0, .eprob = eprob ))),
  last = eval(substitute(expression({
    misc$link =    c("pstr0" = .lpstr0 , "prob" = .lprob )
    misc$earg = list("pstr0" = .epstr0 , "prob" = .eprob )
    misc$imethod = .imethod


    if (intercept.only && all(w == w[1])) {
        phi   = eta2theta(eta[1, 1], .lpstr0 , earg = .epstr0 )
        mubin = eta2theta(eta[1, 2], .lprob  , earg = .eprob  )
        misc$pobs0 = phi + (1-phi) * (1-mubin)^w[1] # P(Y=0)
    }
  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob,
            .imethod = imethod ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    cbind(theta2eta(mu[, 1], .lpstr0 , earg = .epstr0 ),
          theta2eta(mu[, 2], .lprob  , earg = .eprob  ))
  }, list( .lpstr0 = lpstr0, .lprob = lprob,
           .epstr0 = epstr0, .eprob = eprob ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    pstr0 = eta2theta(eta[, 1], .lpstr0 , earg = .epstr0 )
    mubin = eta2theta(eta[, 2], .lprob  , earg = .eprob  )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(dzibinom(x = round(w * y), size = w, prob = mubin,
                   log = TRUE, pstr0 = pstr0))
    }
  }, list( .lpstr0 = lpstr0, .lprob = lprob,
           .epstr0 = epstr0, .eprob = eprob ))),
  vfamily = c("zibinomial"),
  deriv = eval(substitute(expression({
    phi   = eta2theta(eta[, 1], .lpstr0 , earg = .epstr0 )
    mubin = eta2theta(eta[, 2], .lprob  , earg = .eprob  )

    prob0 = (1 - mubin)^w    # Actually q^w
    tmp8 = phi + (1 - phi) * prob0
    index = (y == 0)
    dl.dphi = (1 - prob0) / tmp8
    dl.dphi[!index] = -1 / (1 - phi[!index])
    dl.dmubin = -w * (1 - phi) * (1 - mubin)^(w - 1) / tmp8
    dl.dmubin[!index] = w[!index] *
        (y[!index] / mubin[!index] -
        (1 - y[!index]) / (1 - mubin[!index]))
    dphi.deta   = dtheta.deta(phi,   .lpstr0 , earg = .epstr0 )
    dmubin.deta = dtheta.deta(mubin, .lprob  , earg = .eprob  )
    ans = cbind(dl.dphi   * dphi.deta,
                dl.dmubin * dmubin.deta)
      if ( .lprob == "logit") {
          ans[!index,2] = w[!index] * (y[!index] - mubin[!index])
      }
      ans
  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob ))),
  weight = eval(substitute(expression({
    wz = matrix(as.numeric(NA), nrow = n, ncol = dimm(M))



    d2l.dphi2 = (1 - prob0) / ((1 - phi) * tmp8)


    d2l.dphimubin = -w * (1 - mubin)^(w - 1) / tmp8




    d2l.dmubin2 = w * (1 - phi) *
                  (1 / (mubin * (1 - mubin)) -
                   (tmp8 * (w - 1) * (1 - mubin)^(w - 2) -
                    (1 - phi) * w * (1 - mubin)^(2*(w - 1))) / tmp8)


    wz[,iam(1,1,M)] = d2l.dphi2     * dphi.deta^2
    wz[,iam(2,2,M)] = d2l.dmubin2   * dmubin.deta^2
    wz[,iam(1,2,M)] = d2l.dphimubin * dphi.deta * dmubin.deta
    if (TRUE) {
      ind6 = (wz[,iam(2,2,M)] < .Machine$double.eps)
      if (any(ind6))
        wz[ind6,iam(2,2,M)] = .Machine$double.eps
    }
    wz
  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob ))))
}










dzibinom = function(x, size, prob, pstr0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  LLL = max(length(x), length(size), length(prob), length(pstr0))
  if (length(x)     != LLL) x     = rep(x,     len = LLL);
  if (length(size)  != LLL) size  = rep(size,  len = LLL);
  if (length(prob)  != LLL) prob  = rep(prob,  len = LLL);
  if (length(pstr0) != LLL) pstr0 = rep(pstr0, len = LLL);

  ans = dbinom(x = x, size = size, prob = prob, log = TRUE)


  ans = if (log.arg) {
    ifelse(x == 0, log(pstr0 + (1-pstr0) * exp(ans)), log1p(-pstr0) + ans)
  } else {
    ifelse(x == 0,     pstr0 + (1-pstr0) * exp(ans) ,
                    (1-pstr0) * exp(ans))
  }


  prob0 = (1 - prob)^size
  deflat_limit = -prob0 / (1 - prob0)
  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN


  ans
}


pzibinom = function(q, size, prob, pstr0 = 0,
                    lower.tail = TRUE, log.p = FALSE) {

  LLL = max(length(pstr0), length(size), length(prob), length(q))
  if (length(q)      != LLL) q      = rep(q,      len = LLL);
  if (length(size)   != LLL) size   = rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   = rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  = rep(pstr0,  len = LLL);

  ans = pbinom(q, size, prob, lower.tail = lower.tail, log.p = log.p)
  ans = ifelse(q < 0, 0, pstr0 + (1 - pstr0) * ans)


  prob0 = (1 - prob)^size
  deflat_limit = -prob0 / (1 - prob0)
  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN

  ans
}


qzibinom = function(p, size, prob, pstr0 = 0,
                    lower.tail = TRUE, log.p = FALSE) {
  LLL = max(length(p), length(size), length(prob), length(pstr0))
  p     = rep(p,     length = LLL)
  size  = rep(size,  length = LLL)
  prob  = rep(prob,  length = LLL)
  pstr0 = rep(pstr0, length = LLL)


  ans = p 
  ans[p <= pstr0] = 0 
  ans[p >  pstr0] =
    qbinom((p[p > pstr0] - pstr0[p > pstr0]) / (1 - pstr0[p > pstr0]),
           size[p > pstr0],
           prob[p > pstr0],
           lower.tail = lower.tail, log.p = log.p)



  prob0 = (1 - prob)^size
  deflat_limit = -prob0 / (1 - prob0)
  ind0 = (deflat_limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 = pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[p[ind0] <= pobs0] = 0 
    pindex = (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 = pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
    ans[pindex] = qposbinom((p[pindex] - Pobs0) / (1 - Pobs0),
                             size = size[pindex],
                             prob = prob[pindex])
  }

  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN




  ans
}


rzibinom = function(n, size, prob, pstr0 = 0) {
  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE,
                          allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  pstr0 = rep(pstr0, len = use.n)
  size  = rep(size,  len = use.n)
  prob  = rep(prob,  len = use.n)

  ans = rbinom(use.n, size, prob)
  ans[runif(use.n) < pstr0] <- 0



  prob0 = (1 - prob)^size
  deflat_limit = -prob0 / (1 - prob0)
  ind0 = (deflat_limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 = pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] = rposbinom(sum(ind0), size = size[ind0], prob = prob[ind0])
    ans[ind0] = ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN


  ans
}












dzinegbin = function(x, size, prob = NULL, munb = NULL, pstr0 = 0,
                     log = FALSE) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)


  LLL = max(length(pstr0), length(size), length(prob), length(x))
  if (length(x)      != LLL) x      = rep(x,      len = LLL);
  if (length(size)   != LLL) size   = rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   = rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  = rep(pstr0,  len = LLL);


  ans = dnbinom(x = x, size = size, prob = prob, log = log.arg)

  ans = if (log.arg)
    ifelse(x == 0, log(pstr0+(1-pstr0)*exp(ans)), log1p(-pstr0) + ans) else
    ifelse(x == 0,     pstr0+(1-pstr0)*    ans,       (1-pstr0) * ans)



  prob0 = prob^size
  deflat_limit = -prob0 / (1 - prob0)
  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN


  ans
}


pzinegbin = function(q, size, prob = NULL, munb = NULL, pstr0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  LLL = max(length(pstr0), length(size), length(prob), length(q))
  if (length(q)      != LLL) q      = rep(q,      len = LLL);
  if (length(size)   != LLL) size   = rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   = rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  = rep(pstr0,  len = LLL);



  ans = pnbinom(q = q, size = size, prob = prob)
  ans = ifelse(q < 0, 0, pstr0 + (1 - pstr0) * ans)



  prob0 = prob^size
  deflat_limit = -prob0 / (1 - prob0)
  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN


  ans
}


qzinegbin = function(p, size, prob = NULL, munb = NULL, pstr0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size/(size + munb)
  }
  LLL = max(length(p), length(prob), length(pstr0), length(size))
  if (length(p)     != LLL) p      = rep(p,     len = LLL)
  if (length(pstr0) != LLL) pstr0  = rep(pstr0, len = LLL);
  if (length(prob)  != LLL) prob   = rep(prob,  len = LLL)
  if (length(size)  != LLL) size   = rep(size,  len = LLL);

  ans = p 
  ind4 = (p > pstr0)
  ans[!ind4] = 0
  ans[ ind4] = qnbinom(p = (p[ind4] - pstr0[ind4]) / (1 - pstr0[ind4]),
                       size = size[ind4], prob = prob[ind4])



  prob0 = prob^size
  deflat_limit = -prob0 / (1 - prob0)
  ind0 = (deflat_limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 = pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[p[ind0] <= pobs0] = 0 
    pindex = (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 = pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
    ans[pindex] = qposnegbin((p[pindex] - Pobs0) / (1 - Pobs0),
                              size = size[pindex],
                              prob = prob[pindex])
  }


  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN



  ans
}


rzinegbin = function(n, size, prob = NULL, munb = NULL, pstr0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE,
                          allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n


  pstr0 = rep(pstr0, len = use.n)
  size  = rep(size,  len = use.n)
  prob  = rep(prob,  len = use.n)


  ans = rnbinom(n = use.n, size = size, prob = prob)
  ans = ifelse(runif(use.n) < pstr0, rep(0, use.n), ans)



  prob0 = rep(prob^size, len = use.n)
  deflat_limit = -prob0 / (1 - prob0)
  ind0 = (deflat_limit <= pstr0) & (pstr0 <  0)
  if (any(ind0, na.rm = TRUE)) {
    pobs0 = pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] = rposnegbin(sum(ind0, na.rm = TRUE), size = size[ind0],
                    prob = prob[ind0])
    ans[ind0] = ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN

  ans
}







zinegbinomial.control <- function(save.weight = TRUE, ...)
{
  list(save.weight = save.weight)
}


 zinegbinomial =
  function(lpstr0 = "logit", lmunb = "loge", lsize = "loge",
           epstr0 = list(),  emunb = list(), esize = list(),
           ipstr0 = NULL,                    isize = NULL,
           zero = c(-1, -3),
           imethod = 1, shrinkage.init = 0.95,
           nsimEIM = 250)
{

  if (mode(lpstr0) != "character" && mode(lpstr0) != "name")
    lpstr0 = as.character(substitute(lpstr0))
  if (mode(lmunb) != "character" && mode(lmunb) != "name")
    lmunb = as.character(substitute(lmunb))
  if (mode(lsize) != "character" && mode(lsize) != "name")
    lsize = as.character(substitute(lsize))


  if (length(ipstr0) &&
     (!is.Numeric(ipstr0, positive = TRUE) ||
      any(ipstr0 >= 1)))
    stop("argument 'ipstr0' must contain values in (0,1)")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("argument 'isize' must contain positive values only")

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1, 2 or 3")

  if (!is.Numeric(nsimEIM, allowable.length = 1, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 50)
    warning("argument 'nsimEIM' should be greater than 50, say")

  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
      shrinkage.init < 0 ||
      shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")

  if (!is.list(epstr0)) epstr0 = list()
  if (!is.list(emunb))  emunb  = list()
  if (!is.list(esize))  esize  = list()

  new("vglmff",
  blurb = c("Zero-inflated negative binomial\n\n",
            "Links:    ",
            namesof("pstr0", lpstr0, earg = epstr0, tag = FALSE), ", ",
            namesof("munb",  lmunb,  earg = emunb,  tag = FALSE), ", ",
            namesof("size",  lsize,  earg = esize,  tag = FALSE), "\n",
            "Mean:     (1 - pstr0) * munb"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 3
    eval(negzero.expression)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    Musual <- 3
    y <- as.matrix(y)
    extra$NOS = NOS = ncoly = ncol(y)  # Number of species
    if (length(dimnames(y)))
      extra$dimnamesy2 = dimnames(y)[[2]]

    mynames1 = if (NOS == 1) "pstr0" else paste("pstr0", 1:NOS, sep = "")
    mynames2 = if (NOS == 1) "munb"  else paste("munb",  1:NOS, sep = "")
    mynames3 = if (NOS == 1) "size"  else paste("size",  1:NOS, sep = "")
    predictors.names =
      c(namesof(mynames1, .lpstr0 , earg = .epstr0 , tag = FALSE),
        namesof(mynames2, .lmunb  , earg = .emunb  , tag = FALSE),
        namesof(mynames3, .lsize  , earg = .esize  , tag = FALSE))[
        interleave.VGAM(Musual*NOS, M = Musual)]

    if (!length(etastart)) {
      mum.init = if ( .imethod == 3) {
        y + 1/16
      } else {
        mum.init = y
        for(iii in 1:ncol(y)) {
          index = (y[, iii] > 0)
          mum.init[, iii] = if ( .imethod == 2)
              weighted.mean(y[index, iii], w     = w[index]) else
                 median(rep(y[index, iii], times = w[index])) + 1/8
        }
        (1 - .sinit) * (y + 1/16) + .sinit * mum.init
      }


      pstr0.init = if (length( .ipstr0 )) {
        matrix( .ipstr0 , n, ncoly, byrow = TRUE)
      } else {
        pstr0.init = y
        for(iii in 1:ncol(y))
          pstr0.init[, iii] = sum(w[y[, iii] == 0]) / sum(w)
        pstr0.init[pstr0.init <= 0.02] = 0.02  # Last resort
        pstr0.init[pstr0.init >= 0.98] = 0.98  # Last resort
        pstr0.init
      }

        kay.init =
        if ( is.Numeric( .isize )) {
          matrix( .isize, nrow = n, ncol = ncoly, byrow = TRUE)
        } else {
          zinegbin.Loglikfun = function(kval, y, x, w, extraargs) {
            index0 = (y == 0)
            pstr0vec = extraargs$pstr0
            muvec = extraargs$mu


            ans1 = 0.0
            if (any( index0))
              ans1 = ans1 + sum(w[ index0] *
                     dzinegbin(x = y[ index0], size = kval,
                               munb = muvec[ index0],
                               pstr0 = pstr0vec[ index0], log = TRUE))
            if (any(!index0))
              ans1 = ans1 + sum(w[!index0] *
                     dzinegbin(x = y[!index0], size = kval,
                               munb = muvec[!index0],
                               pstr0 = pstr0vec[!index0], log = TRUE))
            ans1
          }
          k.grid = 2^((-6):6)
          kay.init = matrix(0, nrow = n, ncol = NOS)
          for(spp. in 1:NOS) {
            kay.init[, spp.] = getMaxMin(k.grid,
                              objfun = zinegbin.Loglikfun,
                              y = y[, spp.], x = x, w = w,
                              extraargs = list(pstr0 = pstr0.init[, spp.],
                                               mu  = mum.init[, spp.]))
          }
          kay.init
        }

        etastart = cbind(theta2eta(pstr0.init, .lpstr0 , earg = .epstr0 ),
                         theta2eta(mum.init,   .lmunb  , earg = .emunb  ),
                         theta2eta(kay.init,   .lsize  , earg = .esize  ))
        etastart =
          etastart[, interleave.VGAM(ncol(etastart), M = Musual)]
    }
  }), list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
            .epstr0 = epstr0, .emunb = emunb, .esize = esize,
            .ipstr0 = ipstr0,                 .isize = isize,
            .sinit = shrinkage.init,
            .imethod = imethod ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Musual = 3
    NOS = extra$NOS
    pstr0 = eta2theta(eta[, Musual*(1:NOS)-2, drop = FALSE],
                      .lpstr0 , earg = .epstr0 )
    munb  = eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                      .lmunb  , earg = .emunb  )
    fv.matrix = (1 - pstr0) * munb
    if (length(extra$dimnamesy2))
      dimnames(fv.matrix) = list(dimnames(pstr0)[[1]], extra$dimnamesy2)
    fv.matrix
  }, list( .lpstr0 = lpstr0, .lsize = lsize, .lmunb = lmunb,
           .epstr0 = epstr0, .esize = esize, .emunb = emunb ))),
  last = eval(substitute(expression({
    misc$link =
      c(rep( .lpstr0 , length = NOS),
        rep( .lmunb  , length = NOS),
        rep( .lsize  , length = NOS))[interleave.VGAM(Musual*NOS,
                                                      M = Musual)]
    temp.names =
      c(mynames1,
        mynames2,
        mynames3)[interleave.VGAM(Musual*NOS, M = Musual)]
    names(misc$link) = temp.names

    misc$earg = vector("list", Musual*NOS)
    names(misc$earg) = temp.names
    for(ii in 1:NOS) {
      misc$earg[[Musual*ii-2]] = .epstr0
      misc$earg[[Musual*ii-1]] = .emunb
      misc$earg[[Musual*ii  ]] = .esize
    }

    misc$imethod = .imethod
    misc$nsimEIM = .nsimEIM
    misc$expected = TRUE
    misc$Musual = Musual
    misc$ipstr0  = .ipstr0
    misc$isize = .isize
    if (intercept.only) {
   pstr0.val = eta2theta(eta[1,Musual*(1:NOS)-2], .lpstr0 , earg= .epstr0 )
   munb.val  = eta2theta(eta[1,Musual*(1:NOS)-1], .lmunb  , earg= .emunb  )
   kval      = eta2theta(eta[1,Musual*(1:NOS)  ], .lsize  , earg= .esize  )
   misc$pobs0 =      pstr0.val +
                (1 - pstr0.val) * (kval / (kval + munb.val))^kval # P(Y=0)
    }
  }), list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
            .epstr0 = epstr0, .emunb = emunb, .esize = esize,
            .ipstr0 = ipstr0,                 .isize = isize,
            .nsimEIM = nsimEIM, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Musual <- 3
    NOS = extra$NOS
    pstr0 = eta2theta(eta[, Musual*(1:NOS)-2, drop = FALSE],
                      .lpstr0 , earg = .epstr0 )
    munb  = eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                      .lmunb , earg = .emunb )
    kmat  = eta2theta(eta[, Musual*(1:NOS)  , drop = FALSE],
                      .lsize , earg = .esize )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(w * dzinegbin(x = y, size = kmat, munb = munb,
                        pstr0 = pstr0, log = TRUE))
    }
  }, list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
           .epstr0 = epstr0, .emunb = emunb, .esize = esize ))),
  vfamily = c("zinegbinomial"),
  deriv = eval(substitute(expression({
    Musual <- 3
    NOS = extra$NOS

    pstr0 = eta2theta(eta[, Musual*(1:NOS)-2, drop = FALSE],
                      .lpstr0 , earg = .epstr0 )
    munb  = eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                      .lmunb  , earg = .emunb  )
    kmat  = eta2theta(eta[, Musual*(1:NOS)  , drop = FALSE],
                      .lsize  , earg = .esize  )

    dpstr0.deta = dtheta.deta(pstr0, .lpstr0 , earg = .epstr0 )
    dmunb.deta  = dtheta.deta(munb , .lmunb  , earg = .emunb  )
    dsize.deta  = dtheta.deta(kmat , .lsize  , earg = .esize  )
    dthetas.detas =
        (cbind(dpstr0.deta,
               dmunb.deta,
               dsize.deta))[, interleave.VGAM(Musual*NOS, M = Musual)]



    dl.dpstr0 = -1 / (1 - pstr0)
    dl.dmunb = y / munb - (y + kmat) / (munb + kmat)
    dl.dsize = digamma(y + kmat) - digamma(kmat) -
               (y + kmat) / (munb + kmat) + 1 +
               log(kmat / (kmat + munb))



    for(spp. in 1:NOS) {
      index0 = (y[, spp.] == 0)
      if (!any(index0) || !any(!index0))
        stop("must have some 0s AND some positive counts in the data")

      kmat.  =  kmat[index0, spp.]
      munb.  =  munb[index0, spp.]
      pstr0. = pstr0[index0, spp.]


      tempk. = kmat. / (kmat. + munb.)
      tempm. = munb. / (kmat. + munb.)
      prob0. = tempk.^kmat.
      df0.dmunb.  = -tempk.* prob0.
      df0.dkmat.  = prob0. * (tempm. + log(tempk.))


      denom. = pstr0. + (1 - pstr0.) * prob0.
     dl.dpstr0[index0, spp.]  = (1 - prob0.) / denom.
      dl.dmunb[index0, spp.]  = (1 - pstr0.) * df0.dmunb. / denom.
      dl.dsize[index0, spp.]  = (1 - pstr0.) * df0.dkmat. / denom.
    } # of spp.


    dl.dthetas =
      cbind(dl.dpstr0,
            dl.dmunb,
            dl.dsize)[, interleave.VGAM(Musual*NOS, M = Musual)]


      c(w) * dl.dthetas * dthetas.detas
  }), list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
            .epstr0 = epstr0, .emunb = emunb, .esize = esize ))),

  weight = eval(substitute(expression({



    wz = matrix(0, n, Musual*M - Musual)

    ind3 = iam(NA, NA, M = Musual, both = TRUE, diag = TRUE)

    run.varcov = array(0.0, c(n, length(ind3$row.index), NOS))

    for(ii in 1:( .nsimEIM )) {
      ysim = rzinegbin(n = n*NOS, pstr0 = pstr0,
                       size = kmat, mu = munb)
      dim(ysim) = c(n, NOS)
      index0 = (ysim[, spp.] == 0)

      dl.dpstr0 = -1 / (1 - pstr0)
      dl.dmunb = ysim / munb - (ysim + kmat) / (munb + kmat)
      dl.dsize = digamma(ysim + kmat) - digamma(kmat) -
                 (ysim + kmat) / (munb + kmat) + 1 +
                 log(kmat / (kmat + munb))


      for(spp. in 1:NOS) {
        index0 = (ysim[, spp.] == 0)
        if (!any(index0) || !any(!index0))
          stop("must have some 0s AND some positive counts in the data")

        kmat.  =  kmat[index0, spp.]
        munb.  =  munb[index0, spp.]
        pstr0. = pstr0[index0, spp.]


        tempk. = kmat. / (kmat. + munb.)
        tempm. = munb. / (kmat. + munb.)
        prob0.  = tempk.^kmat.
        df0.dmunb.  = -tempk.* prob0.
        df0.dkmat.  = prob0. * (tempm. + log(tempk.))


        denom. = pstr0. + (1 - pstr0.) * prob0.
       dl.dpstr0[index0, spp.] = (1 - prob0.) / denom.
        dl.dmunb[index0, spp.] = (1 - pstr0.) * df0.dmunb. / denom.
        dl.dsize[index0, spp.] = (1 - pstr0.) * df0.dkmat. / denom.


        sdl.dthetas = cbind(dl.dpstr0[, spp.],
                             dl.dmunb[, spp.],
                             dl.dsize[, spp.])

        temp3 = sdl.dthetas
        run.varcov[,, spp.] = run.varcov[,, spp.] +
                              temp3[, ind3$row.index] *
                              temp3[, ind3$col.index]


      } # End of for(spp.) loop
    } # End of ii nsimEIM loop

    run.varcov = run.varcov / .nsimEIM

    wz1 = if (intercept.only) {
      for(spp. in 1:NOS) {
        for(jay in 1:length(ind3$row.index)) {
          run.varcov[, jay, spp.] = mean(run.varcov[, jay, spp.])
        }
      }
      run.varcov
    } else {
      run.varcov
    }

    for(spp. in 1:NOS) {
      wz1[,, spp.] = wz1[,, spp.] *
                     dthetas.detas[, Musual * (spp. - 1) + ind3$row] *
                     dthetas.detas[, Musual * (spp. - 1) + ind3$col]
    }

    for(spp. in 1:NOS) {
      for(jay in 1:Musual) {
        for(kay in jay:Musual) {
          cptr = iam((spp. - 1) * Musual + jay,
                     (spp. - 1) * Musual + kay, M = M)
          temp.wz1 = wz1[,, spp.]
          wz[, cptr] = temp.wz1[, iam(jay, kay, M = Musual)]
        }
      }
    }
    c(w) * wz
  }), list( .lpstr0 = lpstr0,
            .epstr0 = epstr0, .nsimEIM = nsimEIM ))))
} # End of zinegbinomial








 zipoissonff <- function(llambda = "loge", lprobp = "logit",
                         elambda = list(), eprobp = list(),
                         ilambda = NULL,   iprobp = NULL, imethod = 1,
                         shrinkage.init = 0.8, zero = -2)
{
  lprobp. <- lprobp
  eprobp. <- eprobp
  iprobp. <- iprobp

  if (mode(llambda) != "character" && mode(llambda) != "name")
      llambda <- as.character(substitute(llambda))
  if (mode(lprobp.) != "character" && mode(lprobp.) != "name")
      lprobp. <- as.character(substitute(lprobp.))

  if (length(ilambda))
    if (!is.Numeric(ilambda, positive = TRUE))
      stop("'ilambda' values must be positive")
  if (length(iprobp.))
    if (!is.Numeric(iprobp., positive = TRUE) ||
      any(iprobp. >= 1))
      stop("'iprobp' values must be inside the interval (0,1)")

  if (!is.list(elambda)) elambda <- list()
  if (!is.list(eprobp.)) eprobp. <- list()

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
    imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.Numeric(shrinkage.init, allowable.length = 1) ||
    shrinkage.init < 0 ||
    shrinkage.init > 1)
    stop("bad input for argument 'shrinkage.init'")

  new("vglmff",
  blurb = c("Zero-inflated Poisson\n\n",
            "Links:    ",
            namesof("lambda", llambda, earg = elambda), ", ",
            namesof("probp",  lprobp., earg = eprobp.), "\n",
            "Mean:     probp * lambda"),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero)
  }, list( .zero = zero ))),
  initialize = eval(substitute(expression({
    y <- as.matrix(y)

    ncoly <- ncol(y)
    Musual <- 2
    extra$ncoly <- ncoly
    extra$Musual <- Musual
    M <- Musual * ncoly

    if (any(round(y) != y))
      stop("responses must be integer-valued")

    mynames1 <- paste("lambda", if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("probp",  if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .llambda, earg = .elambda, tag = FALSE),
          namesof(mynames2, .lprobp., earg = .eprobp., tag = FALSE))[
          interleave.VGAM(M, M = Musual)]


      if (!length(etastart)) {

        matL <- matrix(if (length( .ilambda )) .ilambda else 0,
                       n, ncoly, byrow = TRUE)
        matP <- matrix(if (length( .iprobp. )) .iprobp. else 0,
                       n, ncoly, byrow = TRUE)

        for (jay in 1:ncoly) {
          yjay <- y[, jay]

          Phi0.init <- 1 - 0.85 * sum(w[yjay > 0]) / sum(w)
          Phi0.init[Phi0.init <= 0.02] = 0.02  # Last resort
          Phi0.init[Phi0.init >= 0.98] = 0.98  # Last resort

          if ( length(mustart)) {
            mustart <- matrix(mustart, n, ncoly) # Make sure right size
            Lambda.init <- mustart / (1 - Phi0.init)
          } else if ( .imethod == 2) {
            mymean <- weighted.mean(yjay[yjay > 0], w[yjay > 0]) + 1/16
            Lambda.init <- (1 - .sinit) * (yjay + 1/8) + .sinit * mymean
          } else {
            use.this <- median(yjay[yjay > 0]) + 1 / 16
            Lambda.init <- (1 - .sinit) * (yjay + 1/8) + .sinit * use.this
          }

          zipois.Loglikfun <- function(phival, y, x, w, extraargs) {
            sum(w * dzipois(x = y, pstr0 = phival,
                            lambda = extraargs$lambda,
                            log = TRUE))
          }
          phi0.grid <- seq(0.02, 0.98, len = 21)
          Phi0mat.init <- getMaxMin(phi0.grid,
                                   objfun = zipois.Loglikfun,
                                   y = y, x = x, w = w,
                                   extraargs = list(lambda = Lambda.init))
          if (length(mustart)) {
            Lambda.init <- Lambda.init / (1 - Phi0mat.init)
          }

        if (!length( .ilambda ))
          matL[, jay] <- Lambda.init
        if (!length( .iprobp. ))
          matP[, jay] <- Phi0mat.init
      }

      etastart <- cbind(theta2eta(    matL, .llambda , earg = .elambda ),
                        theta2eta(1 - matP, .lprobp. , earg = .eprobp. ))[,
                        interleave.VGAM(M, M = Musual)]

      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .lprobp. = lprobp., .llambda = llambda,
            .eprobp. = eprobp., .elambda = elambda,
            .iprobp. = iprobp., .ilambda = ilambda,
            .imethod = imethod, .sinit = shrinkage.init ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Musual <- 2
    ncoly <- extra$ncoly
    lambda <- eta2theta(eta[, Musual*(1:ncoly) - 1], .llambda,
                        earg = .elambda )
    probp. <- eta2theta(eta[, Musual*(1:ncoly)    ], .lprobp.,
                        earg = .eprobp. )
    probp. * lambda
  }, list( .lprobp. = lprobp., .llambda = llambda,
           .eprobp. = eprobp., .elambda = elambda ))),
  last = eval(substitute(expression({
    Musual <- extra$Musual
    misc$link <-
      c(rep( .llambda, length = ncoly),
        rep( .lprobp., length = ncoly))[interleave.VGAM(M, M = Musual)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = Musual)]
    names(misc$link) <- temp.names


    misc$earg <- vector("list", Musual * ncoly)
    names(misc$earg) <- temp.names
    for(ii in 1:ncoly) {
        misc$earg[[Musual*ii-1]] <- .elambda
        misc$earg[[Musual*ii  ]] <- .eprobp.
    }

    misc$Musual <- Musual
    misc$imethod <- .imethod
    misc$expected = TRUE

      misc$pobs0 <- (1 - probp.) + probp. * exp(-lambda) # P(Y=0)
      misc$pobs0 <- as.matrix(misc$pobs0)
      if (length(dimnames(y)[[2]]) > 0)
        dimnames(misc$pobs0) = dimnames(y)
  }), list( .lprobp. = lprobp., .llambda = llambda,
            .eprobp. = eprobp., .elambda = elambda,
            .imethod = imethod ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Musual <- 2
    ncoly <- extra$ncoly
    lambda <- eta2theta(eta[, Musual*(1:ncoly) - 1], .llambda,
                        earg = .elambda )
    probp. <- eta2theta(eta[, Musual*(1:ncoly)    ], .lprobp.,
                        earg = .eprobp. )


    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(w * dzipois(x = y, pstr0 = 1 - probp., lambda = lambda,
                      log = TRUE))
    }
  }, list( .lprobp. = lprobp., .llambda = llambda,
           .eprobp. = eprobp., .elambda = elambda ))),
  vfamily = c("zipoissonff"),
  deriv = eval(substitute(expression({
    Musual <- 2
    ncoly <- extra$ncoly
    lambda <- eta2theta(eta[, Musual*(1:ncoly) - 1], .llambda,
                        earg = .elambda )
    probp. <- eta2theta(eta[, Musual*(1:ncoly)    ], .lprobp.,
                        earg = .eprobp. )


    dlambda.deta <- dtheta.deta(lambda, .llambda, earg = .elambda )
    dprobp..deta <- dtheta.deta(probp., .lprobp., earg = .eprobp. )

    denom <- 1 + probp. * expm1(-lambda)
    ind0 <- (y == 0)
    dl.dlambda <- -probp. * exp(-lambda) / denom
    dl.dlambda[!ind0] <- (y[!ind0] - lambda[!ind0]) / lambda[!ind0]
    dl.dprobp. <- expm1(-lambda) / denom
    dl.dprobp.[!ind0] <- 1 / probp.[!ind0]

    ans <- c(w) * cbind(dl.dlambda * dlambda.deta,
                        dl.dprobp. * dprobp..deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M = Musual)]


    if ( .llambda == "loge" && is.empty.list( .elambda ) &&
       any(lambda[!ind0] < .Machine$double.eps)) {
      for(spp. in 1:ncoly) {
        ans[!ind0[, spp.], Musual * spp.] =
          w[!ind0[, spp.]] *
         (y[!ind0[, spp.], spp.] - lambda[!ind0[, spp.], spp.])
      }
    }



    ans
  }), list( .lprobp. = lprobp., .llambda = llambda,
            .eprobp. = eprobp., .elambda = elambda ))),
  weight = eval(substitute(expression({


    wz <- matrix(0, nrow = n, ncol = M + M-1)
    d2l.dlambda2 <-  (    probp.) / lambda -
                    probp. * (1 - probp.) * exp(-lambda) / denom
    d2l.dprobp.2 <- -expm1(-lambda) / ((  probp.) * denom)
    d2l.dphilambda <- +exp(-lambda) / denom


    if (ncoly == 1) {  # Make sure these are matrices
      d2l.dlambda2 <- cbind(d2l.dlambda2)
      d2l.dprobp.2 <- cbind(d2l.dprobp.2)
      dlambda.deta <- cbind(dlambda.deta)
      dprobp..deta <- cbind(dprobp..deta)
      d2l.dphilambda <- cbind(d2l.dphilambda)
    }

    for (ii in 1:ncoly) {
      wz[, iam(Musual*ii - 1, Musual*ii - 1, M)] <-
        d2l.dlambda2[, ii] *
        dlambda.deta[, ii]^2
      wz[, iam(Musual*ii    , Musual*ii    , M)] <-
        d2l.dprobp.2[, ii] *
        dprobp..deta[, ii]^2
      wz[, iam(Musual*ii - 1, Musual*ii    , M)] <-
       d2l.dphilambda[, ii] *
         dprobp..deta[, ii] *
         dlambda.deta[, ii]



    } # ii


    c(w) * wz
  }), list( .llambda = llambda ))))
}







dzigeom = function(x, prob, pstr0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  LLL = max(length(x), length(prob), length(pstr0))
  if (length(x)      != LLL) x      = rep(x,      len = LLL);
  if (length(prob)   != LLL) prob   = rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  = rep(pstr0,  len = LLL);


  ans = dgeom(x = x, prob = prob, log = TRUE)


  ans = if (log.arg) {
    ifelse(x == 0, log(pstr0 + (1 - pstr0) * exp(ans)),
                   log1p(-pstr0) + ans)
  } else {
    ifelse(x == 0,     pstr0 + (1 - pstr0) * exp(ans) ,
                               (1 - pstr0) * exp(ans))
  }



  prob0 = prob
  deflat_limit = -prob0 / (1 - prob0)
  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN

  ans
}



pzigeom = function(q, prob, pstr0 = 0) {


  LLL = max(length(q), length(prob), length(pstr0))
  if (length(q)      != LLL) q      = rep(q,      len = LLL);
  if (length(prob)   != LLL) prob   = rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  = rep(pstr0,  len = LLL);

  ans = pgeom(q, prob)
  ans = ifelse(q < 0, 0, pstr0 + (1-pstr0) * ans)


  prob0 = prob
  deflat_limit = -prob0 / (1 - prob0)
  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN

  ans
}



qzigeom = function(p, prob, pstr0 = 0) {
  LLL = max(length(p), length(prob), length(pstr0))
  ans = p = rep(p,     len = LLL)
  prob    = rep(prob,  len = LLL)
  pstr0   = rep(pstr0, len = LLL)
  ans[p <= pstr0] = 0 
  ind1 = (p > pstr0)
  ans[ind1] =
    qgeom((p[ind1] - pstr0[ind1]) / (1 - pstr0[ind1]),
          prob = prob[ind1])


  prob0 = prob
  deflat_limit = -prob0 / (1 - prob0)
  ind0 = (deflat_limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 = pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[p[ind0] <= pobs0] = 0 
    pindex = (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 = pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
    ans[pindex] = 1 + qgeom((p[pindex] - Pobs0) / (1 - Pobs0),
                            prob = prob[pindex])
  }

  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN

  ans
}



rzigeom = function(n, prob, pstr0 = 0) {
  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE,
                          allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n


  pstr0 = rep(pstr0, len = use.n)
  prob  = rep(prob,  len = use.n)


  ans = rgeom(use.n, prob)
  ans[runif(use.n) < pstr0] = 0


  prob0 = prob
  deflat_limit = -prob0 / (1 - prob0)
  ind0 = (deflat_limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 = pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] = 1 + rgeom(sum(ind0), prob = prob[ind0])
    ans[ind0] = ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat_limit] = NaN
  ans[pstr0 > 1] = NaN


  ans
}





 zigeometric = function(lprob = "logit", eprob = list(),
                        lpstr0  = "logit", epstr0  = list(),
                        iprob = NULL,    ipstr0  = NULL,
                        imethod = 1,
                        bias.red = 0.5,
                        zero = 2)
{


  expected = TRUE

  if (mode(lprob) != "character" && mode(lprob) != "name")
    lprob = as.character(substitute(lprob))
  if (mode(lpstr0) != "character" && mode(lpstr0) != "name")
    lpstr0 = as.character(substitute(lpstr0))

  if (!is.list(eprob))    eprob    = list()
  if (!is.list(epstr0))  epstr0  = list()


  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")
  if (length(ipstr0))
    if (!is.Numeric(ipstr0, positive = TRUE) ||
        ipstr0 >= 1)
      stop("argument 'ipstr0' is out of range")

  if (!is.Numeric(bias.red, allowable.length = 1, positive = TRUE) ||
     bias.red > 1)
    stop("argument 'bias.red' must be between 0 and 1")


  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Zero-inflated geometric distribution,\n",
            "P[Y = 0] = pstr0 + (1 - pstr0) * prob,\n",
            "P[Y = y] = (1 - pstr0) * prob * (1 - prob)^y, ",
            "y = 1, 2, ...\n\n",
            "Link:     ",
            namesof("prob",   lprob,   earg = eprob ), ", ",
            namesof("pstr0",  lpstr0,  earg = epstr0), "\n",
            "Mean:     (1 - pstr0) * (1 - prob) / prob"),
  constraints = eval(substitute(expression({
      constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero)
  }, list( .zero = zero ))),
  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a 1-column matrix")

    if (any(y < 0))
      stop("all responses must be >= 0")
    if (any(y != round(y)))
      stop("response should be integer-valued")

    predictors.names =
            c(namesof("prob", .lprob, earg = .earg, tag = FALSE),
              namesof("pstr0",  .lpstr0,  earg = .epstr0, tag = FALSE))

    if (!length(etastart)) {
      prob.init = if ( .imethod == 3)
                      .bias.red / (1 + y + 1/8) else
                  if ( .imethod == 2)
                      .bias.red / (1 +    mean(y) + 1/8) else
                      .bias.red / (1 + weighted.mean(y, w)  + 1/8)
      prob.init = if (length( .iprob )) {
        rep( .iprob, len = n)
      } else {
        rep(prob.init, len = n)
      }


      prob0.est = sum(w[y == 0]) / sum(w)
      psze.init = if ( .imethod == 3)
                      prob0.est / 2 else
                  if ( .imethod == 1)
                      max(0.05, (prob0.est - median(prob.init))) else
                      prob0.est / 5
      psze.init = if (length( .ipstr0 )) {
        rep( .ipstr0 , len = n)
      } else {
        rep( psze.init, len = n)
      }



      etastart =
        cbind(theta2eta(prob.init, .lprob, earg = .eprob),
              theta2eta(psze.init, .lpstr0,  earg = .epstr0))

    }
  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .eprob = eprob, .epstr0 = epstr0,
            .iprob = iprob, .ipstr0 = ipstr0,
            .bias.red = bias.red,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob    = eta2theta(eta[, 1], .lprob,    earg = .eprob)
    pstr0  = eta2theta(eta[, 2], .lpstr0 , earg = .epstr0 )
    (1 - pstr0) * (1 - prob) / prob
  }, list( .lprob = lprob, .lpstr0 = lpstr0,
           .eprob = eprob, .epstr0 = epstr0 ))),
  last = eval(substitute(expression({
    misc$link =    c(prob = .lprob, pstr0 = .lpstr0 )
    misc$earg = list(prob = .eprob, pstr0 = .epstr0 )
    misc$imethod = .imethod
    misc$zero = .zero
    misc$bias.red = .bias.red
    misc$expected = .expected
    misc$ipstr0 = .ipstr0


  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .eprob = eprob, .epstr0 = epstr0,
                            .ipstr0 = ipstr0,
            .zero = zero,
            .expected = expected,
            .bias.red = bias.red,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    prob = eta2theta(eta[, 1], .lprob, earg = .eprob)
    pstr0  = eta2theta(eta[, 2], .lpstr0 , earg = .epstr0 )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(w * dzigeom(x = y, prob = prob, pstr0 = pstr0, log = TRUE))
    }
  }, list( .lprob = lprob, .lpstr0 = lpstr0,
           .eprob = eprob, .epstr0 = epstr0 ))),
  vfamily = c("zigeometric"),

  deriv = eval(substitute(expression({
    prob = eta2theta(eta[, 1], .lprob, earg = .eprob)
    pstr0  = eta2theta(eta[, 2], .lpstr0 , earg = .epstr0 )


    prob0 = prob  # P(Y == 0)
    tmp8 = pstr0 + (1 - pstr0) * prob0
    index0 = (y == 0)

    dl.dpstr0 = (1 - prob0) / tmp8
    dl.dpstr0[!index0] = -1 / (1 - pstr0[!index0])

    dl.dprob = (1 - pstr0) / tmp8
    dl.dprob[!index0]   = 1 / prob[!index0] -
                          y[!index0] / (1 - prob[!index0])

    dprob.deta = dtheta.deta(prob, .lprob, earg = .eprob )
    dpstr0.deta  = dtheta.deta(pstr0 , .lpstr0 , earg = .epstr0 )

    dl.deta12 = 
    c(w) * cbind(dl.dprob * dprob.deta,
                 dl.dpstr0  *  dpstr0.deta)
    dl.deta12
  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .eprob = eprob, .epstr0 = epstr0 ))),
  weight = eval(substitute(expression({
    ed2l.dprob2 = (1 - pstr0) * (1 / (prob^2 * (1 - prob)) +
                              (1 - pstr0) / tmp8)
    ed2l.dpstr0.prob = 1 / tmp8
    ed2l.dpstr02 = (1 - prob0) / ((1 - pstr0) * tmp8)

    od2l.dprob2 = ((1 - pstr0) / tmp8)^2
    od2l.dprob2[!index0] = 1 / (prob[!index0])^2 +
                           y[!index0] / (1 - prob[!index0])^2
    od2l.dpstr0.prob = (tmp8 + (1 - prob0) * (1 - pstr0)) / tmp8^2
    od2l.dpstr0.prob[!index0] = 0


    od2l.dpstr02 = ((1 - prob0) / tmp8)^2
    od2l.dpstr02[!index0] = 1 / (1 - pstr0[!index0])^2


    wz = matrix(as.numeric(NA), nrow = n, ncol = dimm(M))
    if ( .expected ) {
      wz[,iam(1,1,M)] = ed2l.dprob2       * dprob.deta^2
      wz[,iam(2,2,M)] = ed2l.dpstr02     * dpstr0.deta^2
      wz[,iam(1,2,M)] = ed2l.dpstr0.prob * dprob.deta * dpstr0.deta
    } else {
      wz[,iam(1,1,M)] = od2l.dprob2       * dprob.deta^2
      wz[,iam(2,2,M)] = od2l.dpstr02     * dpstr0.deta^2
      wz[,iam(1,2,M)] = od2l.dpstr0.prob * dprob.deta * dpstr0.deta
    }


    c(w) * wz
  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .expected = expected,
            .eprob = eprob, .epstr0 = epstr0 ))))
}






dzageom = function(x, prob, pobs0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  LLL = max(length(x), length(prob), length(pobs0))
  if (length(x)      != LLL) x      = rep(x,       len = LLL);
  if (length(prob)   != LLL) prob   = rep(prob,    len = LLL);
  if (length(pobs0)  != LLL) pobs0    = rep(pobs0, len = LLL);
  ans = rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  index0 = (x == 0)

  if (log.arg) {
    ans[ index0] = log(pobs0[index0])
    ans[!index0] = log1p(-pobs0[!index0]) +
                   dposgeom(x[!index0],
                            prob = prob[!index0], log = TRUE)
  } else {
    ans[ index0] = pobs0[index0]
    ans[!index0] = (1-pobs0[!index0]) *
                   dposgeom(x[!index0],
                            prob = prob[!index0])
  }
  ans
}



pzageom = function(q, prob, pobs0 = 0) {

  LLL = max(length(q), length(prob), length(pobs0))
  if (length(q)      != LLL) q      = rep(q,      len = LLL);
  if (length(prob)   != LLL) prob   = rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  = rep(pobs0,  len = LLL);
  ans = rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans[q >  0] = pobs0[q > 0] +
                (1 - pobs0[q > 0]) *
                pposgeom(q[q > 0], prob = prob[q > 0])
  ans[q <  0] = 0
  ans[q == 0] = pobs0[q == 0]
  ans
}


qzageom = function(p, prob, pobs0 = 0) {

  LLL = max(length(p), length(prob), length(pobs0))
  if (length(p)      != LLL) p      = rep(p,      len = LLL);
  if (length(prob)   != LLL) prob   = rep(prob,   len = LLL);
  if (length(pobs0)    != LLL) pobs0    = rep(pobs0,    len = LLL);

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans = p
  ind4 = (p > pobs0)
  ans[!ind4] = 0.0
  ans[ ind4] = qposgeom((p[ind4] - pobs0[ind4]) / (1 - pobs0[ind4]),
                         prob = prob[ind4])
  ans
}


rzageom = function(n, prob, pobs0 = 0) {
  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE,
                          allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  ans = rposgeom(use.n, prob)
  if (length(pobs0) != use.n) pobs0 = rep(pobs0, len = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ifelse(runif(use.n) < pobs0, 0, ans)
}










dzabinom = function(x, size, prob, pobs0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log))
    stop("bad input for argument 'log'")
  rm(log)

  LLL = max(length(x), length(size), length(prob), length(pobs0))
  if (length(x)      != LLL) x      = rep(x,      len = LLL);
  if (length(size)   != LLL) size   = rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   = rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  = rep(pobs0,  len = LLL);
  ans = rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  index0 = (x == 0)

  if (log.arg) {
    ans[ index0] = log(pobs0[index0])
    ans[!index0] = log1p(-pobs0[!index0]) +
                   dposbinom(x[!index0], size = size[!index0],
                             prob = prob[!index0], log = TRUE)
  } else {
    ans[ index0] = pobs0[index0]
    ans[!index0] = (1-pobs0[!index0]) *
                   dposbinom(x[!index0], size = size[!index0],
                             prob = prob[!index0])
  }
  ans
}



pzabinom = function(q, size, prob, pobs0 = 0) {

  LLL = max(length(q), length(size), length(prob), length(pobs0))
  if (length(q)      != LLL) q      = rep(q,      len = LLL);
  if (length(size)   != LLL) size   = rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   = rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  = rep(pobs0,  len = LLL);
  ans = rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) ||
      any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans[q >  0] = pobs0[q > 0] +
                (1 - pobs0[q > 0]) *
                pposbinom(q[q > 0], size = size[q > 0], prob = prob[q > 0])
  ans[q <  0] = 0
  ans[q == 0] = pobs0[q == 0]
  ans
}


qzabinom = function(p, size, prob, pobs0 = 0) {

  LLL = max(length(p), length(size), length(prob), length(pobs0))
  if (length(p)      != LLL) p      = rep(p,      len = LLL);
  if (length(size)   != LLL) size   = rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   = rep(prob,   len = LLL);
  if (length(pobs0)    != LLL) pobs0    = rep(pobs0,    len = LLL);

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans = p
  ind4 = (p > pobs0)
  ans[!ind4] = 0.0
  ans[ ind4] = qposbinom((p[ind4] - pobs0[ind4]) / (1 - pobs0[ind4]),
                         size = size[ind4],
                         prob = prob[ind4])
  ans
}


rzabinom = function(n, size, prob, pobs0 = 0) {
  use.n = if ((length.n <- length(n)) > 1) length.n else
          if (!is.Numeric(n, integer.valued = TRUE,
                          allowable.length = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  ans = rposbinom(use.n, size, prob)
  if (length(pobs0) != use.n) pobs0 = rep(pobs0, len = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ifelse(runif(use.n) < pobs0, 0, ans)
}





 zabinomial = function(lprob  = "logit", eprob  = list(),
                       lpobs0 = "logit", epobs0 = list(),
                       iprob = NULL, ipobs0 = NULL,
                       imethod = 1,
                       zero = 2)
{




  if (mode(lprob) != "character" && mode(lprob) != "name")
    lprob = as.character(substitute(lprob))
  if (mode(lpobs0) != "character" && mode(lpobs0) != "name")
    lpobs0 = as.character(substitute(lpobs0))

  if (!is.list(eprob))   eprob   = list()
  if (!is.list(epobs0))  epobs0  = list()

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")
  if (length(ipobs0))
    if (!is.Numeric(ipobs0, positive = TRUE) ||
        ipobs0 >= 1)
      stop("argument 'ipobs0' is out of range")



  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Zero-altered binomial distribution ",
            "(Bernoulli and positive-binomial conditional model)\n\n",
            "P[Y = 0] = pobs0,\n",
            "P[Y = y] = (1 - pobs0) * dposbinom(x = y, size, prob), ",
            "y = 1, 2, ..., size,\n\n",
            "Link:     ",
            namesof("prob" ,   lprob,  earg = eprob), ", ",
            namesof("pobs0",   lpobs0, earg = epobs0), "\n",
            "Mean:     (1 - pobs0) * prob / (1 - (1 - prob)^size)"),
  constraints = eval(substitute(expression({
      constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(Musual = 2,
         zero = .zero)
  }, list( .zero = zero ))),
  initialize = eval(substitute(expression({
            if (!all(w == 1))
                extra$orig.w = w



    {
        NCOL = function (x)
            if (is.array(x) && length(dim(x)) > 1 ||
            is.data.frame(x)) ncol(x) else as.integer(1)

        if (NCOL(y) == 1) {
            if (is.factor(y)) y <- y != levels(y)[1]
            nn = rep(1, n)
            if (!all(y >= 0 & y <= 1))
                stop("response values must be in [0, 1]")
            if (!length(mustart) && !length(etastart))
                mustart = (0.5 + w * y) / (1.0 + w)


            no.successes = y
            if (min(y) < 0)
                stop("Negative data not allowed!")
            if (any(abs(no.successes - round(no.successes)) > 1.0e-8))
                stop("Number of successes must be integer-valued")

        } else if (NCOL(y) == 2) {
            if (min(y) < 0)
                stop("Negative data not allowed!")
            if (any(abs(y - round(y)) > 1.0e-8))
                stop("Count data must be integer-valued")
            y = round(y)
            nvec = y[, 1] + y[, 2]
            y = ifelse(nvec > 0, y[, 1] / nvec, 0)
            w = w * nvec
            if (!length(mustart) && !length(etastart))
              mustart = (0.5 + nvec * y) / (1 + nvec)
        } else {
            stop("for the binomialff family, response 'y' must be a ",
                 "vector of 0 and 1's\n",
                 "or a factor ",
                 "(first level = fail, other levels = success),\n",
                 "or a 2-column matrix where col 1 is the no. of ",
                 "successes and col 2 is the no. of failures")
        }

    }
    if (!all(w == 1))
      extra$new.w = w


    y = as.matrix(y)
    extra$y0 = y0 = ifelse(y == 0, 1, 0)
    extra$NOS = NOS = ncoly = ncol(y)  # Number of species
    extra$skip.these = skip.these = matrix(as.logical(y0), n, NOS)



    predictors.names =
        c(namesof("prob" , .lprob  , earg = .eprob  , tag = FALSE),
          namesof("pobs0", .lpobs0 , earg = .epobs0 , tag = FALSE))


    orig.w = if (length(extra$orig.w)) extra$orig.w else 1
    new.w  = if (length(extra$new.w))  extra$new.w  else 1
    Size = new.w / orig.w

    phi.init = if (length( .ipobs0 )) .ipobs0 else {
        prob0.est = sum(Size[y == 0]) / sum(Size)
        if ( .imethod == 1) {
          (prob0.est - (1 - mustart)^Size) / (1 - (1 - mustart)^Size)
        } else
        if ( .imethod == 2) {
          prob0.est
        } else {
          prob0.est * 0.5
        }
    }

    phi.init[phi.init <= -0.10] = 0.50  # Lots of sample variation
    phi.init[phi.init <=  0.01] = 0.05  # Last resort
    phi.init[phi.init >=  0.99] = 0.95  # Last resort




    if (!length(etastart)) {
      etastart =
        cbind(theta2eta( mustart, .lprob,  earg = .eprob  ),
              theta2eta(phi.init, .lpobs0, earg = .epobs0 ))

      mustart <- NULL
    }
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0,
            .iprob = iprob, .ipobs0 = ipobs0,
            .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob  = eta2theta(eta[, 1], .lprob,  earg = .eprob  )
    phi0  = eta2theta(eta[, 2], .lpobs0, earg = .epobs0 )
    orig.w = if (length(extra$orig.w)) extra$orig.w else 1
    new.w  = if (length(extra$new.w))  extra$new.w  else 1
    Size = new.w / orig.w
    (1 - phi0) * prob / (1 - (1 - prob)^Size)
  }, list( .lprob = lprob, .lpobs0 = lpobs0,
           .eprob = eprob, .epobs0 = epobs0 ))),

  last = eval(substitute(expression({
    misc$link =    c(prob = .lprob, pobs0 = .lpobs0 )
    misc$earg = list(prob = .eprob, pobs0 = .epobs0 )
    misc$imethod = .imethod
    misc$zero = .zero
    misc$expected = TRUE
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0,
            .zero = zero,
            .imethod = imethod ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    orig.w = if (length(extra$orig.w)) extra$orig.w else 1
    new.w  = if (length(extra$new.w))  extra$new.w  else 1
    Size = new.w / orig.w
    prob  = eta2theta(eta[, 1], .lprob  , earg = .eprob  )
    pobs0 = eta2theta(eta[, 2], .lpobs0 , earg = .epobs0 )
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(orig.w * dzabinom(x = round(y * Size), size = Size,
                            prob = prob, pobs0 = pobs0,
                            log = TRUE))
    }
  }, list( .lprob = lprob, .lpobs0 = lpobs0,
           .eprob = eprob, .epobs0 = epobs0 ))),
  vfamily = c("zabinomial"),

  deriv = eval(substitute(expression({
    NOS = if (length(extra$NOS)) extra$NOS else 1
    Musual = 2

    orig.w = if (length(extra$orig.w)) extra$orig.w else 1
    new.w  = if (length(extra$new.w))  extra$new.w  else 1
    Size = new.w / orig.w

    prob = eta2theta(eta[, 1], .lprob  , earg = .eprob  )
    phi0 = eta2theta(eta[, 2], .lpobs0 , earg = .epobs0 )

    dprob.deta = dtheta.deta(prob, .lprob , earg = .eprob  )
    dphi0.deta = dtheta.deta(phi0, .lpobs0, earg = .epobs0 )

    df0.dprob   = -Size *              (1 -  prob)^(Size - 1)
    df02.dprob2 =  Size * (Size - 1) * (1 -  prob)^(Size - 2)
    prob0  = (1 -  prob)^(Size)
    oneminusf0  = 1 - prob0


    dl.dprob =  c(w)      * (y / prob - (1 - y) / (1 - prob)) +
                c(orig.w) * df0.dprob / oneminusf0
    dl.dphi0 = -1 / (1 - phi0)


    dl.dphi0[y == 0] = 1 / phi0[y == 0]  # Do it in one line
    skip = extra$skip.these
    for(spp. in 1:NOS) {
      dl.dprob[skip[, spp.], spp.] = 0
    }


    ans <- cbind(            dl.dprob * dprob.deta,
                 c(orig.w) * dl.dphi0 * dphi0.deta)
                 
    ans
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0 ))),


  weight = eval(substitute(expression({
    wz = matrix(0.0, n, Musual)

    usualmeanY =  prob
    meanY = (1 - phi0) * usualmeanY / oneminusf0


    term1 =  c(Size) * (meanY /      prob^2 -
                        meanY / (1 - prob)^2) +
             c(Size) * (1 - phi0) / (1 - prob)^2

    term2 =  -(1 - phi0) * df02.dprob2 / oneminusf0
    term3 =  -(1 - phi0) * (df0.dprob  / oneminusf0)^2
    ed2l.dprob2 = term1 + term2 + term3
    wz[, iam(1,1,M)] = ed2l.dprob2 * dprob.deta^2


    mu.phi0 = phi0
    tmp100 = mu.phi0 * (1.0 - mu.phi0)
    tmp200 = if ( .lpobs0 == "logit" && is.empty.list( .epobs0 )) {
      tmp100
    } else {
      (dphi0.deta^2) / tmp100
    }
    wz[, iam(2,2,M)] = tmp200

    c(orig.w) * wz
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0 ))))
}





 zageometric = function(lpobs0 = "logit", lprob = "logit",
                        epobs0 = list(),  eprob = list(),
                        imethod = 1,
                        ipobs0 = NULL, iprob = NULL,
                        zero = NULL) {


  if (mode(lpobs0) != "character" && mode(lpobs0) != "name")
    lpobs0 = as.character(substitute(lpobs0))
  if (mode(lprob) != "character" && mode(lprob) != "name")
    lprob = as.character(substitute(lprob))

  if (!is.list(epobs0)) epobs0 = list()
  if (!is.list(eprob))  eprob  = list()

  if (!is.Numeric(imethod, allowable.length = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")
  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
       max(iprob) >= 1)
    stop("argument 'iprob' out of range")
  if (length(ipobs0))
    if (!is.Numeric(ipobs0, positive = TRUE) ||
       max(ipobs0) >= 1)
      stop("argument 'ipobs0' out of range")


  new("vglmff",
  blurb = c("Zero-altered geometric ",
            "(Bernoulli and positive-geometric conditional model)\n\n",
            "Links:    ",
         namesof("pobs0", lpobs0, earg = epobs0, tag = FALSE), ", ",
         namesof("prob" , lprob , earg = eprob , tag = FALSE), "\n",
            "Mean:     (1 - pobs0) / prob"),

  constraints = eval(substitute(expression({

    dotzero <- .zero
    Musual <- 2
    eval(negzero.expression)
  }), list( .zero = zero ))),
  initialize = eval(substitute(expression({
    Musual <- 2
    y <- as.matrix(y)
    if (any(y != round(y )))
      stop("the response must be integer-valued")
    if (any(y < 0))
      stop("the response must not have negative values")

    extra$y0 = y0 = ifelse(y == 0, 1, 0)
    extra$NOS = NOS = ncoly = ncol(y)  # Number of species
    extra$skip.these = skip.these = matrix(as.logical(y0), n, NOS)

    mynames1 = if (ncoly == 1) "pobs0"  else paste("pobs0",  1:ncoly, sep = "")
    mynames2 = if (ncoly == 1) "prob" else paste("prob", 1:ncoly, sep = "")
    predictors.names = 
        c(namesof(mynames1, .lpobs0 , earg = .epobs0 , tag = FALSE),
          namesof(mynames2, .lprob  , earg = .eprob  , tag = FALSE))[
          interleave.VGAM(Musual*NOS, M = Musual)]

    if (!length(etastart)) {

      foo = function(x) mean(as.numeric(x == 0))
      phi0.init = matrix(apply(y, 2, foo), n, ncoly, byrow = TRUE)
      if (length( .ipobs0 ))
        phi0.init = matrix( .ipobs0 , n, ncoly, byrow = TRUE)


      prob.init =
        if ( .imethod == 2)
          1 / (1 + y + 1/16) else
        if ( .imethod == 1)
          (1 - phi0.init) / (1 + matrix(apply(y, 2, weighted.mean, w = w),
                                    n, ncoly, byrow = TRUE) + 1/16) else
          (1 - phi0.init) / (1 + matrix(apply(y, 2, median),
                                    n, ncoly, byrow = TRUE) + 1/16)
      if (length( .iprob ))
        prob.init = matrix( .iprob , n, ncoly, byrow = TRUE)



      etastart = cbind(theta2eta(phi0.init, .lpobs0 , earg = .epobs0 ),
                       theta2eta(prob.init, .lprob , earg = .eprob ))
      etastart = etastart[, interleave.VGAM(ncol(etastart), M = Musual)]
    }
  }), list( .lpobs0 = lpobs0, .lprob = lprob,
            .epobs0 = epobs0, .eprob = eprob,
            .ipobs0 = ipobs0, .iprob = iprob,
            .imethod = imethod ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS = extra$NOS
    Musual <- 2

    phi0 = cbind(eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                             .lpobs0 , earg = .epobs0 ))
    prob = cbind(eta2theta(eta[, Musual*(1:NOS)-0, drop = FALSE],
                             .lprob  , earg = .eprob ))

    (1 - phi0) / prob
  }, list( .lpobs0 = lpobs0, .lprob = lprob,
           .epobs0 = epobs0, .eprob = eprob ))),
  last = eval(substitute(expression({
    temp.names = c(rep( .lpobs0 , len = NOS),
                   rep( .lprob  , len = NOS))
    temp.names = temp.names[interleave.VGAM(Musual*NOS, M = Musual)]
    misc$link  = temp.names
    misc$expected = TRUE
    misc$earg = vector("list", Musual * NOS)
    misc$imethod = .imethod
    misc$ipobs0  = .ipobs0
    misc$iprob   = .iprob

    names(misc$link) <-
    names(misc$earg) <-
        c(mynames1, mynames2)[interleave.VGAM(Musual*NOS, M = Musual)]

    for(ii in 1:NOS) {
      misc$earg[[Musual*ii-1]] = .epobs0
      misc$earg[[Musual*ii  ]] = .eprob
    }
  }), list( .lpobs0 = lpobs0, .lprob = lprob,
            .epobs0 = epobs0, .eprob = eprob,
            .ipobs0 = ipobs0, .iprob = iprob,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    NOS = extra$NOS
    Musual <- 2

    phi0 = cbind(eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                           .lpobs0 , earg = .epobs0 ))
    prob = cbind(eta2theta(eta[, Musual*(1:NOS)-0, drop = FALSE],
                           .lprob  , earg = .eprob  ))

    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(w * dzageom(x = y, pobs0 = phi0, prob = prob, log = TRUE))
    }
  }, list( .lpobs0 = lpobs0, .lprob = lprob,
           .epobs0 = epobs0, .eprob = eprob ))),
  vfamily = c("zageometric"),
  deriv = eval(substitute(expression({
    Musual <- 2
    NOS = extra$NOS
    y0 = extra$y0
    skip = extra$skip.these

    phi0 = cbind(eta2theta(eta[, Musual*(1:NOS)-1, drop = FALSE],
                           .lpobs0 , earg = .epobs0 ))
    prob = cbind(eta2theta(eta[, Musual*(1:NOS)-0, drop = FALSE],
                           .lprob  , earg = .eprob  ))


    dl.dprob =  1 / prob - (y - 1) / (1 - prob)
    dl.dphi0 = -1 / (1 - phi0)


    for(spp. in 1:NOS) {
      dl.dphi0[skip[, spp.], spp.] = 1 / phi0[skip[, spp.], spp.]
      dl.dprob[skip[, spp.], spp.] = 0
    }
    dphi0.deta = dtheta.deta(phi0, .lpobs0 , earg = .epobs0 )
    dprob.deta = dtheta.deta(prob, .lprob  , earg = .eprob  )


    ans <- c(w) * cbind(dl.dphi0 * dphi0.deta,
                        dl.dprob * dprob.deta)
    ans = ans[, interleave.VGAM(ncol(ans), M = Musual)]
    ans
  }), list( .lpobs0 = lpobs0, .lprob = lprob,
            .epobs0 = epobs0, .eprob = eprob ))),
  weight = eval(substitute(expression({

    wz = matrix(0.0, n, Musual*NOS)


    ed2l.dprob2 = (1 - phi0) / (prob^2 * (1 - prob))

    wz[, NOS+(1:NOS)] = c(w) * ed2l.dprob2 * dprob.deta^2


    mu.phi0 = phi0
    tmp100 = mu.phi0 * (1.0 - mu.phi0)
    tmp200 = if ( .lpobs0 == "logit" && is.empty.list( .epobs0 )) {
      cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (dphi0.deta^2) / tmp100)
    }
    wz[, 1:NOS] =  tmp200


    wz = wz[, interleave.VGAM(ncol(wz), M = Musual)]
    wz
  }), list( .lpobs0 = lpobs0,
            .epobs0 = epobs0 ))))
} #   End of zageometric







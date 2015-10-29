# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.












dzanegbin <- function(x, size, prob = NULL, munb = NULL, pobs0 = 0,
                      log = FALSE) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(pobs0), length(prob), length(size))
  if (length(x)     != LLL) x     <- rep(x,     len = LLL)
  if (length(pobs0) != LLL) pobs0 <- rep(pobs0, len = LLL)
  if (length(prob)  != LLL) prob  <- rep(prob,  len = LLL)
  if (length(size)  != LLL) size  <- rep(size,  len = LLL)

  ans <- rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  if (!is.Numeric(prob, positive = TRUE))
    stop("argument 'prob' must be in (0,Inf)")
  if (!is.Numeric(size, positive = TRUE))
    stop("argument 'size' must be in (0,Inf)")
  index0 <- x == 0

  if (log.arg) {
    ans[ index0] <- log(pobs0[index0])
    ans[!index0] <- log1p(-pobs0[!index0]) +
                    dposnegbin(x[!index0], prob = prob[!index0],
                               size = size[!index0], log = TRUE)
  } else {
    ans[ index0] <- pobs0[index0]
    ans[!index0] <- (1 - pobs0[!index0]) * dposnegbin(x[!index0],
                      prob = prob[!index0], size = size[!index0])
  }
  ans
}



pzanegbin <- function(q, size, prob = NULL, munb = NULL, pobs0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  LLL <- max(length(q), length(pobs0), length(prob), length(size))
  if (length(q)     != LLL) q     <- rep(q,     len = LLL)
  if (length(pobs0) != LLL) pobs0 <- rep(pobs0, len = LLL)
  if (length(prob)  != LLL) prob  <- rep(prob,  len = LLL)
  if (length(size)  != LLL) size  <- rep(size,  len = LLL)
  ans <- rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  qindex <- (q >  0)
  ans[ qindex] <- pobs0[qindex] + (1 - pobs0[qindex]) *
                  pposnegbin(q[qindex], size = size[qindex],
                                        prob = prob[qindex])
  ans[q <  0] <- 0
  ans[q == 0] <- pobs0[q == 0]
  ans
}


qzanegbin <- function(p, size, prob = NULL, munb = NULL, pobs0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size/(size + munb)
  }

  LLL <- max(length(p), length(pobs0), length(prob), length(size))
  if (length(p)     != LLL) p      <- rep(p,     len = LLL)
  if (length(pobs0) != LLL) pobs0  <- rep(pobs0, len = LLL)
  if (length(prob)  != LLL) prob   <- rep(prob,  len = LLL)
  if (length(size)  != LLL) size   <- rep(size,  len = LLL)
  ans <- rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ans <- p
  ans[p <= pobs0] <- 0
  pindex <- (p > pobs0)
  ans[pindex] <-
    qposnegbin((p[pindex] - pobs0[pindex]) / (1 - pobs0[pindex]),
               prob = prob[pindex],
               size = size[pindex])
  ans
}


rzanegbin <- function(n, size, prob = NULL, munb = NULL, pobs0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  ans <- rposnegbin(n = use.n, prob = prob, size = size)
  if (length(pobs0) != use.n)
    pobs0 <- rep(pobs0, len = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")

  ifelse(runif(use.n) < pobs0, 0, ans)
}





dzapois <- function(x, lambda, pobs0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(lambda), length(pobs0))
  if (length(x)      != LLL) x      <- rep(x,      len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL)
  ans <- rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  index0 <- (x == 0)

  if (log.arg) {
    ans[ index0] <- log(pobs0[index0])
    ans[!index0] <- log1p(-pobs0[!index0]) +
                    dpospois(x[!index0], lambda[!index0], log = TRUE)
  } else {
    ans[ index0] <- pobs0[index0]
    ans[!index0] <- (1 - pobs0[!index0]) *
                    dpospois(x[!index0], lambda[!index0])
  }
  ans
}



pzapois <- function(q, lambda, pobs0 = 0) {
  LLL <- max(length(q), length(lambda), length(pobs0))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL)
  ans <- rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  ans[q >  0] <-    pobs0[q > 0] +
                 (1-pobs0[q > 0]) * ppospois(q[q > 0], lambda[q > 0])
  ans[q <  0] <- 0
  ans[q == 0] <- pobs0[q == 0]
  ans
}


qzapois <- function(p, lambda, pobs0 = 0) {
  LLL <- max(length(p), length(lambda), length(pobs0))
  if (length(p)      != LLL) p      <- rep(p,      len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ans <- p
  ind4 <- (p > pobs0)
  ans[!ind4] <- 0
  ans[ ind4] <- qpospois((p[ind4] - pobs0[ind4]) / (1 - pobs0[ind4]),
                         lambda = lambda[ind4])
  ans
}


rzapois <- function(n, lambda, pobs0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  ans <- rpospois(use.n, lambda)
  if (length(pobs0) != use.n)
    pobs0 <- rep(pobs0, length = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must in [0,1]")

  ifelse(runif(use.n) < pobs0, 0, ans)
}





dzipois <- function(x, lambda, pstr0 = 0, log = FALSE) {



  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(lambda), length(pstr0))
  if (length(x)      != LLL) x      <- rep(x,      len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL)

  ans <- x + lambda + pstr0



  index0 <- (x == 0)
  if (log.arg) {
    ans[ index0] <- log(pstr0[ index0] + (1 - pstr0[ index0]) *
                       dpois(x[ index0], lambda[ index0]))
    ans[!index0] <- log1p(-pstr0[!index0]) +
                   dpois(x[!index0], lambda[!index0], log = TRUE)
  } else {
    ans[ index0] <-      pstr0[ index0] + (1 - pstr0[ index0]) *
                       dpois(x[ index0], lambda[ index0])
    ans[!index0] <- (1 - pstr0[!index0]) *
                    dpois(x[!index0], lambda[!index0])
  }


  deflat.limit <- -1 / expm1(lambda)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}


pzipois <- function(q, lambda, pstr0 = 0) {

  LLL <- max(length(pstr0), length(lambda), length(q))
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(q)      != LLL) q      <- rep(q,      len = LLL)

  ans <- ppois(q, lambda)
  ans <- ifelse(q < 0, 0, pstr0 + (1 - pstr0) * ans)


  deflat.limit <- -1 / expm1(lambda)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans
}


qzipois <- function(p, lambda, pstr0 = 0) {

  LLL <- max(length(p), length(lambda), length(pstr0))
  if (length(p)      != LLL) p      <- rep(p,      len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL)
  ans    <- p

  ans[p <= pstr0] <- 0 
  pindex <- (p > pstr0)
  ans[pindex] <-
    qpois((p[pindex] - pstr0[pindex]) / (1 - pstr0[pindex]),
          lambda = lambda[pindex])


  deflat.limit <- -1 / expm1(lambda)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * exp(-lambda[ind0])
    ans[p[ind0] <= pobs0] <- 0 
    pindex <- (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * exp(-lambda[pindex])
    ans[pindex] <- qpospois((p[pindex] - Pobs0) / (1 - Pobs0),
                            lambda = lambda[pindex])
  }


  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans
}


rzipois <- function(n, lambda, pstr0 = 0) {

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (length(pstr0)  != use.n) pstr0  <- rep(pstr0,  len = use.n)
  if (length(lambda) != use.n) lambda <- rep(lambda, len = use.n)
 
  ans <- rpois(use.n, lambda)
  ans <- ifelse(runif(use.n) < pstr0, 0, ans)



  prob0 <- exp(-lambda)
  deflat.limit <- -1 / expm1(lambda)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] <- rpospois(sum(ind0), lambda[ind0]) 
    ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}






 yip88 <- function(link = "loge", n.arg = NULL) {








  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Zero-inflated Poisson (based on Yip (1988))\n\n",
            "Link:     ",
            namesof("lambda", link, earg), "\n",
            "Variance: (1 - pstr0) * lambda"),
  first = eval(substitute(expression({
    zero <- y == 0
    if (any(zero)) {
      if (length(extra)) extra$sumw <- sum(w) else
        extra <- list(sumw=sum(w))
      if (is.numeric(.n.arg) && extra$sumw != .n.arg) 
        stop("value of 'n.arg' conflicts with data ",
             "(it need not be specified anyway)")
      warning("trimming out the zero observations")


      axa.save <-  attr(x, "assign")
      x <- x[!zero,, drop = FALSE]
      attr(x, "assign") <- axa.save    # Don't lose these!!
      w <- w[!zero]
      y <- y[!zero]
    } else {
      if (!is.numeric(.n.arg)) 
        stop("n.arg must be supplied")
    }
        
  }), list( .n.arg = n.arg ))),

  initialize = eval(substitute(expression({
    narg <- if (is.numeric(.n.arg)) .n.arg else extra$sumw
    if (sum(w) > narg)
      stop("sum(w) > narg")

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)


    predictors.names <-
      namesof("lambda", .link, list(theta = NULL), tag = FALSE)

    if (!length(etastart)) {
      lambda.init <- rep(median(y), length = length(y))
      etastart <- theta2eta(lambda.init, .link , earg = .earg )
    }
    if (length(extra)) {
      extra$sumw <- sum(w)
      extra$narg <- narg   # For @linkinv
    } else {
      extra <- list(sumw = sum(w), narg = narg)
    }
  }), list( .link = link, .earg = earg, .n.arg = n.arg ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    lambda <- eta2theta(eta, .link, .earg)
    temp5 <- exp(-lambda)
    pstr0 <- (1 - temp5 - extra$sumw/extra$narg) / (1 - temp5)
    if (any(pstr0 <= 0))
      stop("non-positive value(s) of pstr0")
    (1 - pstr0) * lambda
  }, list( .link = link, .earg = earg ))),

  last = eval(substitute(expression({
    misc$link <-    c(lambda = .link )

    misc$earg <- list(lambda = .earg )

    if (intercept.only) {
      suma <- extra$sumw
      pstr0 <- (1 - temp5[1] - suma / narg) / (1 - temp5[1])
      pstr0 <- if (pstr0 < 0 || pstr0 > 1) NA else pstr0
      misc$pstr0 <- pstr0
    }
  }), list( .link = link, .earg = earg ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    lambda <- eta2theta(eta, .link)
    temp5 <- exp(-lambda)
    pstr0 <- (1 - temp5 - extra$sumw / extra$narg) / (1 - temp5)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
                 dzipois(x = y, pstr0 = pstr0, lambda = lambda, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),

  vfamily = c("yip88"),
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta, .link , earg = .earg )
    temp5 <- exp(-lambda)
    dl.dlambda <- -1 + y/lambda - temp5/(1-temp5)
    dlambda.deta <- dtheta.deta(lambda, .link , earg = .earg )
    w * dl.dlambda * dlambda.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    d2lambda.deta2 <- d2theta.deta2(lambda, .link , earg = .earg )
    d2l.dlambda2 <- -y / lambda^2 + temp5 / (1 - temp5)^2
    -w * (d2l.dlambda2*dlambda.deta^2 + dl.dlambda*d2lambda.deta2)
  }), list( .link = link, .earg = earg ))))
}





 zapoisson <-
  function(lpobs0 = "logit", llambda = "loge",
           type.fitted = c("mean", "pobs0", "onempobs0"),
           zero = NULL) {




  lpobs.0 <- as.list(substitute(lpobs0))
  epobs.0 <- link2list(lpobs.0)
  lpobs.0 <- attr(epobs.0, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "onempobs0"))[1]



  new("vglmff",
  blurb = c("Zero-altered Poisson ",
            "(Bernoulli and positive-Poisson conditional model)\n\n",
            "Links:    ",
            namesof("pobs0",  lpobs.0, earg = epobs.0, tag = FALSE), ", ",
            namesof("lambda", llambda, earg = elambda, tag = FALSE), "\n",
            "Mean:     (1 - pobs0) * lambda / (1 - exp(-lambda))"),

  constraints = eval(substitute(expression({

    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 2
    if (any(y < 0))
      stop("the response must not have negative values")

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)
    extra$dimnamesy <- dimnames(y)
    extra$type.fitted      <- .type.fitted

    mynames1 <- if (ncoly == 1) "pobs0"    else
                paste("pobs0",    1:ncoly, sep = "")
    mynames2 <- if (ncoly == 1) "lambda" else
                paste("lambda", 1:ncoly, sep = "")
    predictors.names <-
        c(namesof(mynames1, .lpobs.0, earg = .epobs.0, tag = FALSE),
          namesof(mynames2, .llambda, earg = .elambda, tag = FALSE))[
          interleave.VGAM(M1*NOS, M = M1)]

    if (!length(etastart)) {
      etastart <-
        cbind(theta2eta((0.5 + w*y0) / (1+w),
                        .lpobs.0, earg = .epobs.0 ),
              matrix(1, n, NOS))  # 1 here is any old value
      for (spp. in 1:NOS) {
        sthese <- skip.these[, spp.]
        etastart[!sthese, NOS+spp.] =
          theta2eta(y[!sthese, spp.] / (-expm1(-y[!sthese, spp.])),
                    .llambda, earg = .elambda )
      }
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = M1)]
    }
  }), list( .lpobs.0 = lpobs.0, .llambda = llambda,
            .epobs.0 = epobs.0, .elambda = elambda,
            .type.fitted = type.fitted ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "onempobs0"))[1]

    NOS <- extra$NOS
    M1 <- 2


    pobs.0 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                              .lpobs.0, earg = .epobs.0 ))
    lambda <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                              .llambda, earg = .elambda ))


    ans <- switch(type.fitted,
                  "mean"      = (1 - pobs.0) * lambda / (-expm1(-lambda)),
                  "pobs0"     =      pobs.0,  # P(Y=0)
                  "onempobs0" =  1 - pobs.0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))       
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpobs.0 = lpobs.0, .llambda = llambda,
           .epobs.0 = epobs.0, .elambda = elambda ))),
  last = eval(substitute(expression({
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE

    temp.names <- c(rep( .lpobs.0 , len = NOS),
                    rep( .llambda , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epobs.0
      misc$earg[[M1*ii  ]] <- .elambda
    }
  }), list( .lpobs.0 = lpobs.0, .llambda = llambda,
            .epobs.0 = epobs.0, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    NOS <- extra$NOS
    M1 <- 2

    pobs0  <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                              .lpobs.0, earg = .epobs.0))
    lambda <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                              .llambda, earg = .elambda ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dzapois(x = y, pobs0 = pobs0, lambda = lambda,
                                log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs.0 = lpobs.0, .llambda = llambda,
           .epobs.0 = epobs.0, .elambda = elambda ))),
  vfamily = c("zapoisson"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs0  <- eta2theta(eta[, c(TRUE, FALSE)], .lpobs.0 , earg = .epobs.0 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )
    rzapois(nsim * length(lambda), lambda = lambda, pobs0 = pobs0)
  }, list( .lpobs.0 = lpobs.0, .llambda = llambda,
           .epobs.0 = epobs.0, .elambda = elambda ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- extra$NOS
    y0 <- extra$y0
    skip <- extra$skip.these

    phimat <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                              .lpobs.0, earg = .epobs.0 ))
    lambda <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                              .llambda, earg = .elambda ))

    dl.dlambda <- y / lambda + 1 / expm1(-lambda)
    dl.dphimat <- -1 / (1 - phimat)  # For y > 0 obsns

    for (spp. in 1:NOS) {
      dl.dphimat[skip[, spp.], spp.] <- 1 / phimat[skip[, spp.], spp.]
      dl.dlambda[skip[, spp.], spp.] <- 0
    }
    dlambda.deta <- dtheta.deta(lambda, .llambda, earg = .elambda)
    mu.phi0 <- phimat

    temp3 <- if (.lpobs.0 == "logit") {
      c(w) * (y0 - mu.phi0)
    } else {
      c(w) * dtheta.deta(mu.phi0, link = .lpobs.0 , earg = .epobs.0 ) *
            dl.dphimat
    }

    ans <- cbind(temp3,
                 c(w) * dl.dlambda * dlambda.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M = M1)]
    ans
  }), list( .lpobs.0 = lpobs.0, .llambda = llambda,
            .epobs.0 = epobs.0, .elambda = elambda ))),
  weight = eval(substitute(expression({

    wz <- matrix(0.0, n, M1 * NOS)



    temp5 <- expm1(lambda)
    ned2l.dlambda2 <- (1 - phimat) * (temp5 + 1) *
                      (1 / lambda - 1 / temp5) / temp5
    wz[, NOS+(1:NOS)] <- w * ned2l.dlambda2 * dlambda.deta^2


    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if ( .lpobs.0 == "logit" && is.empty.list( .epobs.0 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi0, link = .lpobs.0, earg = .epobs.0)^2)
    }


  if (FALSE)
    for (ii in 1:NOS) {
      index200 <- abs(tmp200[, ii]) < .Machine$double.eps
      if (any(index200)) {
        tmp200[index200, ii] <- 10.0 * .Machine$double.eps^(3/4)
      }
    }


    wz[, 1:NOS] <-  tmp200

    wz <- wz[, interleave.VGAM(ncol(wz), M = M1)]



    wz
  }), list( .lpobs.0 = lpobs.0,
            .epobs.0 = epobs.0 ))))
}  # End of zapoisson





 zapoissonff <-
  function(llambda = "loge", lonempobs0 = "logit",
           type.fitted = c("mean", "pobs0", "onempobs0"),
           zero = -2) {



  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  lonempobs0 <- as.list(substitute(lonempobs0))
  eonempobs0 <- link2list(lonempobs0)
  lonempobs0 <- attr(eonempobs0, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "onempobs0"))[1]


  new("vglmff",
  blurb = c("Zero-altered Poisson ",
            "(Bernoulli and positive-Poisson conditional model)\n\n",
            "Links:    ",
            namesof("lambda",     llambda,    earg = elambda,
                    tag = FALSE), ", ",
            namesof("onempobs0",  lonempobs0, earg = eonempobs0,
                    tag = FALSE), "\n",
            "Mean:     onempobs0 * lambda / (1 - exp(-lambda))"),

  constraints = eval(substitute(expression({

    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 2
    if (any(y < 0))
      stop("the response must not have negative values")

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)

    extra$dimnamesy   <- dimnames(y)
    extra$type.fitted <- .type.fitted

    mynames1 <- if (ncoly == 1) "lambda"    else
                paste("lambda",    1:ncoly, sep = "")
    mynames2 <- if (ncoly == 1) "onempobs0" else
                paste("onempobs0", 1:ncoly, sep = "")

    predictors.names <-
        c(namesof(mynames1, .llambda,     earg = .elambda    , tag = FALSE),
          namesof(mynames2, .lonempobs0 , earg = .eonempobs0 , tag = FALSE))[
          interleave.VGAM(M1*NOS, M = M1)]

    if (!length(etastart)) {
      etastart <-
        cbind(matrix(1, n, NOS),  # 1 here is any old value
              theta2eta(1 - (0.5 + w * y0) / (1 + w),
                        .lonempobs0 , earg = .eonempobs0 ))
      for (spp. in 1:NOS) {
        sthese <- skip.these[, spp.]
        etastart[!sthese, 0 * NOS + spp.] <-
          theta2eta(y[!sthese, spp.] / (-expm1(-y[!sthese, spp.])),
                    .llambda, earg = .elambda )
      }
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = M1)]
    }
  }), list( .lonempobs0 = lonempobs0, .llambda = llambda,
            .eonempobs0 = eonempobs0, .elambda = elambda,
            .type.fitted = type.fitted ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "onempobs0"))[1]

    NOS <- extra$NOS
    M1 <- 2

    lambda    <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                                 .llambda    , earg = .elambda    ))
    onempobs0 <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                                 .lonempobs0 , earg = .eonempobs0 ))


    ans <- switch(type.fitted,
                  "mean"      =    (onempobs0) * lambda / (-expm1(-lambda)),
                  "pobs0"     = 1 - onempobs0,  # P(Y=0)
                  "onempobs0" =     onempobs0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempobs0 = lonempobs0, .llambda = llambda,
           .eonempobs0 = eonempobs0, .elambda = elambda ))),
  last = eval(substitute(expression({
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE

    temp.names <- c(rep( .llambda    , len = NOS),
                    rep( .lonempobs0 , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .elambda
      misc$earg[[M1*ii  ]] <- .eonempobs0
    }
  }), list( .lonempobs0 = lonempobs0, .llambda = llambda,
            .eonempobs0 = eonempobs0, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    NOS <- extra$NOS
    M1 <- 2

    lambda     <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                                  .llambda    , earg = .elambda    ))
    onempobs0  <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                                  .lonempobs0 , earg = .eonempobs0 ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzapois(x = y, lambda = lambda, pobs0 = 1 - onempobs0,
                       log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempobs0 = lonempobs0, .llambda = llambda,
           .eonempobs0 = eonempobs0, .elambda = elambda ))),
  vfamily = c("zapoissonff"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    lambda    <- eta2theta(eta[, c(TRUE, FALSE)], .llambda    ,
                           earg = .elambda    )
    onempobs0 <- eta2theta(eta[, c(FALSE, TRUE)], .lonempobs0 ,
                           earg = .eonempobs0 )
    rzapois(nsim * length(lambda), lambda = lambda, pobs0 = 1 - onempobs0)
  }, list( .lonempobs0 = lonempobs0, .llambda = llambda,
           .eonempobs0 = eonempobs0, .elambda = elambda ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- extra$NOS
    y0 <- extra$y0
    skip <- extra$skip.these

    lambda   <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                                .llambda, earg = .elambda ))
    omphimat <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                                .lonempobs0, earg = .eonempobs0 ))
    phimat <- 1 - omphimat


    dl.dlambda <- y / lambda + 1 / expm1(-lambda)
    dl.dPHImat <- +1 / (omphimat)  # For y > 0 obsns

    for (spp. in 1:NOS) {
      dl.dPHImat[skip[, spp.], spp.] <- -1 / phimat[skip[, spp.], spp.]
      dl.dlambda[skip[, spp.], spp.] <-  0
    }
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    mu.phi0 <- omphimat

    temp3 <- if ( FALSE && .lonempobs0 == "logit") {
    } else {
      c(w) * dtheta.deta(mu.phi0, link = .lonempobs0 , earg = .eonempobs0 ) *
            dl.dPHImat
    }

    ans <- cbind(c(w) * dl.dlambda * dlambda.deta,
                 temp3)
    ans <- ans[, interleave.VGAM(ncol(ans), M = M1)]
    ans
  }), list( .lonempobs0 = lonempobs0, .llambda = llambda,
            .eonempobs0 = eonempobs0, .elambda = elambda ))),
  weight = eval(substitute(expression({

    wz <- matrix(0.0, n, M1 * NOS)

    temp5 <- expm1(lambda)

    ned2l.dlambda2 <- (1 - phimat) * (temp5 + 1) *
                      (1 / lambda - 1 / temp5) / temp5


    wz[, 0 * NOS + (1:NOS)] <- c(w) * ned2l.dlambda2 * dlambda.deta^2


    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if ( .lonempobs0 == "logit" && is.empty.list( .eonempobs0 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi0, link = .lonempobs0, earg = .eonempobs0)^2)
    }


    wz[, 1 * NOS + (1:NOS)] <-  tmp200

    wz <- wz[, interleave.VGAM(ncol(wz), M = M1)]



    wz
  }), list( .lonempobs0 = lonempobs0,
            .eonempobs0 = eonempobs0 ))))
}  # End of zapoissonff







zanegbinomial.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 zanegbinomial <-
  function(lpobs0 = "logit", lmunb = "loge", lsize = "loge",
           type.fitted = c("mean", "pobs0"),
           ipobs0 = NULL,                    isize = NULL,
           zero = -3,  # Prior to 20130917 the default was: c(-1, -3),
           imethod = 1,
           nsimEIM = 250,
           ishrinkage = 0.95) {






  if (!is.Numeric(nsimEIM, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 30)
    warning("argument 'nsimEIM' should be greater than 30, say")


  if (length(ipobs0) && (!is.Numeric(ipobs0, positive = TRUE) ||
     max(ipobs0) >= 1))
    stop("If given, argument 'ipobs0' must contain values in (0,1) only")

  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("If given, argument 'isize' must contain positive values only")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")


  lpobs0 <- as.list(substitute(lpobs0))
  epobs0 <- link2list(lpobs0)
  lpobs0 <- attr(epobs0, "function.name")

  lmunb <- as.list(substitute(lmunb))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")


  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0"))[1]


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
    M1 <- 3
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 3

    if (any(y < 0))
      stop("the response must not have negative values")

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    M <- M1 * ncoly

    extra$dimnamesy   <- dimnames(y)
    extra$type.fitted <- .type.fitted

    mynames1 <- if (NOS == 1) "pobs0" else paste("pobs0", 1:NOS, sep = "")
    mynames2 <- if (NOS == 1) "munb"  else paste("munb",  1:NOS, sep = "")
    mynames3 <- if (NOS == 1) "size"  else paste("size",  1:NOS, sep = "")
    predictors.names <-
        c(namesof(mynames1, .lpobs0 , earg = .epobs0 , tag = FALSE),
          namesof(mynames2, .lmunb  , earg = .emunb  , tag = FALSE),
          namesof(mynames3, .lsize  , earg = .esize  , tag = FALSE))[
          interleave.VGAM(M1*NOS, M = M1)]


    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)


    if (!length(etastart)) {
      mu.init <- y
      for (iii in 1:ncol(y)) {
        index.posy <- (y[, iii] > 0)
        if ( .imethod == 1) {
          use.this <- weighted.mean(y[index.posy, iii],
                                    w[index.posy, iii])
          mu.init[ index.posy, iii] <- (1 - .ishrinkage ) * y[index.posy, iii] +
                                            .ishrinkage   * use.this
          mu.init[!index.posy, iii] <- use.this
        } else {
          use.this <-
          mu.init[, iii] <- (y[, iii] +
            weighted.mean(y[index.posy, iii],
                          w[index.posy, iii])) / 2
        }
        max.use.this <-  7 * use.this + 10
        vecTF <- (mu.init[, iii] > max.use.this)
        if (any(vecTF))
          mu.init[vecTF, iii] <- max.use.this
      }

      pnb0 <- matrix(if (length( .ipobs0 )) .ipobs0 else -1,
                     nrow = n, ncol = NOS, byrow = TRUE)
      for (spp. in 1:NOS) {
        if (any(pnb0[, spp.] < 0)) {
          index.y0 <- y[, spp.] < 0.5
          pnb0[, spp.] <- max(min(sum(index.y0) / n, 0.97), 0.03)
        }
      }


      if ( is.Numeric( .isize )) {
        kmat0 <- matrix( .isize , nrow = n, ncol = ncoly, byrow = TRUE)
      } else {
        posnegbinomial.Loglikfun <- function(kmat, y, x, w, extraargs) {
         munb <- extraargs
         sum(c(w) * dposnegbin(x = y, munb = munb, size = kmat,
                               log = TRUE))
        }
        k.grid <- 2^((-6):6)
        kmat0 <- matrix(0, nrow = n, ncol = NOS) 
        for (spp. in 1:NOS) {
          index.posy <- (y[, spp.] > 0)
          posy <- y[index.posy, spp.]
          kmat0[, spp.] <-
            grid.search(k.grid, objfun = posnegbinomial.Loglikfun,
                        y = posy, x = x[index.posy, ],
                        w = w[index.posy, spp.],
                        extraargs = mu.init[index.posy, spp.])
        }
      }

      etastart <- cbind(theta2eta(pnb0,    .lpobs0 , earg = .epobs0 ),
                        theta2eta(mu.init, .lmunb  , earg = .emunb  ),
                        theta2eta(kmat0,   .lsize  , earg = .esize  ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = M1)]
    }  # End of if (!length(etastart))


  }), list( .lpobs0 = lpobs0, .lmunb = lmunb, .lsize = lsize,
            .epobs0 = epobs0, .emunb = emunb, .esize = esize,
            .ipobs0 = ipobs0,                 .isize = isize,
            .imethod = imethod, .ishrinkage = ishrinkage,
            .type.fitted = type.fitted ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0"))[1]

    M1 <- 3
    NOS <- extra$NOS
    phi0 <- eta2theta(eta[, M1*(1:NOS)-2], .lpobs0 , earg = .epobs0 )
    munb <- eta2theta(eta[, M1*(1:NOS)-1], .lmunb  , earg = .emunb  )
    kmat <- eta2theta(eta[, M1*(1:NOS)  ], .lsize  , earg = .esize  )
    pnb0 <- (kmat / (kmat + munb))^kmat # p(0) from negative binomial


    ans <- switch(type.fitted,
                  "mean"      = (1 - phi0) * munb / (1 - pnb0),
                  "pobs0"     = phi0)  # P(Y=0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpobs0 = lpobs0, .lsize = lsize, .lmunb = lmunb,
           .epobs0 = epobs0, .emunb = emunb, .esize = esize ))),
  last = eval(substitute(expression({
    misc$link =
      c(rep( .lpobs0 , length = NOS),
        rep( .lmunb  , length = NOS),
        rep( .lsize  , length = NOS))[interleave.VGAM(M1*NOS,
                                                      M = M1)]
    temp.names <- c(mynames1,
                   mynames2,
                   mynames3)[interleave.VGAM(M1*NOS, M = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M1*NOS)
    names(misc$earg) <- temp.names
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-2]] <- .epobs0
      misc$earg[[M1*ii-1]] <- .emunb
      misc$earg[[M1*ii  ]] <- .esize
    }

    misc$nsimEIM <- .nsimEIM
    misc$imethod <- .imethod
    misc$ipobs0  <- .ipobs0
    misc$isize <- .isize
    misc$multipleResponses <- TRUE
  }), list( .lpobs0 = lpobs0, .lmunb = lmunb, .lsize = lsize,
            .epobs0 = epobs0, .emunb = emunb, .esize = esize,
            .ipobs0 = ipobs0, .isize = isize,
            .nsimEIM = nsimEIM,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    NOS <- extra$NOS
    M1 <- 3
    phi0 <- eta2theta(eta[, M1*(1:NOS)-2], .lpobs0 , earg = .epobs0 )
    munb <- eta2theta(eta[, M1*(1:NOS)-1], .lmunb  , earg = .emunb  )
    kmat <- eta2theta(eta[, M1*(1:NOS)  ], .lsize  , earg = .esize  )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzanegbin(x = y, pobs0 = phi0, munb = munb, size = kmat,
                         log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs0 = lpobs0, .lmunb = lmunb, .lsize = lsize,
           .epobs0 = epobs0, .emunb = emunb, .esize = esize ))),
  vfamily = c("zanegbinomial"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    phi0 <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lpobs0 , earg = .epobs0 )
    munb <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lmunb  , earg = .emunb  )
    kmat <- eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lsize  , earg = .esize  )
    rzanegbin(nsim * length(munb),
              pobs0 = phi0, munb = munb, size = kmat)
  }, list( .lpobs0 = lpobs0, .lmunb = lmunb, .lsize = lsize,
           .epobs0 = epobs0, .emunb = emunb, .esize = esize ))),





  deriv = eval(substitute(expression({
    M1 <- 3
    NOS <- extra$NOS
    y0 <- extra$y0

    phi0 <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                     .lpobs0 , earg = .epobs0 )
    munb <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                     .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                     .lsize , earg = .esize )
    skip <- extra$skip.these

    dphi0.deta <- dtheta.deta(phi0, .lpobs0 , earg = .epobs0 )
    dmunb.deta <- dtheta.deta(munb, .lmunb  , earg = .emunb  )
    dsize.deta <- dtheta.deta(kmat, .lsize  , earg = .esize  )


    tempk <- kmat / (kmat + munb)
    tempm <- munb / (kmat + munb)
    prob0  <- tempk^kmat
    oneminusf0  <- 1 - prob0
    df0.dmunb   <- -tempk * prob0
    df0.dkmat   <- prob0 * (tempm + log(tempk))


    dl.dphi0 <- -1 / (1 - phi0)
    dl.dmunb <- y / munb - (y + kmat) / (munb + kmat) +
                df0.dmunb / oneminusf0
    dl.dsize <- digamma(y + kmat) - digamma(kmat) -
                (y + kmat)/(munb + kmat) + 1 + log(tempk) +
                df0.dkmat / oneminusf0



    dl.dphi0[y == 0] <- 1 / phi0[y == 0]  # Do it in one line
    skip <- extra$skip.these
    for (spp. in 1:NOS) {
      dl.dsize[skip[, spp.], spp.] <-
      dl.dmunb[skip[, spp.], spp.] <- 0
    }

    dl.deta23 <- c(w) * cbind(dl.dmunb * dmunb.deta,
                              dl.dsize * dsize.deta)


    muphi0 <- phi0
    dl.deta1 <- if ( .lpobs0 == "logit") {
      c(w) * (y0 - muphi0)
    } else {
      c(w) * dphi0.deta * (y0 / muphi0 - 1) / (1 - muphi0)
    }
    ans <- cbind(dl.deta1, dl.deta23)
    ans <- ans[, interleave.VGAM(ncol(ans), M = M1)]
    ans
  }), list( .lpobs0 = lpobs0 , .lmunb = lmunb , .lsize = lsize ,
            .epobs0 = epobs0 , .emunb = emunb , .esize = esize  ))),

  weight = eval(substitute(expression({

    six <- dimm(M1)
    wz <- run.varcov <- matrix(0.0, n, six*NOS-1)
    M1m1 <- M1 - 1






    ind2 <- iam(NA, NA, M = M1 - 1, both = TRUE, diag = TRUE)


    for (ii in 1:( .nsimEIM )) {
      ysim <- rzanegbin(n = n*NOS, pobs0 = phi0,
                        size = kmat, mu = munb)
      dim(ysim) <- c(n, NOS)




      dl.dphi0 <- -1 / (1 - phi0)
      dl.dmunb <- ysim / munb - (ysim + kmat) / (munb + kmat) +
                  df0.dmunb / oneminusf0
      dl.dsize <- digamma(ysim + kmat) - digamma(kmat) -
                  (ysim + kmat)/(munb + kmat) + 1 + log(tempk) +
                  df0.dkmat / oneminusf0




      dl.dphi0[ysim == 0] <- 1 / phi0[ysim == 0]  # Do it in one line
      ysim0 <- ifelse(ysim == 0, 1, 0)
      skip.sim <- matrix(as.logical(ysim0), n, NOS)
      for (spp. in 1:NOS) {
        dl.dsize[skip.sim[, spp.], spp.] <-
        dl.dmunb[skip.sim[, spp.], spp.] <- 0
      }


      for (kk in 1:NOS) {
        temp2 <- cbind(dl.dmunb[, kk] * dmunb.deta[, kk],
                       dl.dsize[, kk] * dsize.deta[, kk])
        small.varcov <- temp2[, ind2$row.index] *
                       temp2[, ind2$col.index]




        run.varcov[, ((kk-1)*M1+2):(kk*M1)] <-
        run.varcov[, ((kk-1)*M1+2):(kk*M1)] +
          c(small.varcov[, 1:M1m1])
        run.varcov[, M + (kk-1)*M1 + 2] <-
        run.varcov[, M + (kk-1)*M1 + 2] +
          c(small.varcov[, M1m1 + 1])
      }  # kk; end of NOS
    }  # ii; end of nsimEIM


    run.varcov <- cbind(run.varcov / .nsimEIM )
    run.varcov <- if (intercept.only)
      matrix(colMeans(run.varcov),
             n, ncol(run.varcov), byrow = TRUE) else run.varcov




    wzind1 <- sort(c(    M1*(1:NOS) - 1,
                         M1*(1:NOS) - 0,
                     M + M1*(1:NOS) - 1))
    wz[, wzind1] <- c(w) * run.varcov[, wzind1]




    tmp100 <- muphi0 * (1 - muphi0)
    tmp200 <- if ( .lpobs0 == "logit") {
      cbind(c(w) * tmp100)
    } else {
      c(w) * cbind(dphi0.deta^2 / tmp100)
    }
    for (ii in 1:NOS) {
      index200 <- abs(tmp200[, ii]) < .Machine$double.eps
      if (any(index200)) {
        tmp200[index200, ii] <- .Machine$double.eps  # Diagonal 0's are bad 
      }
    }
    wz[, M1*(1:NOS)-2] <- tmp200



    wz
  }), list( .lpobs0 = lpobs0,
            .epobs0 = epobs0,
            .nsimEIM = nsimEIM ))))
}  # End of zanegbinomial()




zanegbinomialff.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 zanegbinomialff <-
  function(lmunb = "loge", lsize = "loge", lonempobs0 = "logit",
           type.fitted = c("mean", "pobs0", "onempobs0"),
           isize = NULL, ionempobs0 = NULL,
           zero = c(-2, -3),
           imethod = 1,
           nsimEIM = 250,
           ishrinkage = 0.95) {



  if (!is.Numeric(nsimEIM, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 30)
    warning("argument 'nsimEIM' should be greater than 30, say")


  if (length(ionempobs0) && (!is.Numeric(ionempobs0, positive = TRUE) ||
     max(ionempobs0) >= 1))
    stop("If given, argument 'ionempobs0' must contain values in (0,1) only")

  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("If given, argument 'isize' must contain positive values only")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")

  lmunb <- as.list(substitute(lmunb))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")

  lonempobs0 <- as.list(substitute(lonempobs0))
  eonempobs0 <- link2list(lonempobs0)
  lonempobs0 <- attr(eonempobs0, "function.name")


  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "onempobs0"))[1]


  new("vglmff",
  blurb = c("Zero-altered negative binomial (Bernoulli and\n",
            "positive-negative binomial conditional model)\n\n",
            "Links:    ",
            namesof("munb",  lmunb,  earg = emunb,  tag = FALSE), ", ",
            namesof("size",  lsize,  earg = esize,  tag = FALSE), ", ",
            namesof("onempobs0", lonempobs0, earg = eonempobs0,
                    tag = FALSE), "\n",
            "Mean:     onempobs0 * munb / (1 - (size / (size + ",
                                                 "munb))^size)"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    M1 <- 3
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 3

    if (any(y < 0))
      stop("the response must not have negative values")

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    M <- M1 * ncoly

    extra$dimnamesy   <- dimnames(y)
    extra$type.fitted <- .type.fitted

    mynames1 <- if (NOS == 1) "munb"  else paste("munb",  1:NOS, sep = "")
    mynames2 <- if (NOS == 1) "size"  else paste("size",  1:NOS, sep = "")
    mynames3 <- if (NOS == 1) "onempobs0" else paste("onempobs0", 1:NOS,
                                                     sep = "")
    predictors.names <-
        c(namesof(mynames1, .lmunb  , earg = .emunb  , tag = FALSE),
          namesof(mynames2, .lsize  , earg = .esize  , tag = FALSE),
          namesof(mynames3, .lonempobs0 , earg = .eonempobs0 ,
                  tag = FALSE))[
          interleave.VGAM(M1*NOS, M = M1)]


    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)


    if (!length(etastart)) {
      mu.init <- y
      for (iii in 1:ncol(y)) {
        index.posy <- (y[, iii] > 0)
        if ( .imethod == 1) {
          use.this <- weighted.mean(y[index.posy, iii],
                                    w[index.posy, iii])
          mu.init[ index.posy, iii] <- (1 - .ishrinkage ) * y[index.posy, iii] +
                                            .ishrinkage   * use.this
          mu.init[!index.posy, iii] <- use.this
        } else {
          use.this <-
          mu.init[, iii] <- (y[, iii] +
            weighted.mean(y[index.posy, iii],
                          w[index.posy, iii])) / 2
        }
        max.use.this <-  7 * use.this + 10
        vecTF <- (mu.init[, iii] > max.use.this)
        if (any(vecTF))
          mu.init[vecTF, iii] <- max.use.this
      }

      pnb0 <- matrix(if (length( .ionempobs0 )) 1 - .ionempobs0 else -1,
                     nrow = n, ncol = NOS, byrow = TRUE)
      for (spp. in 1:NOS) {
        if (any(pnb0[, spp.] < 0)) {
          index.y0 <- y[, spp.] < 0.5
          pnb0[, spp.] <- max(min(sum(index.y0) / n, 0.97), 0.03)
        }
      }


      if ( is.Numeric( .isize )) {
        kmat0 <- matrix( .isize , nrow = n, ncol = ncoly, byrow = TRUE)
      } else {
        posnegbinomial.Loglikfun <- function(kmat, y, x, w, extraargs) {
         munb <- extraargs
         sum(c(w) * dposnegbin(x = y, munb = munb, size = kmat,
                               log = TRUE))
        }
        k.grid <- 2^((-6):6)
        kmat0 <- matrix(0, nrow = n, ncol = NOS) 
        for (spp. in 1:NOS) {
          index.posy <- (y[, spp.] > 0)
          posy <- y[index.posy, spp.]
          kmat0[, spp.] <-
            grid.search(k.grid, objfun = posnegbinomial.Loglikfun,
                        y = posy, x = x[index.posy, ],
                        w = w[index.posy, spp.],
                        extraargs = mu.init[index.posy, spp.])
        }
      }

      etastart <-
        cbind(theta2eta(mu.init , .lmunb      , earg = .emunb      ),
              theta2eta(kmat0   , .lsize      , earg = .esize      ),
              theta2eta(1 - pnb0, .lonempobs0 , earg = .eonempobs0 ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = M1)]
    }  # End of if (!length(etastart))


  }), list( .lonempobs0 = lonempobs0, .lmunb = lmunb, .lsize = lsize,
            .eonempobs0 = eonempobs0, .emunb = emunb, .esize = esize,
            .ionempobs0 = ionempobs0,                 .isize = isize,
            .imethod = imethod, .ishrinkage = ishrinkage,
            .type.fitted = type.fitted ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "onempobs0"))[1]

    M1 <- 3
    NOS <- extra$NOS
    munb <- eta2theta(eta[, M1*(1:NOS)-2], .lmunb  , earg = .emunb  )
    kmat <- eta2theta(eta[, M1*(1:NOS)-1], .lsize  , earg = .esize  )
    onempobs0 <- eta2theta(eta[, M1*(1:NOS)  ], .lonempobs0 ,
                           earg = .eonempobs0 )
    pnb0 <- (kmat / (kmat + munb))^kmat  # p(0) from negative binomial


    ans <- switch(type.fitted,
                  "mean"      =    (onempobs0) * munb / (1 - pnb0),
                  "pobs0"     = 1 - onempobs0,  # P(Y=0)
                  "onempobs0" =     onempobs0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))        
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempobs0 = lonempobs0, .lsize = lsize, .lmunb = lmunb,
           .eonempobs0 = eonempobs0, .emunb = emunb, .esize = esize ))),
  last = eval(substitute(expression({
    misc$link <-
      c(rep( .lmunb      , length = NOS),
        rep( .lsize      , length = NOS),
        rep( .lonempobs0 , length = NOS))[
        interleave.VGAM(M1*NOS, M = M1)]
    temp.names <- c(mynames1,
                    mynames2,
                    mynames3)[interleave.VGAM(M1*NOS, M = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M1*NOS)
    names(misc$earg) <- temp.names
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-2]] <- .emunb
      misc$earg[[M1*ii-1]] <- .esize
      misc$earg[[M1*ii  ]] <- .eonempobs0
    }

    misc$nsimEIM <- .nsimEIM
    misc$imethod <- .imethod
    misc$ionempobs0  <- .ionempobs0
    misc$isize <- .isize
    misc$multipleResponses <- TRUE
  }), list( .lonempobs0 = lonempobs0, .lmunb = lmunb, .lsize = lsize,
            .eonempobs0 = eonempobs0, .emunb = emunb, .esize = esize,
            .ionempobs0 = ionempobs0, .isize = isize,
            .nsimEIM = nsimEIM,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    NOS <- extra$NOS
    M1 <- 3
    munb <- eta2theta(eta[, M1*(1:NOS)-2], .lmunb  , earg = .emunb  )
    kmat <- eta2theta(eta[, M1*(1:NOS)-1], .lsize  , earg = .esize  )
    onempobs0 <- eta2theta(eta[, M1*(1:NOS)  ], .lonempobs0 ,
                           earg = .eonempobs0 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzanegbin(x = y, pobs0 = 1 - onempobs0,
                         munb = munb, size = kmat,
                         log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempobs0 = lonempobs0, .lmunb = lmunb, .lsize = lsize,
           .eonempobs0 = eonempobs0, .emunb = emunb, .esize = esize ))),
  vfamily = c("zanegbinomialff"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    munb      <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lmunb  , earg = .emunb  )
    kmat      <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lsize  , earg = .esize  )
    onempobs0 <- eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lonempobs0 ,
                           earg = .eonempobs0 )

    rzanegbin(nsim * length(munb),
              pobs0 = 1 - onempobs0, munb = munb, size = kmat)
  }, list( .lonempobs0 = lonempobs0, .lmunb = lmunb, .lsize = lsize,
           .eonempobs0 = eonempobs0, .emunb = emunb, .esize = esize ))),




  deriv = eval(substitute(expression({
    M1 <- 3
    NOS <- extra$NOS
    y0 <- extra$y0

    munb      <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                           .lmunb      , earg = .emunb )
    kmat      <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                           .lsize      , earg = .esize )
    onempobs0 <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                           .lonempobs0 , earg = .eonempobs0 )
    skip <- extra$skip.these
    phi0 <- 1 - onempobs0

    dmunb.deta      <- dtheta.deta(munb, .lmunb  , earg = .emunb  )
    dsize.deta      <- dtheta.deta(kmat, .lsize  , earg = .esize  )
    donempobs0.deta <- dtheta.deta(onempobs0, .lonempobs0 ,
                                   earg = .eonempobs0 )


    tempk <- kmat / (kmat + munb)
    tempm <- munb / (kmat + munb)
    prob0  <- tempk^kmat
    oneminusf0  <- 1 - prob0
    df0.dmunb   <- -tempk * prob0
    df0.dkmat   <- prob0 * (tempm + log(tempk))


    dl.dmunb <- y / munb - (y + kmat) / (munb + kmat) +
                df0.dmunb / oneminusf0
    dl.dsize <- digamma(y + kmat) - digamma(kmat) -
                (y + kmat)/(munb + kmat) + 1 + log(tempk) +
                df0.dkmat / oneminusf0
    dl.donempobs0 <- +1 / (onempobs0)



    dl.donempobs0[y == 0] <-
      -1 / (1 - onempobs0[y == 0])  # Do it in 1 line
    skip <- extra$skip.these
    for (spp. in 1:NOS) {
      dl.dsize[skip[, spp.], spp.] <-
      dl.dmunb[skip[, spp.], spp.] <- 0
    }

    dl.deta12 <- c(w) * cbind(dl.dmunb * dmunb.deta,
                              dl.dsize * dsize.deta)


    muphi0 <- onempobs0  # Originally: phi0
    dl.deta3 <- if (FALSE &&
                    .lonempobs0 == "logit") {
    } else {

      c(w) * donempobs0.deta * dl.donempobs0
    }
    ans <- cbind(dl.deta12, dl.deta3)
    ans <- ans[, interleave.VGAM(ncol(ans), M = M1)]
    ans
  }), list( .lonempobs0 = lonempobs0 , .lmunb = lmunb , .lsize = lsize ,
            .eonempobs0 = eonempobs0 , .emunb = emunb , .esize = esize  ))),

  weight = eval(substitute(expression({

    six <- dimm(M1)
    wz <- run.varcov <- matrix(0.0, n, six*NOS-1)
    M1m1 <- M1 - 1





    ind2 <- iam(NA, NA, M = M1 - 1, both = TRUE, diag = TRUE)


    for (ii in 1:( .nsimEIM )) {
      ysim <- rzanegbin(n = n*NOS, pobs0 = phi0,
                        size = kmat, mu = munb)
      dim(ysim) <- c(n, NOS)


      dl.dmunb <- ysim / munb - (ysim + kmat) / (munb + kmat) +
                  df0.dmunb / oneminusf0
      dl.dsize <- digamma(ysim + kmat) - digamma(kmat) -
                  (ysim + kmat)/(munb + kmat) + 1 + log(tempk) +
                  df0.dkmat / oneminusf0
      dl.donempobs0 <- +1 / (onempobs0)



      dl.donempobs0[ysim == 0] <-
        -1 / (1 - onempobs0[ysim == 0])  # Do it in 1 line
      ysim0 <- ifelse(ysim == 0, 1, 0)
      skip.sim <- matrix(as.logical(ysim0), n, NOS)
      for (spp. in 1:NOS) {
        dl.dsize[skip.sim[, spp.], spp.] <-
        dl.dmunb[skip.sim[, spp.], spp.] <- 0
      }


      for (kk in 1:NOS) {
        temp2 <- cbind(dl.dmunb[, kk] * dmunb.deta[, kk],
                       dl.dsize[, kk] * dsize.deta[, kk])
        small.varcov <- temp2[, ind2$row.index] *
                        temp2[, ind2$col.index]


        run.varcov[, ((kk-1)*M1+2-1):(kk*M1-1)] <-
        run.varcov[, ((kk-1)*M1+2-1):(kk*M1-1)] +
          c(small.varcov[, 1:M1m1])
        run.varcov[, M + (kk-1)*M1 + 2-1] <-
        run.varcov[, M + (kk-1)*M1 + 2-1] +
          c(small.varcov[, M1m1 + 1])
      }  # kk; end of NOS
    }  # ii; end of nsimEIM


    run.varcov <- cbind(run.varcov / .nsimEIM )
    run.varcov <- if (intercept.only)
      matrix(colMeans(run.varcov),
             n, ncol(run.varcov), byrow = TRUE) else run.varcov



    wzind1 <- sort(c(    M1*(1:NOS) - 1 - 1,
                         M1*(1:NOS) - 0 - 1,
                     M + M1*(1:NOS) - 1 - 1))
    wz[, wzind1] <- c(w) * run.varcov[, wzind1]


    tmp100 <- muphi0 * (1 - muphi0)
    tmp200 <- if (FALSE &&
                  .lpobs0 == "logit") {
    } else {
      c(w) * cbind(donempobs0.deta^2 / tmp100)
    }
    for (ii in 1:NOS) {
      index200 <- abs(tmp200[, ii]) < .Machine$double.eps
      if (any(index200)) {
        tmp200[index200, ii] <- .Machine$double.eps  # Diagonal 0's are bad 
      }
    }
    wz[, M1*(1:NOS)  ] <- tmp200



    wz
  }), list( .lonempobs0 = lonempobs0,
            .eonempobs0 = eonempobs0,
            .nsimEIM = nsimEIM ))))
}  # End of zanegbinomialff()











 zipoisson <-
  function(lpstr0 = "logit", llambda = "loge",
           type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
           ipstr0 = NULL,    ilambda = NULL,
           imethod = 1,
           ishrinkage = 0.8, zero = NULL) {
  ipstr00 <- ipstr0


  lpstr0 <- as.list(substitute(lpstr0))
  epstr00 <- link2list(lpstr0)
  lpstr00 <- attr(epstr00, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")



  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "pstr0", "onempstr0"))[1]


  if (length(ipstr00))
    if (!is.Numeric(ipstr00, positive = TRUE) ||
        any(ipstr00 >= 1))
      stop("argument 'ipstr0' values must be inside the interval (0,1)")
  if (length(ilambda))
    if (!is.Numeric(ilambda, positive = TRUE))
      stop("argument 'ilambda' values must be positive")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")


  new("vglmff",
  blurb = c("Zero-inflated Poisson\n\n",
            "Links:    ",
            namesof("pstr0",  lpstr00, earg = epstr00 ), ", ",
            namesof("lambda", llambda, earg = elambda ), "\n",
            "Mean:     (1 - pstr0) * lambda"),

  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    ncoly <- ncol(y)
    M1 <- 2
    extra$ncoly <- ncoly
    extra$M1 <- M1
    extra$dimnamesy <- dimnames(y)
    M <- M1 * ncoly
    extra$type.fitted      <- .type.fitted


    if (any(round(y) != y))
      stop("integer-valued responses only allowed for ",
           "the 'zipoisson' family")

    mynames1 <- paste("pstr0",   if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("lambda",  if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lpstr00 , earg = .epstr00 , tag = FALSE),
          namesof(mynames2, .llambda , earg = .elambda , tag = FALSE))[
          interleave.VGAM(M, M = M1)]



    if (!length(etastart)) {

      matL <- matrix(if (length( .ilambda )) .ilambda else 0,
                     n, ncoly, byrow = TRUE)
      matP <- matrix(if (length( .ipstr00 )) .ipstr00 else 0,
                     n, ncoly, byrow = TRUE)


      for (spp. in 1:ncoly) {
        yvec <- y[, spp.]

        Phi.init <- 1 - 0.85 * sum(w[yvec > 0]) / sum(w)
        Phi.init[Phi.init <= 0.02] <- 0.02 # Last resort
        Phi.init[Phi.init >= 0.98] <- 0.98 # Last resort

        if ( length(mustart)) {
          mustart <- matrix(mustart, n, ncoly)  # Make sure right size
          Lambda.init <- mustart / (1 - Phi.init)
        } else if ( .imethod == 2) {
          mymean <- weighted.mean(yvec[yvec > 0],
                                     w[yvec > 0]) + 1/16
          Lambda.init <- (1 - .ishrinkage ) * (yvec + 1/8) + .ishrinkage * mymean
        } else {
          use.this <- median(yvec[yvec > 0]) + 1 / 16
          Lambda.init <- (1 - .ishrinkage ) * (yvec + 1/8) + .ishrinkage * use.this
        }

        zipois.Loglikfun <- function(phival, y, x, w, extraargs) {
          sum(c(w) * dzipois(x = y, pstr0 = phival,
                          lambda = extraargs$lambda,
                          log = TRUE))
        }
        phi.grid <- seq(0.02, 0.98, len = 21)
        Phimat.init <- grid.search(phi.grid, objfun = zipois.Loglikfun,
                                   y = y, x = x, w = w,
                                   extraargs = list(lambda = Lambda.init))

        if (length(mustart)) {
          Lambda.init <- Lambda.init / (1 - Phimat.init)
        }

        if (!length( .ipstr00 ))
          matP[, spp.] <- Phimat.init
        if (!length( .ilambda ))
          matL[, spp.] <- Lambda.init
      }  # spp.

      etastart <- cbind(theta2eta(matP, .lpstr00, earg = .epstr00 ),
                        theta2eta(matL, .llambda, earg = .elambda ))[,
                        interleave.VGAM(M, M = M1)]
      mustart <- NULL  # Since etastart has been computed.
    }  # End of !length(etastart)
  }), list( .lpstr00 = lpstr00, .llambda = llambda,
            .epstr00 = epstr00, .elambda = elambda,
            .ipstr00 = ipstr00, .ilambda = ilambda,
            .imethod = imethod,
            .type.fitted = type.fitted,
            .ishrinkage = ishrinkage ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "pstr0", "onempstr0"))[1]

    phimat <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr00 , earg = .epstr00 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )

    
    ans <- switch(type.fitted,
                  "mean"      = (1 - phimat) * lambda,
                  "pobs0"     = phimat + (1-phimat)*exp(-lambda),  # P(Y=0)
                  "pstr0"     =     phimat,
                  "onempstr0" = 1 - phimat)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans)) 
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpstr00 = lpstr00, .llambda = llambda,
           .epstr00 = epstr00, .elambda = elambda,
           .type.fitted = type.fitted
         ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <-
      c(rep( .lpstr00 , length = ncoly),
        rep( .llambda , length = ncoly))[interleave.VGAM(M, M = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .epstr00
      misc$earg[[M1*ii  ]] <- .elambda
    }

    misc$M1 <- M1
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE

      misc$pobs0 <- phimat + (1 - phimat) * exp(-lambda)  # P(Y=0)
      if (length(dimnames(y)[[2]]) > 0)
        dimnames(misc$pobs0) <- dimnames(y)

      misc$pstr0 <- phimat
      if (length(dimnames(y)[[2]]) > 0)
        dimnames(misc$pstr0) <- dimnames(y)
  }), list( .lpstr00 = lpstr00, .llambda = llambda,
            .epstr00 = epstr00, .elambda = elambda,
            .imethod = imethod ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    phimat <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr00 , earg = .epstr00 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dzipois(x = y, pstr0 = phimat, lambda = lambda,
                                log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr00 = lpstr00, .llambda = llambda,
           .epstr00 = epstr00, .elambda = elambda ))),
  vfamily = c("zipoisson"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    phimat <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr00 , earg = .epstr00 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )
    rzipois(nsim * length(lambda), lambda = lambda, pstr0 = phimat)
  }, list( .lpstr00 = lpstr00, .llambda = llambda,
           .epstr00 = epstr00, .elambda = elambda ))),




  deriv = eval(substitute(expression({
    M1 <- 2
    phimat <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr00 ,
                        earg = .epstr00 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .llambda ,
                        earg = .elambda )

    prob0 <- exp(-lambda)
    pobs0 <- phimat + (1 - phimat) * prob0
    index0 <- as.matrix(y == 0)

    dl.dphimat <- -expm1(-lambda) / pobs0
    dl.dphimat[!index0] <- -1 / (1 - phimat[!index0])

    dl.dlambda <- -(1 - phimat) * exp(-lambda) / pobs0
    dl.dlambda[!index0] <- (y[!index0] - lambda[!index0]) / lambda[!index0]

    dphimat.deta <- dtheta.deta(phimat, .lpstr00 , earg = .epstr00 )
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )

    ans <- c(w) * cbind(dl.dphimat * dphimat.deta,
                        dl.dlambda * dlambda.deta)
    ans <- ans[, interleave.VGAM(M, M = M1)]


    if ( .llambda == "loge" && is.empty.list( .elambda ) &&
       any(lambda[!index0] < .Machine$double.eps)) {
      for (spp. in 1:(M / M1)) {
        ans[!index0[, spp.], M1 * spp.] <-
          w[!index0[, spp.]] *
         (y[!index0[, spp.], spp.] - lambda[!index0[, spp.], spp.])
      }
    }

    ans
  }), list( .lpstr00 = lpstr00, .llambda = llambda,
            .epstr00 = epstr00, .elambda = elambda ))),
  weight = eval(substitute(expression({

    ned2l.dphimat2 <- -expm1(-lambda) / ((1 - phimat) * pobs0)
    ned2l.dphimatlambda <- -exp(-lambda) / pobs0
    ned2l.dlambda2 <- (1 - phimat) / lambda -
                      phimat * (1 - phimat) * exp(-lambda) / pobs0




    wz <- array(c(c(w) * ned2l.dphimat2 * dphimat.deta^2,
                  c(w) * ned2l.dlambda2 * dlambda.deta^2,
                  c(w) * ned2l.dphimatlambda * dphimat.deta * dlambda.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)





    wz
  }), list( .llambda = llambda, .elambda = elambda ))))
}  # zipoisson









 zibinomial <-
  function(lpstr0 = "logit", lprob = "logit",
           type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
           ipstr0 = NULL,
           zero = NULL,  # 20130917; was originally zero = 1,
           multiple.responses = FALSE, imethod = 1) {
  if (as.logical(multiple.responses))
    stop("argument 'multiple.responses' must be FALSE")

  lpstr0 <- as.list(substitute(lpstr0))
  epstr0 <- link2list(lpstr0)
  lpstr0 <- attr(epstr0, "function.name")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "pstr0", "onempstr0"))[1]


  if (is.Numeric(ipstr0))
    if (!is.Numeric(ipstr0, positive = TRUE) || any(ipstr0 >= 1))
      stop("'ipstr0' values must be inside the interval (0,1)")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")



  new("vglmff",
  blurb = c("Zero-inflated binomial\n\n",
            "Links:    ",
            namesof("pstr0", lpstr0, earg = epstr0), ", ",
            namesof("prob" , lprob , earg = eprob ), "\n",
            "Mean:     (1 - pstr0) * prob"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         type.fitted  = .type.fitted ,
         expected = TRUE,
         multiple.responses  = FALSE,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),
      
  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w



    if (NCOL(y) == 1) {
      if (is.factor(y))
        y <- y != levels(y)[1]
      nn <- rep(1, n)
      if (!all(y >= 0 & y <= 1))
        stop("response values must be in [0, 1]")
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + w * y) / (1.0 + w)


      no.successes <- y
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(no.successes - round(no.successes)) > 1.0e-8))
        stop("Number of successes must be integer-valued")

    } else if (NCOL(y) == 2) {
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(y - round(y)) > 1.0e-8))
        stop("Count data must be integer-valued")
      y <- round(y)
      nvec <- y[, 1] + y[, 2]
      y <- ifelse(nvec > 0, y[, 1] / nvec, 0)
      w <- w * nvec
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + nvec * y) / (1 + nvec)
    } else {
      stop("for the binomialff family, response 'y' must be a ",
           "vector of 0 and 1's\n",
           "or a factor ",
           "(first level = fail, other levels = success),\n",
           "or a 2-column matrix where col 1 is the no. of ",
           "successes and col 2 is the no. of failures")
    }


    if ( .imethod == 1)
      mustart <- (mustart + y) / 2


    extra$type.fitted <- .type.fitted
    extra$dimnamesy   <- dimnames(y)




    predictors.names <-
        c(namesof("pstr0", .lpstr0 , earg = .epstr0 , tag = FALSE),
          namesof("prob" , .lprob  , earg = .eprob  , tag = FALSE))


    extra$w <- w  # Needed for @linkinv
    phi.init <- if (length( .ipstr0 )) .ipstr0 else {
        prob0.est <- sum(w[y == 0]) / sum(w)
        if ( .imethod == 1) {
          (prob0.est - (1 - mustart)^w) / (1 - (1 - mustart)^w)
        } else {
          prob0.est
        }
    }

    phi.init[phi.init <= -0.10] <- 0.10  # Lots of sample variation
    phi.init[phi.init <=  0.05] <- 0.15  # Last resort
    phi.init[phi.init >=  0.80] <- 0.80  # Last resort

    if ( length(mustart) && !length(etastart))
      mustart <- cbind(rep(phi.init, len = n),
                       mustart)  # 1st coln not a real mu
  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob,
            .ipstr0 = ipstr0,
            .type.fitted = type.fitted,          
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    pstr0 <- eta2theta(eta[, 1], .lpstr0 , earg = .epstr0 )
    mubin <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )


    orig.w <- if (length(tmp3 <- extra$orig.w)) tmp3 else
              rep(1, len = nrow(eta))
    priorw <- extra$w
    nvec <- priorw / orig.w


    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "pstr0", "onempstr0"))[1]

    ans <- switch(type.fitted,
                  "mean"      = (1 - pstr0) * mubin,
                  "pobs0"     = pstr0 + (1-pstr0)*(1-mubin)^nvec,  # P(Y=0)
                  "pstr0"     =     pstr0,
                  "onempstr0" = 1 - pstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpstr0 = lpstr0, .lprob = lprob,
           .epstr0 = epstr0, .eprob = eprob,
           .type.fitted = type.fitted ))),
  last = eval(substitute(expression({
    misc$link <-    c("pstr0" = .lpstr0 , "prob" = .lprob )

    misc$earg <- list("pstr0" = .epstr0 , "prob" = .eprob )

    misc$imethod <- .imethod


  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob,
            .imethod = imethod ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    cbind(theta2eta(mu[, 1], .lpstr0 , earg = .epstr0 ),
          theta2eta(mu[, 2], .lprob  , earg = .eprob  ))
  }, list( .lpstr0 = lpstr0, .lprob = lprob,
           .epstr0 = epstr0, .eprob = eprob ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pstr0 <- eta2theta(eta[, 1], .lpstr0 , earg = .epstr0 )
    mubin <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        dzibinom(x = round(w * y), size = w, prob = mubin,
                 log = TRUE, pstr0 = pstr0)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr0 = lpstr0, .lprob = lprob,
           .epstr0 = epstr0, .eprob = eprob ))),
  vfamily = c("zibinomial"),
  deriv = eval(substitute(expression({
    phi   <- eta2theta(eta[, 1], .lpstr0 , earg = .epstr0 )
    mubin <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )

    prob0 <- (1 - mubin)^w  # Actually q^w
    pobs0 <- phi + (1 - phi) * prob0
    index <- (y == 0)
    dl.dphi <- (1 - prob0) / pobs0
    dl.dphi[!index] <- -1 / (1 - phi[!index])

    dl.dmubin <- -w * (1 - phi) * (1 - mubin)^(w - 1) / pobs0
    dl.dmubin[!index] <- w[!index] *
        (    y[!index]  /      mubin[!index]   -
        (1 - y[!index]) / (1 - mubin[!index]))

    dphi.deta   <- dtheta.deta(phi,   .lpstr0 , earg = .epstr0 )
    dmubin.deta <- dtheta.deta(mubin, .lprob  , earg = .eprob  )

    ans <- cbind(dl.dphi   * dphi.deta,
                 dl.dmubin * dmubin.deta)

      if ( .lprob == "logit") {
        ans[!index, 2] <- w[!index] * (y[!index] - mubin[!index])
      }

      ans
  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob ))),
  weight = eval(substitute(expression({
    wz <- matrix(as.numeric(NA), nrow = n, ncol = dimm(M))



    ned2l.dphi2 <- (1 - prob0) / ((1 - phi) * pobs0)


    ned2l.dphimubin <- -w * ((1 - mubin)^(w - 1)) / pobs0







    ned2l.dmubin2 <- (w * (1 - phi) / (mubin * (1 - mubin)^2)) *
                     (1 - mubin - w * mubin *
                     (1 - mubin)^w * phi / pobs0)





    wz[,iam(1, 1, M)] <- ned2l.dphi2     * dphi.deta^2
    wz[,iam(2, 2, M)] <- ned2l.dmubin2   * dmubin.deta^2
    wz[,iam(1, 2, M)] <- ned2l.dphimubin * dphi.deta * dmubin.deta
    if (TRUE) {
      ind6 <- (wz[, iam(2, 2, M)] < .Machine$double.eps)
      if (any(ind6))
        wz[ind6, iam(2, 2, M)] <- .Machine$double.eps
    }
    wz
  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob ))))
}






 zibinomialff <-
  function(lprob = "logit", lonempstr0 = "logit",
           type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
           ionempstr0 = NULL,
           zero = 2,
           multiple.responses = FALSE, imethod = 1) {






  if (as.logical(multiple.responses))
    stop("argument 'multiple.responses' must be FALSE")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lonempstr0 <- as.list(substitute(lonempstr0))
  eonempstr0 <- link2list(lonempstr0)
  lonempstr0 <- attr(eonempstr0, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "pstr0", "onempstr0"))[1]


  if (is.Numeric(ionempstr0))
    if (!is.Numeric(ionempstr0, positive = TRUE) || any(ionempstr0 >= 1))
      stop("'ionempstr0' values must be inside the interval (0,1)")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")



  new("vglmff",
  blurb = c("Zero-inflated binomial\n\n",
            "Links:    ",
            namesof("prob" ,     lprob     , earg = eprob     ), ", ",
            namesof("onempstr0", lonempstr0, earg = eonempstr0), "\n",
            "Mean:     onempstr0 * prob"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),
      
  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w



    if (NCOL(y) == 1) {
      if (is.factor(y))
        y <- y != levels(y)[1]
      nn <- rep(1, n)
      if (!all(y >= 0 & y <= 1))
        stop("response values must be in [0, 1]")
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + w * y) / (1.0 + w)


      no.successes <- y
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(no.successes - round(no.successes)) > 1.0e-8))
        stop("Number of successes must be integer-valued")

    } else if (NCOL(y) == 2) {
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(y - round(y)) > 1.0e-8))
        stop("Count data must be integer-valued")
      y <- round(y)
      nvec <- y[, 1] + y[, 2]
      y <- ifelse(nvec > 0, y[, 1] / nvec, 0)
      w <- w * nvec
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + nvec * y) / (1 + nvec)
    } else {
      stop("for the binomialff family, response 'y' must be a ",
           "vector of 0 and 1's\n",
           "or a factor ",
           "(first level = fail, other levels = success),\n",
           "or a 2-column matrix where col 1 is the no. of ",
           "successes and col 2 is the no. of failures")
    }


    if ( .imethod == 1)
      mustart <- (mustart + y) / 2


    extra$type.fitted <- .type.fitted
    extra$dimnamesy   <- dimnames(y)




    predictors.names <-
        c(namesof("prob"     , .lprob      , earg = .eprob      , tag = FALSE),
          namesof("onempstr0", .lonempstr0 , earg = .eonempstr0 , tag = FALSE))


    extra$w <- w  # Needed for @linkinv
    onemphi.init <- if (length( .ionempstr0 )) .ionempstr0 else {
        prob0.est <- sum(w[y == 0]) / sum(w)
        if ( .imethod == 1) {
          1 - (prob0.est - (1 - mustart)^w) / (1 - (1 - mustart)^w)
        } else {
          1 - prob0.est
        }
    }

    onemphi.init[onemphi.init <= -0.10] <- 0.10  # Lots of sample variation
    onemphi.init[onemphi.init <=  0.05] <- 0.15  # Last resort
    onemphi.init[onemphi.init >=  0.80] <- 0.80  # Last resort

    if ( length(mustart) && !length(etastart))
      mustart <- cbind(mustart,
                       rep(onemphi.init, len = n))  # 1st coln not a real mu

  }), list( .lonempstr0 = lonempstr0, .lprob = lprob,
            .eonempstr0 = eonempstr0, .eprob = eprob,
            .ionempstr0 = ionempstr0,
            .type.fitted = type.fitted,          
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mubin     <- eta2theta(eta[, 1], .lprob      , earg = .eprob      )
    onempstr0 <- eta2theta(eta[, 2], .lonempstr0 , earg = .eonempstr0 )


    orig.w <- if (length(tmp3 <- extra$orig.w)) tmp3 else
              rep(1, len = nrow(eta))
    priorw <- extra$w
    nvec <- priorw / orig.w


    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "pstr0", "onempstr0"))[1]

    ans <- switch(type.fitted,
                  "mean"      = (onempstr0) * mubin,
                  "pobs0"     = 1 - onempstr0 + (onempstr0)*(1-mubin)^nvec,  # P(Y=0)
                  "pstr0"     = 1 - onempstr0,
                  "onempstr0" =     onempstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempstr0 = lonempstr0, .lprob = lprob,
           .eonempstr0 = eonempstr0, .eprob = eprob,
           .type.fitted = type.fitted ))),
  last = eval(substitute(expression({
    misc$link <-    c("prob" = .lprob , "onempstr0" = .lonempstr0 )

    misc$earg <- list("prob" = .eprob , "onempstr0" = .eonempstr0 )

    misc$imethod <- .imethod


      misc$pobs0 <- phi + (1 - phi) * (1 - mubin)^w  # [1]  # P(Y=0)
      misc$pstr0 <- phi
  }), list( .lonempstr0 = lonempstr0, .lprob = lprob,
            .eonempstr0 = eonempstr0, .eprob = eprob,
            .imethod = imethod ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    cbind(theta2eta(mu[, 1], .lprob      , earg = .eprob      ),
          theta2eta(mu[, 2], .lonempstr0 , earg = .eonempstr0 ))
  }, list( .lonempstr0 = lonempstr0, .lprob = lprob,
           .eonempstr0 = eonempstr0, .eprob = eprob ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mubin     <- eta2theta(eta[, 1], .lprob      , earg = .eprob      )
    onempstr0 <- eta2theta(eta[, 2], .lonempstr0 , earg = .eonempstr0 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        dzibinom(x = round(w * y), size = w, prob = mubin,
                 log = TRUE, pstr0 = 1 - onempstr0)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempstr0 = lonempstr0, .lprob = lprob,
           .eonempstr0 = eonempstr0, .eprob = eprob ))),
  vfamily = c("zibinomialff"),
  deriv = eval(substitute(expression({
    mubin     <- eta2theta(eta[, 1], .lprob      , earg = .eprob      )
    onempstr0 <- eta2theta(eta[, 2], .lonempstr0 , earg = .eonempstr0 )
    omphi     <-     onempstr0
    phi       <- 1 - onempstr0


    prob0 <- (1 - mubin)^w  # Actually q^w
    pobs0 <- phi + (omphi) * prob0
    index <- (y == 0)
    dl.domphi <- -(1 - prob0) / pobs0  # Note "-"
    dl.domphi[!index] <- +1 / (omphi[!index])  # Note "+"

    dl.dmubin <- -w * (omphi) * (1 - mubin)^(w - 1) / pobs0
    dl.dmubin[!index] <- w[!index] *
        (    y[!index]  /      mubin[!index]   -
        (1 - y[!index]) / (1 - mubin[!index]))

    dmubin.deta <- dtheta.deta(mubin, .lprob      , earg = .eprob      )
    domphi.deta <- dtheta.deta(omphi, .lonempstr0 , earg = .eonempstr0 )

    ans <- cbind(dl.dmubin * dmubin.deta,
                 dl.domphi * domphi.deta)

      if ( .lprob == "logit") {
        ans[!index, 1] <- w[!index] * (y[!index] - mubin[!index])
      }

      ans
  }), list( .lonempstr0 = lonempstr0, .lprob = lprob,
            .eonempstr0 = eonempstr0, .eprob = eprob ))),
  weight = eval(substitute(expression({
    wz <- matrix(as.numeric(NA), nrow = n, ncol = dimm(M))



    ned2l.domphi2 <- (1 - prob0) / ((omphi) * pobs0)


    ned2l.domphimubin <- +w * ((1 - mubin)^(w - 1)) / pobs0  # Note "+"






    ned2l.dmubin2 <- (w * (omphi) / (mubin * (1 - mubin)^2)) *
                     (1 - mubin - w * mubin *
                     (1 - mubin)^w * phi / pobs0)





    wz[,iam(1, 1, M)] <- ned2l.dmubin2     * dmubin.deta^2
    wz[,iam(2, 2, M)] <- ned2l.domphi2     * domphi.deta^2
    wz[,iam(1, 2, M)] <- ned2l.domphimubin * domphi.deta * dmubin.deta
    if (TRUE) {
      ind6 <- (wz[, iam(1, 1, M)] < .Machine$double.eps)
      if (any(ind6))
        wz[ind6, iam(1, 1, M)] <- .Machine$double.eps
    }
    wz
  }), list( .lonempstr0 = lonempstr0, .lprob = lprob,
            .eonempstr0 = eonempstr0, .eprob = eprob ))))
}










dzibinom <- function(x, size, prob, pstr0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(size), length(prob), length(pstr0))
  if (length(x)     != LLL) x     <- rep(x,     len = LLL);
  if (length(size)  != LLL) size  <- rep(size,  len = LLL);
  if (length(prob)  != LLL) prob  <- rep(prob,  len = LLL);
  if (length(pstr0) != LLL) pstr0 <- rep(pstr0, len = LLL);

  ans <- dbinom(x = x, size = size, prob = prob, log = TRUE)


  ans <- if (log.arg) {
    ifelse(x == 0, log(pstr0 + (1-pstr0) * exp(ans)), log1p(-pstr0) + ans)
  } else {
    ifelse(x == 0,     pstr0 + (1-pstr0) * exp(ans) ,
                    (1-pstr0) * exp(ans))
  }


  prob0 <- (1 - prob)^size
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans
}


pzibinom <- function(q, size, prob, pstr0 = 0,
                    lower.tail = TRUE, log.p = FALSE) {

  LLL <- max(length(pstr0), length(size), length(prob), length(q))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL);

  ans <- pbinom(q, size, prob, lower.tail = lower.tail, log.p = log.p)
  ans <- ifelse(q < 0, 0, pstr0 + (1 - pstr0) * ans)


  prob0 <- (1 - prob)^size
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}


qzibinom <- function(p, size, prob, pstr0 = 0,
                    lower.tail = TRUE, log.p = FALSE) {
  LLL <- max(length(p), length(size), length(prob), length(pstr0))
  p     <- rep(p,     length = LLL)
  size  <- rep(size,  length = LLL)
  prob  <- rep(prob,  length = LLL)
  pstr0 <- rep(pstr0, length = LLL)


  ans <- p 
  ans[p <= pstr0] <- 0 
  ans[p >  pstr0] <-
    qbinom((p[p > pstr0] - pstr0[p > pstr0]) / (1 - pstr0[p > pstr0]),
           size[p > pstr0],
           prob[p > pstr0],
           lower.tail = lower.tail, log.p = log.p)



  prob0 <- (1 - prob)^size
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[p[ind0] <= pobs0] <- 0 
    pindex <- (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
    ans[pindex] <- qposbinom((p[pindex] - Pobs0) / (1 - Pobs0),
                             size = size[pindex],
                             prob = prob[pindex])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans
}


rzibinom <- function(n, size, prob, pstr0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  pstr0 <- rep(pstr0, len = use.n)
  size  <- rep(size,  len = use.n)
  prob  <- rep(prob,  len = use.n)

  ans <- rbinom(use.n, size, prob)
  ans[runif(use.n) < pstr0] <- 0



  prob0 <- (1 - prob)^size
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] <- rposbinom(sum(ind0), size = size[ind0], prob = prob[ind0])
    ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}












dzinegbin <- function(x, size, prob = NULL, munb = NULL, pstr0 = 0,
                     log = FALSE) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  LLL <- max(length(pstr0), length(size), length(prob), length(x))
  if (length(x)      != LLL) x      <- rep(x,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL);


  ans <- dnbinom(x = x, size = size, prob = prob, log = log.arg)

  ans <- if (log.arg)
    ifelse(x == 0, log(pstr0+(1-pstr0)*exp(ans)), log1p(-pstr0) + ans) else
    ifelse(x == 0,     pstr0+(1-pstr0)*    ans,       (1-pstr0) * ans)



  prob0 <- prob^size
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans
}


pzinegbin <- function(q, size, prob = NULL, munb = NULL, pstr0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  LLL <- max(length(pstr0), length(size), length(prob), length(q))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL);



  ans <- pnbinom(q = q, size = size, prob = prob)
  ans <- ifelse(q < 0, 0, pstr0 + (1 - pstr0) * ans)



  prob0 <- prob^size
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}


qzinegbin <- function(p, size, prob = NULL, munb = NULL, pstr0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size/(size + munb)
  }
  LLL <- max(length(p), length(prob), length(pstr0), length(size))
  if (length(p)     != LLL) p      <- rep(p,     len = LLL)
  if (length(pstr0) != LLL) pstr0  <- rep(pstr0, len = LLL);
  if (length(prob)  != LLL) prob   <- rep(prob,  len = LLL)
  if (length(size)  != LLL) size   <- rep(size,  len = LLL);

  ans <- p 
  ind4 <- (p > pstr0)
  ans[!ind4] <- 0
  ans[ ind4] <- qnbinom(p = (p[ind4] - pstr0[ind4]) / (1 - pstr0[ind4]),
                       size = size[ind4], prob = prob[ind4])



  prob0 <- prob^size
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[p[ind0] <= pobs0] <- 0 
    pindex <- (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
    ans[pindex] <- qposnegbin((p[pindex] - Pobs0) / (1 - Pobs0),
                              size = size[pindex],
                              prob = prob[pindex])
  }


  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN



  ans
}



rzinegbin <- function(n, size, prob = NULL, munb = NULL, pstr0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
               stop("bad input for argument 'n'") else n


  pstr0 <- rep(pstr0, len = use.n)
  size  <- rep(size,  len = use.n)
  prob  <- rep(prob,  len = use.n)


  ans <- rnbinom(n = use.n, size = size, prob = prob)
  ans <- ifelse(runif(use.n) < pstr0, rep(0, use.n), ans)



  prob0 <- rep(prob^size, len = use.n)
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0, na.rm = TRUE)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] <- rposnegbin(sum(ind0, na.rm = TRUE), size = size[ind0],
                    prob = prob[ind0])
    ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}








zinegbinomial.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 zinegbinomial <-
  function(lpstr0 = "logit", lmunb = "loge", lsize = "loge",
           type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
           ipstr0 = NULL,                    isize = NULL,
           zero = -3,  # 20130917; used to be c(-1, -3)
           imethod = 1, ishrinkage = 0.95,
           nsimEIM = 250) {


  lpstr0 <- as.list(substitute(lpstr0))
  epstr0 <- link2list(lpstr0)
  lpstr0 <- attr(epstr0, "function.name")

  lmunb <- as.list(substitute(lmunb))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")


  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "pstr0", "onempstr0"))[1]



  if (length(ipstr0) &&
     (!is.Numeric(ipstr0, positive = TRUE) ||
      any(ipstr0 >= 1)))
    stop("argument 'ipstr0' must contain values in (0,1)")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("argument 'isize' must contain positive values only")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1, 2 or 3")

  if (!is.Numeric(nsimEIM, length.arg = 1, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 50)
    warning("argument 'nsimEIM' should be greater than 50, say")

  if (!is.Numeric(ishrinkage, length.arg = 1) ||
      ishrinkage < 0 ||
      ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")




  new("vglmff",
  blurb = c("Zero-inflated negative binomial\n\n",
            "Links:    ",
            namesof("pstr0", lpstr0, earg = epstr0, tag = FALSE), ", ",
            namesof("munb",  lmunb,  earg = emunb,  tag = FALSE), ", ",
            namesof("size",  lsize,  earg = esize,  tag = FALSE), "\n",
            "Mean:     (1 - pstr0) * munb"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    M1 <- 3
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

      
  initialize = eval(substitute(expression({
    M1 <- 3

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted      <- .type.fitted
    extra$dimnamesy <- dimnames(y)


    
    mynames1 <- if (NOS == 1) "pstr0" else paste("pstr0", 1:NOS, sep = "")
    mynames2 <- if (NOS == 1) "munb"  else paste("munb",  1:NOS, sep = "")
    mynames3 <- if (NOS == 1) "size"  else paste("size",  1:NOS, sep = "")
    predictors.names <-
      c(namesof(mynames1, .lpstr0 , earg = .epstr0 , tag = FALSE),
        namesof(mynames2, .lmunb  , earg = .emunb  , tag = FALSE),
        namesof(mynames3, .lsize  , earg = .esize  , tag = FALSE))[
        interleave.VGAM(M1*NOS, M = M1)]

    if (!length(etastart)) {
      mum.init <- if ( .imethod == 3) {
        y + 1/16
      } else {
        mum.init <- y
        for (iii in 1:ncol(y)) {
          index <- (y[, iii] > 0)
          mum.init[, iii] <- if ( .imethod == 2)
              weighted.mean(y[index, iii], w     = w[index, iii]) else
                 median(rep(y[index, iii], times = w[index, iii])) + 1/8
        }
        (1 - .ishrinkage ) * (y + 1/16) + .ishrinkage * mum.init
      }


      pstr0.init <- if (length( .ipstr0 )) {
        matrix( .ipstr0 , n, ncoly, byrow = TRUE)
      } else {
        pstr0.init <- y
        for (iii in 1:ncol(y))
          pstr0.init[, iii] <- sum(w[y[, iii] == 0, iii]) / sum(w[, iii])
        pstr0.init[pstr0.init <= 0.02] <- 0.02 # Last resort
        pstr0.init[pstr0.init >= 0.98] <- 0.98 # Last resort
        pstr0.init
      }

        kay.init <-
        if ( is.Numeric( .isize )) {
          matrix( .isize, nrow = n, ncol = ncoly, byrow = TRUE)
        } else {
          zinegbin.Loglikfun <- function(kval, y, x, w, extraargs) {
            index0 <- (y == 0)
            pstr0vec <- extraargs$pstr0
            muvec <- extraargs$mu

            ans1 <- 0.0
            if (any( index0))
              ans1 <- ans1 + sum(w[ index0] *
                     dzinegbin(x = y[ index0], size = kval,
                               munb = muvec[ index0],
                               pstr0 = pstr0vec[ index0], log = TRUE))
            if (any(!index0))
              ans1 <- ans1 + sum(w[!index0] *
                     dzinegbin(x = y[!index0], size = kval,
                               munb = muvec[!index0],
                               pstr0 = pstr0vec[!index0], log = TRUE))
            ans1
          }
          k.grid <- 2^((-6):6)
          kay.init <- matrix(0, nrow = n, ncol = NOS)
          for (spp. in 1:NOS) {
            kay.init[, spp.] <-
              grid.search(k.grid, objfun = zinegbin.Loglikfun,
                          y = y[, spp.], x = x, w = w[, spp.],
                          extraargs = list(pstr0 = pstr0.init[, spp.],
                          mu  = mum.init[, spp.]))
          }
          kay.init
        }

        etastart <-
          cbind(theta2eta(pstr0.init, .lpstr0 , earg = .epstr0 ),
                theta2eta(mum.init,   .lmunb  , earg = .emunb  ),
                theta2eta(kay.init,   .lsize  , earg = .esize  ))
        etastart <-
          etastart[, interleave.VGAM(ncol(etastart), M = M1)]
    }
  }), list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
            .epstr0 = epstr0, .emunb = emunb, .esize = esize,
            .ipstr0 = ipstr0,                 .isize = isize,
            .type.fitted = type.fitted,
            .ishrinkage = ishrinkage,
            .imethod = imethod ))),
      
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "pstr0", "onempstr0"))[1]

    M1 <- 3
    NOS <- extra$NOS
    pstr0 <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                       .lpstr0 , earg = .epstr0 )
    if (type.fitted %in% c("mean", "pobs0"))
      munb  <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                         .lmunb  , earg = .emunb  )
    if (type.fitted %in% c("pobs0"))
      kmat  <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                        .lsize , earg = .esize )

    ans <- switch(type.fitted,
                  "mean"      = (1 - pstr0) * munb,
                  "pobs0"     = pstr0 + (1 - pstr0) *
                                (kmat / (kmat + munb))^kmat,  # P(Y=0)
                  "pstr0"     =     pstr0,
                  "onempstr0" = 1 - pstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))        
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpstr0 = lpstr0, .lsize = lsize, .lmunb = lmunb,
           .epstr0 = epstr0, .esize = esize, .emunb = emunb,
           .type.fitted = type.fitted ))),
      
  last = eval(substitute(expression({
    misc$link <-
      c(rep( .lpstr0 , length = NOS),
        rep( .lmunb  , length = NOS),
        rep( .lsize  , length = NOS))[interleave.VGAM(M1*NOS,
                                                      M = M1)]
    temp.names <-
      c(mynames1,
        mynames2,
        mynames3)[interleave.VGAM(M1*NOS, M = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M1*NOS)
    names(misc$earg) <- temp.names
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-2]] <- .epstr0
      misc$earg[[M1*ii-1]] <- .emunb
      misc$earg[[M1*ii  ]] <- .esize
    }

    misc$imethod <- .imethod
    misc$nsimEIM <- .nsimEIM
    misc$expected <- TRUE
    misc$M1 <- M1
    misc$ipstr0  <- .ipstr0
    misc$isize <- .isize
    misc$multipleResponses <- TRUE


  }), list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
            .epstr0 = epstr0, .emunb = emunb, .esize = esize,
            .ipstr0 = ipstr0,                 .isize = isize,
            .nsimEIM = nsimEIM, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    M1 <- 3
    NOS <- extra$NOS
    pstr0 <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                       .lpstr0 , earg = .epstr0 )
    munb  <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                       .lmunb , earg = .emunb )
    kmat  <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                       .lsize , earg = .esize )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzinegbin(x = y, size = kmat, munb = munb,
                         pstr0 = pstr0, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
           .epstr0 = epstr0, .emunb = emunb, .esize = esize ))),
  vfamily = c("zinegbinomial"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    pstr0 <- eta2theta(eta[, c(TRUE, FALSE, FALSE)],
                       .lpstr0 , earg = .epstr0 )
    munb  <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lmunb  , earg = .emunb  )
    kmat  <- eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lsize  , earg = .esize  )
    rzinegbin(nsim * length(munb),
              size = kmat, munb = munb, pstr0 = pstr0)
  }, list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
           .epstr0 = epstr0, .emunb = emunb, .esize = esize ))),






  deriv = eval(substitute(expression({
    M1 <- 3
    NOS <- extra$NOS

    pstr0 <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                      .lpstr0 , earg = .epstr0 )
    munb  <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                      .lmunb  , earg = .emunb  )
    kmat  <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                      .lsize  , earg = .esize  )

    dpstr0.deta <- dtheta.deta(pstr0, .lpstr0 , earg = .epstr0 )
    dmunb.deta  <- dtheta.deta(munb , .lmunb  , earg = .emunb  )
    dsize.deta  <- dtheta.deta(kmat , .lsize  , earg = .esize  )
    dthetas.detas <-
        (cbind(dpstr0.deta,
               dmunb.deta,
               dsize.deta))[, interleave.VGAM(M1*NOS, M = M1)]



    dl.dpstr0 <- -1 / (1 - pstr0)
    dl.dmunb <- y / munb - (y + kmat) / (munb + kmat)
    dl.dsize <- digamma(y + kmat) - digamma(kmat) -
               (y + kmat) / (munb + kmat) + 1 +
               log(kmat / (kmat + munb))



    for (spp. in 1:NOS) {
      index0 <- (y[, spp.] == 0)
      if (all(index0) || all(!index0))
        stop("must have some 0s AND some positive counts in the data")

      kmat.  <-  kmat[index0, spp.]
      munb.  <-  munb[index0, spp.]
      pstr0. <- pstr0[index0, spp.]


      tempk. <- kmat. / (kmat. + munb.)
      tempm. <- munb. / (kmat. + munb.)
      prob0. <- tempk.^kmat.
      df0.dmunb.  <- -tempk.* prob0.
      df0.dkmat.  <- prob0. * (tempm. + log(tempk.))

      denom. <- pstr0. + (1 - pstr0.) * prob0.
     dl.dpstr0[index0, spp.]  <- (1 - prob0.) / denom.
      dl.dmunb[index0, spp.]  <- (1 - pstr0.) * df0.dmunb. / denom.
      dl.dsize[index0, spp.]  <- (1 - pstr0.) * df0.dkmat. / denom.
    }  # of spp.


    dl.dthetas <-
      cbind(dl.dpstr0,
            dl.dmunb,
            dl.dsize)[, interleave.VGAM(M1*NOS, M = M1)]


      c(w) * dl.dthetas * dthetas.detas
  }), list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
            .epstr0 = epstr0, .emunb = emunb, .esize = esize ))),

  weight = eval(substitute(expression({



    wz <- matrix(0, n, M1*M - M1)

    ind3 <- iam(NA, NA, M = M1, both = TRUE, diag = TRUE)

    run.varcov <- array(0.0, c(n, length(ind3$row.index), NOS))

    for (ii in 1:( .nsimEIM )) {
      ysim <- rzinegbin(n = n*NOS, pstr0 = pstr0,
                        size = kmat, mu = munb)
      dim(ysim) <- c(n, NOS)
      index0 <- (ysim[, spp.] == 0)

      dl.dpstr0 <- -1 / (1 - pstr0)
      dl.dmunb <- ysim / munb - (ysim + kmat) / (munb + kmat)
      dl.dsize <- digamma(ysim + kmat) - digamma(kmat) -
                 (ysim + kmat) / (munb + kmat) + 1 +
                 log(kmat / (kmat + munb))


      for (spp. in 1:NOS) {
        index0 <- (ysim[, spp.] == 0)
        if (all(index0) || all(!index0)) {
          repeat {
            ysim[, spp.] <- rzinegbin(n = n,
                                      pstr0 = pstr0[, spp.],
                                      size  =  kmat[, spp.],
                                      mu    =  munb[, spp.])
            index0 <- (ysim[, spp.] == 0)
            if (any(!index0) && any(index0))
              break
          }
        }

        kmat.  <-  kmat[index0, spp.]
        munb.  <-  munb[index0, spp.]
        pstr0. <- pstr0[index0, spp.]


        tempk. <- kmat. / (kmat. + munb.)
        tempm. <- munb. / (kmat. + munb.)
        prob0.  <- tempk.^kmat.
        df0.dmunb.  <- -tempk.* prob0.
        df0.dkmat.  <- prob0. * (tempm. + log(tempk.))

        denom. <- pstr0. + (1 - pstr0.) * prob0.
       dl.dpstr0[index0, spp.] <- (1 - prob0.) / denom.
        dl.dmunb[index0, spp.] <- (1 - pstr0.) * df0.dmunb. / denom.
        dl.dsize[index0, spp.] <- (1 - pstr0.) * df0.dkmat. / denom.


        sdl.dthetas <- cbind(dl.dpstr0[, spp.],
                             dl.dmunb[, spp.],
                             dl.dsize[, spp.])

        temp3 <- sdl.dthetas
        run.varcov[,, spp.] <- run.varcov[,, spp.] +
                              temp3[, ind3$row.index] *
                              temp3[, ind3$col.index]


      }  # End of for (spp.) loop
    }  # End of ii nsimEIM loop

    run.varcov <- run.varcov / .nsimEIM

    wz1 <- if (intercept.only) {
      for (spp. in 1:NOS) {
        for (jay in 1:length(ind3$row.index)) {
          run.varcov[, jay, spp.] <- mean(run.varcov[, jay, spp.])
        }
      }
      run.varcov
    } else {
      run.varcov
    }

    for (spp. in 1:NOS) {
      wz1[,, spp.] <- wz1[,, spp.] *
                      dthetas.detas[, M1 * (spp. - 1) + ind3$row] *
                      dthetas.detas[, M1 * (spp. - 1) + ind3$col]
    }

    for (spp. in 1:NOS) {
      for (jay in 1:M1) {
        for (kay in jay:M1) {
          cptr <- iam((spp. - 1) * M1 + jay,
                     (spp. - 1) * M1 + kay, M = M)
          temp.wz1 <- wz1[,, spp.]
          wz[, cptr] <- temp.wz1[, iam(jay, kay, M = M1)]
        }
      }
    }


    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / M1)
  }), list( .lpstr0 = lpstr0,
            .epstr0 = epstr0, .nsimEIM = nsimEIM ))))
}  # End of zinegbinomial








zinegbinomialff.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 zinegbinomialff <-
  function(lmunb = "loge", lsize = "loge", lonempstr0 = "logit", 
           type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
           isize = NULL, ionempstr0 = NULL,  
           zero = c(-2, -3),
           imethod = 1, ishrinkage = 0.95,
           nsimEIM = 250) {



  lmunb <- as.list(substitute(lmunb))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")

  lonempstr0 <- as.list(substitute(lonempstr0))
  eonempstr0 <- link2list(lonempstr0)
  lonempstr0 <- attr(eonempstr0, "function.name")


  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "pstr0", "onempstr0"))[1]



  if (length(ionempstr0) &&
     (!is.Numeric(ionempstr0, positive = TRUE) ||
      any(ionempstr0 >= 1)))
    stop("argument 'ionempstr0' must contain values in (0,1)")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("argument 'isize' must contain positive values only")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1, 2 or 3")

  if (!is.Numeric(nsimEIM, length.arg = 1, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 50)
    warning("argument 'nsimEIM' should be greater than 50, say")

  if (!is.Numeric(ishrinkage, length.arg = 1) ||
      ishrinkage < 0 ||
      ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")




  new("vglmff",
  blurb = c("Zero-inflated negative binomial\n\n",
            "Links:    ",
            namesof("munb",  lmunb,  earg = emunb,  tag = FALSE), ", ",
            namesof("size",  lsize,  earg = esize,  tag = FALSE), ", ",
            namesof("onempstr0", lonempstr0, earg = eonempstr0, tag = FALSE),
            "\n",
            "Mean:     (1 - pstr0) * munb"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    M1 <- 3
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

      
  initialize = eval(substitute(expression({
    M1 <- 3

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- .type.fitted
    extra$dimnamesy   <- dimnames(y)


    
    mynames1 <- if (NOS == 1) "munb"  else paste("munb",  1:NOS, sep = "")
    mynames2 <- if (NOS == 1) "size"  else paste("size",  1:NOS, sep = "")
    mynames3 <- if (NOS == 1) "onempstr0" else paste("onempstr0", 1:NOS,
                                                     sep = "")
    predictors.names <-
      c(namesof(mynames1, .lmunb  , earg = .emunb  , tag = FALSE),
        namesof(mynames2, .lsize  , earg = .esize  , tag = FALSE),
        namesof(mynames3, .lonempstr0 , earg = .eonempstr0 , tag = FALSE))[
        interleave.VGAM(M1*NOS, M = M1)]

    if (!length(etastart)) {
      mum.init <- if ( .imethod == 3) {
        y + 1/16
      } else {
        mum.init <- y
        for (iii in 1:ncol(y)) {
          index <- (y[, iii] > 0)
          mum.init[, iii] <- if ( .imethod == 2)
              weighted.mean(y[index, iii], w     = w[index, iii]) else
                 median(rep(y[index, iii], times = w[index, iii])) + 1/8
        }
        (1 - .ishrinkage ) * (y + 1/16) + .ishrinkage * mum.init
      }


      onempstr0.init <- if (length( .ionempstr0 )) {
        matrix( .ionempstr0 , n, ncoly, byrow = TRUE)
      } else {
        pstr0.init <- y
        for (iii in 1:ncol(y))
          pstr0.init[, iii] <- sum(w[y[, iii] == 0, iii]) / sum(w[, iii])
        pstr0.init[pstr0.init <= 0.02] <- 0.02  # Last resort
        pstr0.init[pstr0.init >= 0.98] <- 0.98  # Last resort
        1 - pstr0.init
      }

        kay.init <-
        if ( is.Numeric( .isize )) {
          matrix( .isize, nrow = n, ncol = ncoly, byrow = TRUE)
        } else {
          zinegbin.Loglikfun <- function(kval, y, x, w, extraargs) {
            index0 <- (y == 0)
            pstr0vec <- extraargs$pstr0
            muvec <- extraargs$mu

            ans1 <- 0.0
            if (any( index0))
              ans1 <- ans1 + sum(w[ index0] *
                     dzinegbin(x = y[ index0], size = kval,
                               munb = muvec[ index0],
                               pstr0 = pstr0vec[ index0], log = TRUE))
            if (any(!index0))
              ans1 <- ans1 + sum(w[!index0] *
                     dzinegbin(x = y[!index0], size = kval,
                               munb = muvec[!index0],
                               pstr0 = pstr0vec[!index0], log = TRUE))
            ans1
          }
          k.grid <- 2^((-6):6)
          kay.init <- matrix(0, nrow = n, ncol = NOS)
          for (spp. in 1:NOS) {
            kay.init[, spp.] <-
              grid.search(k.grid, objfun = zinegbin.Loglikfun,
                          y = y[, spp.], x = x, w = w[, spp.],
                         extraargs = list(pstr0 = 1 - onempstr0.init[, spp.],
                                          mu    = mum.init[, spp.]))
          }
          kay.init
        }

        etastart <-
          cbind(theta2eta(mum.init,   .lmunb  , earg = .emunb  ),
                theta2eta(kay.init,   .lsize  , earg = .esize  ),
                theta2eta(onempstr0.init, .lonempstr0 ,
                          earg = .eonempstr0 ))
        etastart <-
          etastart[, interleave.VGAM(ncol(etastart), M = M1)]
    }
  }), list( .lonempstr0 = lonempstr0, .lmunb = lmunb, .lsize = lsize,
            .eonempstr0 = eonempstr0, .emunb = emunb, .esize = esize,
            .ionempstr0 = ionempstr0,                 .isize = isize,
            .type.fitted = type.fitted,
            .ishrinkage = ishrinkage,
            .imethod = imethod ))),
      
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "pstr0", "onempstr0"))[1]

    M1 <- 3
    NOS <- extra$NOS
    if (type.fitted %in% c("mean", "pobs0"))
      munb    <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                           .lmunb  , earg = .emunb  )
    if (type.fitted %in% c("pobs0"))
      kmat    <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                           .lsize , earg = .esize )
    onempstr0 <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                           .lonempstr0 , earg = .eonempstr0 )

    ans <- switch(type.fitted,
                  "mean"      = (onempstr0) * munb,
                  "pobs0"     = 1 -  onempstr0 + (onempstr0) *
                                (kmat / (kmat + munb))^kmat,  # P(Y=0)
                  "pstr0"     = 1 - onempstr0,
                  "onempstr0" =     onempstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))        
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempstr0 = lonempstr0, .lsize = lsize, .lmunb = lmunb,
           .eonempstr0 = eonempstr0, .esize = esize, .emunb = emunb,
           .type.fitted = type.fitted ))),
      
  last = eval(substitute(expression({
    misc$link <-
      c(rep( .lmunb      , length = NOS),
        rep( .lsize      , length = NOS),
        rep( .lonempstr0 , length = NOS))[interleave.VGAM(M1*NOS,
                                                          M = M1)]
    temp.names <-
      c(mynames1,
        mynames2,
        mynames3)[interleave.VGAM(M1*NOS, M = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M1*NOS)
    names(misc$earg) <- temp.names
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-2]] <- .emunb
      misc$earg[[M1*ii-1]] <- .esize
      misc$earg[[M1*ii  ]] <- .eonempstr0
    }

    misc$imethod <- .imethod
    misc$nsimEIM <- .nsimEIM
    misc$expected <- TRUE
    misc$M1 <- M1
    misc$ionempstr0  <- .ionempstr0
    misc$isize <- .isize
    misc$multipleResponses <- TRUE

  }), list( .lonempstr0 = lonempstr0, .lmunb = lmunb, .lsize = lsize,
            .eonempstr0 = eonempstr0, .emunb = emunb, .esize = esize,
            .ionempstr0 = ionempstr0,                 .isize = isize,
            .nsimEIM = nsimEIM, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    M1 <- 3
    NOS <- extra$NOS
    munb      <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                           .lmunb , earg = .emunb )
    kmat      <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                           .lsize , earg = .esize )
    onempstr0 <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                           .lonempstr0 , earg = .eonempstr0 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzinegbin(x = y, size = kmat, munb = munb,
                         pstr0 = 1 - onempstr0, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempstr0 = lonempstr0, .lmunb = lmunb, .lsize = lsize,
           .eonempstr0 = eonempstr0, .emunb = emunb, .esize = esize ))),
  vfamily = c("zinegbinomialff"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    munb <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lsize , earg = .esize )
    onempstr0 <- eta2theta(eta[, c(FALSE, FALSE, TRUE)],
                       .lpstr0 , earg = .epstr0 )
    rzinegbin(nsim * length(munb),
              size = kmat, munb = munb, pstr0 = 1 - onempstr0)
  }, list( .lonempstr0 = lonempstr0, .lmunb = lmunb, .lsize = lsize,
           .eonempstr0 = eonempstr0, .emunb = emunb, .esize = esize ))),




  deriv = eval(substitute(expression({
    M1 <- 3
    NOS <- extra$NOS

    munb      <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                           .lmunb  , earg = .emunb  )
    kmat      <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                           .lsize  , earg = .esize  )
    onempstr0 <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                           .lonempstr0 , earg = .eonempstr0 )

    donempstr0.deta <- dtheta.deta(onempstr0, .lonempstr0 ,
                                   earg = .eonempstr0 )
    dmunb.deta  <- dtheta.deta(munb , .lmunb  , earg = .emunb  )
    dsize.deta  <- dtheta.deta(kmat , .lsize  , earg = .esize  )
    dthetas.detas <-
        (cbind(dmunb.deta,
               dsize.deta,
               donempstr0.deta))[, interleave.VGAM(M1*NOS,
                                                   M = M1)]



    dl.dmunb <- y / munb - (y + kmat) / (munb + kmat)
    dl.dsize <- digamma(y + kmat) - digamma(kmat) -
               (y + kmat) / (munb + kmat) + 1 +
               log(kmat / (kmat + munb))
    dl.donempstr0 <- +1 / (onempstr0)



    for (spp. in 1:NOS) {
      index0 <- (y[, spp.] == 0)
      if (all(index0) || all(!index0))
        stop("must have some 0s AND some positive counts in the data")

      kmat.      <-      kmat[index0, spp.]
      munb.      <-      munb[index0, spp.]
      onempstr0. <- onempstr0[index0, spp.]


      tempk. <- kmat. / (kmat. + munb.)
      tempm. <- munb. / (kmat. + munb.)
      prob0. <- tempk.^kmat.
      df0.dmunb.  <- -tempk.* prob0.
      df0.dkmat.  <- prob0. * (tempm. + log(tempk.))

      denom. <- 1 - onempstr0. + (onempstr0.) * prob0.
     dl.donempstr0[index0, spp.]  <- -(1 - prob0.) / denom.  # note "-"
          dl.dmunb[index0, spp.]  <- (onempstr0.) * df0.dmunb. / denom.
          dl.dsize[index0, spp.]  <- (onempstr0.) * df0.dkmat. / denom.
    }  # of spp.


    dl.dthetas <-
      cbind(dl.dmunb,
            dl.dsize,
            dl.donempstr0)[, interleave.VGAM(M1*NOS, M = M1)]


      c(w) * dl.dthetas * dthetas.detas
  }), list( .lonempstr0 = lonempstr0, .lmunb = lmunb, .lsize = lsize,
            .eonempstr0 = eonempstr0, .emunb = emunb, .esize = esize ))),

  weight = eval(substitute(expression({



    wz <- matrix(0, n, M1*M - M1)

    ind3 <- iam(NA, NA, M = M1, both = TRUE, diag = TRUE)

    run.varcov <- array(0.0, c(n, length(ind3$row.index), NOS))

    for (ii in 1:( .nsimEIM )) {
      ysim <- rzinegbin(n = n*NOS, pstr0 = 1 - onempstr0,
                        size = kmat, mu = munb)
      dim(ysim) <- c(n, NOS)
      index0 <- (ysim[, spp.] == 0)

      dl.dmunb <- ysim / munb - (ysim + kmat) / (munb + kmat)
      dl.dsize <- digamma(ysim + kmat) - digamma(kmat) -
                  (ysim + kmat) / (munb + kmat) + 1 +
                  log(kmat / (kmat + munb))
      dl.donempstr0 <- +1 / (onempstr0)


      for (spp. in 1:NOS) {
        index0 <- (ysim[, spp.] == 0)
        if (all(index0) || all(!index0)) {
          repeat {
            ysim[, spp.] <- rzinegbin(n = n,
                                      pstr0 = 1 - onempstr0[, spp.],
                                      size  =  kmat[, spp.],
                                      mu    =  munb[, spp.])
            index0 <- (ysim[, spp.] == 0)
            if (any(!index0) && any(index0))
              break
          }
        }

        munb.      <-      munb[index0, spp.]
        kmat.      <-      kmat[index0, spp.]
        onempstr0. <- onempstr0[index0, spp.]


        tempk. <- kmat. / (kmat. + munb.)
        tempm. <- munb. / (kmat. + munb.)
        prob0.  <- tempk.^kmat.
        df0.dmunb.  <- -tempk.* prob0.
        df0.dkmat.  <- prob0. * (tempm. + log(tempk.))

        denom. <- 1 - onempstr0. + (onempstr0.) * prob0.
       dl.donempstr0[index0, spp.] <- -(1 - prob0.) / denom.  # note "-"
        dl.dmunb[index0, spp.] <- (onempstr0.) * df0.dmunb. / denom.
        dl.dsize[index0, spp.] <- (onempstr0.) * df0.dkmat. / denom.


        sdl.dthetas <- cbind(dl.dmunb[, spp.],
                             dl.dsize[, spp.],
                             dl.donempstr0[, spp.])

        temp3 <- sdl.dthetas
        run.varcov[,, spp.] <- run.varcov[,, spp.] +
                              temp3[, ind3$row.index] *
                              temp3[, ind3$col.index]


      }  # End of for (spp.) loop
    }  # End of ii nsimEIM loop

    run.varcov <- run.varcov / .nsimEIM

    wz1 <- if (intercept.only) {
      for (spp. in 1:NOS) {
        for (jay in 1:length(ind3$row.index)) {
          run.varcov[, jay, spp.] <- mean(run.varcov[, jay, spp.])
        }
      }
      run.varcov
    } else {
      run.varcov
    }

    for (spp. in 1:NOS) {
      wz1[,, spp.] <- wz1[,, spp.] *
                     dthetas.detas[, M1 * (spp. - 1) + ind3$row] *
                     dthetas.detas[, M1 * (spp. - 1) + ind3$col]
    }

    for (spp. in 1:NOS) {
      for (jay in 1:M1) {
        for (kay in jay:M1) {
          cptr <- iam((spp. - 1) * M1 + jay,
                      (spp. - 1) * M1 + kay, M = M)
          temp.wz1 <- wz1[,, spp.]
          wz[, cptr] <- temp.wz1[, iam(jay, kay, M = M1)]
        }
      }
    }


    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / M1)
  }), list( .lonempstr0 = lonempstr0,
            .eonempstr0 = eonempstr0, .nsimEIM = nsimEIM ))))
}  # End of zinegbinomialff










 zipoissonff <-
  function(llambda = "loge", lonempstr0 = "logit",
           type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
           ilambda = NULL,   ionempstr0 = NULL, imethod = 1,
           ishrinkage = 0.8, zero = -2) {


  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "pstr0", "onempstr0"))[1]



  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  lonempstr0 <- as.list(substitute(lonempstr0))
  eonempstr0 <- link2list(lonempstr0)
  lonempstr0 <- attr(eonempstr0, "function.name")



  if (length(ilambda))
    if (!is.Numeric(ilambda, positive = TRUE))
      stop("'ilambda' values must be positive")
  if (length(ionempstr0))
    if (!is.Numeric(ionempstr0, positive = TRUE) ||
      any(ionempstr0 >= 1))
      stop("'ionempstr0' values must be inside the interval (0,1)")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
    imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.Numeric(ishrinkage, length.arg = 1) ||
    ishrinkage < 0 ||
    ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")



  new("vglmff",
  blurb = c("Zero-inflated Poisson\n\n",
            "Links:    ",
            namesof("lambda",    llambda,    earg = elambda), ", ",
            namesof("onempstr0", lonempstr0, earg = eonempstr0), "\n",
            "Mean:     onempstr0 * lambda"),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    ncoly <- ncol(y)
    M1 <- 2
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly
    extra$type.fitted      <- .type.fitted
    extra$dimnamesy <- dimnames(y)


    mynames1 <- paste("lambda",    if (ncoly > 1) 1:ncoly else "", sep = "")
    mynames2 <- paste("onempstr0", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
      c(namesof(mynames1, .llambda    , earg = .elambda    , tag = FALSE),
        namesof(mynames2, .lonempstr0 , earg = .eonempstr0 , tag = FALSE))[
          interleave.VGAM(M, M = M1)]


      if (!length(etastart)) {

        matL <- matrix(if (length( .ilambda )) .ilambda else 0,
                       n, ncoly, byrow = TRUE)
        matP <- matrix(if (length( .ionempstr0 )) .ionempstr0 else 0,
                       n, ncoly, byrow = TRUE)

        for (jay in 1:ncoly) {
          yjay <- y[, jay]

          Phi0.init <- 1 - 0.85 * sum(w[yjay > 0]) / sum(w)
          Phi0.init[Phi0.init <= 0.02] <- 0.02  # Last resort
          Phi0.init[Phi0.init >= 0.98] <- 0.98  # Last resort

          if ( length(mustart)) {
            mustart <- matrix(mustart, n, ncoly)  # Make sure right size
            Lambda.init <- mustart / (1 - Phi0.init)
          } else if ( .imethod == 2) {
            mymean <- weighted.mean(yjay[yjay > 0],
                                       w[yjay > 0]) + 1/16
            Lambda.init <- (1 - .ishrinkage ) * (yjay + 1/8) + .ishrinkage * mymean
          } else {
            use.this <- median(yjay[yjay > 0]) + 1 / 16
            Lambda.init <- (1 - .ishrinkage ) * (yjay + 1/8) + .ishrinkage * use.this
          }

          zipois.Loglikfun <- function(phival, y, x, w, extraargs) {
            sum(c(w) * dzipois(x = y, pstr0 = phival,
                            lambda = extraargs$lambda,
                            log = TRUE))
          }
          phi0.grid <- seq(0.02, 0.98, len = 21)
          Phi0mat.init <- grid.search(phi0.grid,
                                      objfun = zipois.Loglikfun,
                                      y = y, x = x, w = w,
                                      extraargs = list(lambda = Lambda.init))
          if (length(mustart)) {
            Lambda.init <- Lambda.init / (1 - Phi0mat.init)
          }

        if (!length( .ilambda ))
          matL[, jay] <- Lambda.init
        if (!length( .ionempstr0 ))
          matP[, jay] <- Phi0mat.init
      }

      etastart <-
        cbind(theta2eta(    matL, .llambda    , earg = .elambda    ),
              theta2eta(1 - matP, .lonempstr0 , earg = .eonempstr0 ))[,
                        interleave.VGAM(M, M = M1)]

      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .lonempstr0 = lonempstr0, .llambda = llambda,
            .eonempstr0 = eonempstr0, .elambda = elambda,
            .ionempstr0 = ionempstr0, .ilambda = ilambda,
            .type.fitted = type.fitted,
            .imethod = imethod, .ishrinkage = ishrinkage ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "pstr0", "onempstr0"))[1]

    M1 <- 2
    ncoly <- extra$ncoly
    lambda    <- eta2theta(eta[, M1*(1:ncoly) - 1], .llambda ,
                           earg = .elambda )
    onempstr0 <- eta2theta(eta[, M1*(1:ncoly)    ], .lonempstr0 ,
                           earg = .eonempstr0 )


    ans <- switch(type.fitted,
                  "mean"      = onempstr0 * lambda,
                  "pobs0"     = 1 + onempstr0 * expm1(-lambda),  # P(Y=0)
                  "pstr0"     = 1 - onempstr0,
                  "onempstr0" =     onempstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempstr0 = lonempstr0, .llambda = llambda,
           .eonempstr0 = eonempstr0, .elambda = elambda,
           .type.fitted = type.fitted ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <-
      c(rep( .llambda    , length = ncoly),
        rep( .lonempstr0 , length = ncoly))[interleave.VGAM(M, M = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M = M1)]
    names(misc$link) <- temp.names


    misc$earg <- vector("list", M1 * ncoly)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .elambda
      misc$earg[[M1*ii  ]] <- .eonempstr0
    }

    misc$M1 <- M1
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE

      misc$pobs0 <- (1 - onempstr0) + onempstr0 * exp(-lambda)  # P(Y=0)
      misc$pobs0 <- as.matrix(misc$pobs0)
      if (length(dimnames(y)[[2]]) > 0)
        dimnames(misc$pobs0) <- dimnames(y)

      misc$pstr0 <- (1 - onempstr0)
      misc$pstr0 <- as.matrix(misc$pstr0)
      if (length(dimnames(y)[[2]]) > 0)
        dimnames(misc$pstr0) <- dimnames(y)
  }), list( .lonempstr0 = lonempstr0, .llambda = llambda,
            .eonempstr0 = eonempstr0, .elambda = elambda,
            .imethod = imethod ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    M1 <- 2
    ncoly <- extra$ncoly
    lambda    <- eta2theta(eta[, M1*(1:ncoly) - 1], .llambda    ,
                           earg = .elambda )
    onempstr0 <- eta2theta(eta[, M1*(1:ncoly)    ], .lonempstr0 ,
                           earg = .eonempstr0 )


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
                 dzipois(x = y, pstr0 = 1 - onempstr0, lambda = lambda,
                         log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempstr0 = lonempstr0, .llambda = llambda,
           .eonempstr0 = eonempstr0, .elambda = elambda ))),
  vfamily = c("zipoissonff"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    lambda    <- eta2theta(eta[, c(TRUE, FALSE)], .llambda ,
                           earg = .elambda    )
    onempstr0 <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                           earg = .eonempstr0 )
    rzipois(nsim * length(lambda), lambda = lambda, pstr0 = 1 - onempstr0)
  }, list( .lonempstr0 = lonempstr0, .llambda = llambda,
           .eonempstr0 = eonempstr0, .elambda = elambda ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    ncoly <- extra$ncoly
    lambda    <- eta2theta(eta[, M1*(1:ncoly) - 1], .llambda    ,
                           earg = .elambda )
    onempstr0 <- eta2theta(eta[, M1*(1:ncoly)    ], .lonempstr0 ,
                           earg = .eonempstr0 )


    dlambda.deta    <- dtheta.deta(lambda   , .llambda    ,
                                   earg = .elambda )
    donempstr0.deta <- dtheta.deta(onempstr0, .lonempstr0 ,
                                   earg = .eonempstr0 )

    denom <- 1 + onempstr0 * expm1(-lambda)
    ind0 <- (y == 0)
    dl.dlambda <- -onempstr0 * exp(-lambda) / denom
    dl.dlambda[!ind0] <- (y[!ind0] - lambda[!ind0]) / lambda[!ind0]
    dl.donempstr0 <- expm1(-lambda) / denom
    dl.donempstr0[!ind0] <- 1 / onempstr0[!ind0]

    ans <- c(w) * cbind(dl.dlambda    * dlambda.deta,
                        dl.donempstr0 * donempstr0.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M = M1)]


    if ( .llambda == "loge" && is.empty.list( .elambda ) &&
       any(lambda[!ind0] < .Machine$double.eps)) {
      for (spp. in 1:ncoly) {
        ans[!ind0[, spp.], M1 * spp.] <-
          w[!ind0[, spp.]] *
         (y[!ind0[, spp.], spp.] - lambda[!ind0[, spp.], spp.])
      }
    }



    ans
  }), list( .lonempstr0 = lonempstr0, .llambda = llambda,
            .eonempstr0 = eonempstr0, .elambda = elambda ))),
  weight = eval(substitute(expression({


    ned2l.dlambda2 <-  (    onempstr0) / lambda -
                    onempstr0 * (1 - onempstr0) * exp(-lambda) / denom
    ned2l.donempstr0.2 <- -expm1(-lambda) / ((onempstr0) * denom)
    ned2l.dphilambda <- +exp(-lambda) / denom


    wz <- array(c(c(w) * ned2l.dlambda2 * dlambda.deta^2,
                  c(w) * ned2l.donempstr0.2 * donempstr0.deta^2,
                  c(w) * ned2l.dphilambda * donempstr0.deta * dlambda.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)

    wz
  }), list( .llambda = llambda ))))
}







dzigeom <- function(x, prob, pstr0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(prob), length(pstr0))
  if (length(x)      != LLL) x      <- rep(x,      len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL);


  ans <- dgeom(x = x, prob = prob, log = TRUE)


  ans <- if (log.arg) {
    ifelse(x == 0, log(pstr0 + (1 - pstr0) * exp(ans)),
                   log1p(-pstr0) + ans)
  } else {
    ifelse(x == 0,     pstr0 + (1 - pstr0) * exp(ans) ,
                               (1 - pstr0) * exp(ans))
  }



  prob0 <- prob
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}



pzigeom <- function(q, prob, pstr0 = 0) {


  LLL <- max(length(q), length(prob), length(pstr0))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL);

  ans <- pgeom(q, prob)
  ans <- ifelse(q < 0, 0, pstr0 + (1-pstr0) * ans)


  prob0 <- prob
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}



qzigeom <- function(p, prob, pstr0 = 0) {
  LLL <- max(length(p), length(prob), length(pstr0))
  ans <- p <- rep(p,     len = LLL)
  prob     <- rep(prob,  len = LLL)
  pstr0    <- rep(pstr0, len = LLL)
  ans[p <= pstr0] <- 0 
  ind1 <- (p > pstr0)
  ans[ind1] <-
    qgeom((p[ind1] - pstr0[ind1]) / (1 - pstr0[ind1]),
          prob = prob[ind1])


  prob0 <- prob
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[p[ind0] <= pobs0] <- 0 
    pindex <- (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
    ans[pindex] <- 1 + qgeom((p[pindex] - Pobs0) / (1 - Pobs0),
                            prob = prob[pindex])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}



rzigeom <- function(n, prob, pstr0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n


  pstr0 <- rep(pstr0, len = use.n)
  prob  <- rep(prob,  len = use.n)


  ans <- rgeom(use.n, prob)
  ans[runif(use.n) < pstr0] <- 0


  prob0 <- prob
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] <- 1 + rgeom(sum(ind0), prob = prob[ind0])
    ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans
}




 zigeometric <-
  function(
           lpstr0 = "logit",
           lprob  = "logit",
           type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
           ipstr0  = NULL, iprob = NULL,
           imethod = 1,
           bias.red = 0.5,
           zero = NULL) {



  expected <- TRUE



  lpstr0 <- as.list(substitute(lpstr0))
  epstr0 <- link2list(lpstr0)
  lpstr0 <- attr(epstr0, "function.name")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "pstr0", "onempstr0"))[1]


  if (length(ipstr0))
    if (!is.Numeric(ipstr0, positive = TRUE) ||
        ipstr0 >= 1)
      stop("argument 'ipstr0' is out of range")

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")

  if (!is.Numeric(bias.red, length.arg = 1, positive = TRUE) ||
     bias.red > 1)
    stop("argument 'bias.red' must be between 0 and 1")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Zero-inflated geometric distribution,\n",
            "P[Y = 0] = pstr0 + (1 - pstr0) * prob,\n",
            "P[Y = y] = (1 - pstr0) * prob * (1 - prob)^y, ",
            "y = 1, 2, ...\n\n",
            "Link:     ",
            namesof("pstr0",  lpstr0,  earg = epstr0), ", ",
            namesof("prob",   lprob,   earg = eprob ), "\n",
            "Mean:     (1 - pstr0) * (1 - prob) / prob"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted  = type.fitted ))),
  initialize = eval(substitute(expression({

    M1 <- 2
    if (any(y < 0))
      stop("the response must not have negative values")

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted      <- .type.fitted
    extra$dimnamesy <- dimnames(y)


    mynames1 <- if (ncoly == 1) "pstr0" else
                paste("pstr0", 1:ncoly, sep = "")
    mynames2 <- if (ncoly == 1) "prob"  else
                paste("prob",  1:ncoly, sep = "")

    predictors.names <-
            c(namesof(mynames1, .lpstr0,  earg = .epstr0, tag = FALSE),
              namesof(mynames2, .lprob,   earg = .eprob,  tag = FALSE))[
          interleave.VGAM(M1 * NOS, M = M1)]


    if (!length(etastart)) {
      prob.init <- if ( .imethod == 3)
                       .bias.red / (1 + y + 1/8) else
                   if ( .imethod == 2)
                       .bias.red / (1 +
                   matrix(colMeans(y) + 1/8,
                          n, ncoly, byrow = TRUE)) else
                       .bias.red / (1 +
                   matrix(colSums(y * w) / colSums(w) + 1/8,
                          n, ncoly, byrow = TRUE))

      prob.init <- if (length( .iprob )) {
        matrix( .iprob , n, ncoly, byrow = TRUE)
      } else {
        prob.init # Already a matrix
      }


      prob0.est <- psze.init <- matrix(0, n, NOS)
      for (jlocal in 1:NOS) {
        prob0.est[, jlocal] <-
          sum(w[y[, jlocal] == 0, jlocal]) / sum(w[, jlocal])
        psze.init[, jlocal] <- if ( .imethod == 3)
                         prob0.est[, jlocal] / 2 else
                     if ( .imethod == 1)
                         pmax(0.05, (prob0.est[, jlocal] -
                                     median(prob.init[, jlocal]))) else
                         prob0.est[, jlocal] / 5
      }
      psze.init <- if (length( .ipstr0 )) {
        matrix( .ipstr0 , n, ncoly, byrow = TRUE)
      } else {
        psze.init # Already a matrix
      }



      etastart <-
        cbind(theta2eta(psze.init, .lpstr0, earg = .epstr0),
              theta2eta(prob.init, .lprob , earg = .eprob ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = M1)]
    }
  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .eprob = eprob, .epstr0 = epstr0,
            .iprob = iprob, .ipstr0 = ipstr0,
            .type.fitted = type.fitted,
            .bias.red = bias.red,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    pstr0  <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr0 , earg = .epstr0 )
    prob   <- eta2theta(eta[, c(FALSE, TRUE)], .lprob  , earg = .eprob  )

    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "pstr0", "onempstr0"))[1]

    ans <- switch(type.fitted,
                  "mean"      = (1 - pstr0) * (1 - prob) / prob,
                  "pobs0"     = pstr0 + (1 - pstr0) * prob,  # P(Y=0)
                  "pstr0"     =     pstr0,
                  "onempstr0" = 1 - pstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lprob = lprob, .lpstr0 = lpstr0,
           .eprob = eprob, .epstr0 = epstr0,
           .type.fitted = type.fitted ))),
  last = eval(substitute(expression({
    temp.names <- c(rep( .lpstr0 , len = NOS),
                    rep( .lprob  , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M = M1)]
    misc$link  <- temp.names


    misc$earg <- vector("list", M1 * NOS)
    names(misc$link) <-
    names(misc$earg) <-
        c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M = M1)]

    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epstr0
      misc$earg[[M1*ii  ]] <- .eprob
    }


    misc$imethod <- .imethod
    misc$zero <- .zero
    misc$bias.red <- .bias.red
    misc$expected <- .expected
    misc$ipstr0 <- .ipstr0
    misc$type.fitted <- .type.fitted


    misc$pobs0 <- pobs0 
    if (length(dimnames(y)[[2]]) > 0)
      dimnames(misc$pobs0) <- dimnames(y)
    misc$pstr0 <- pstr0
    if (length(dimnames(y)[[2]]) > 0)
      dimnames(misc$pstr0) <- dimnames(y)
  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .eprob = eprob, .epstr0 = epstr0,
                            .ipstr0 = ipstr0,
            .zero = zero,
            .expected = expected,
            .type.fitted = type.fitted,
            .bias.red = bias.red,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pstr0  <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr0 , earg = .epstr0 )
    prob   <- eta2theta(eta[, c(FALSE, TRUE)], .lprob  , earg = .eprob  )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzigeom(x = y, prob = prob, pstr0 = pstr0, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .lpstr0 = lpstr0,
           .eprob = eprob, .epstr0 = epstr0 ))),
  vfamily = c("zigeometric"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    pstr0  <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr0 , earg = .epstr0 )
    prob   <- eta2theta(eta[, c(FALSE, TRUE)], .lprob  , earg = .eprob  )
    rzigeom(nsim * length(pstr0), prob = prob, pstr0 = pstr0)
  }, list( .lprob = lprob, .lpstr0 = lpstr0,
           .eprob = eprob, .epstr0 = epstr0 ))),




  deriv = eval(substitute(expression({
    M1 <- 2
    pstr0  <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr0 , earg = .epstr0 )
    prob   <- eta2theta(eta[, c(FALSE, TRUE)], .lprob  , earg = .eprob  )


    prob0 <- prob  # P(Y == 0) from parent distribution, aka f(0)
    pobs0 <- pstr0 + (1 - pstr0) * prob0  # P(Y == 0)
    index0 <- (y == 0)

    dl.dpstr0 <- (1 - prob0) / pobs0
    dl.dpstr0[!index0] <- -1 / (1 - pstr0[!index0])

    dl.dprob <- (1 - pstr0) / pobs0
    dl.dprob[!index0]   <- 1 / prob[!index0] -
                           y[!index0] / (1 - prob[!index0])

    dpstr0.deta  <- dtheta.deta(pstr0 , .lpstr0 , earg = .epstr0 )
    dprob.deta   <- dtheta.deta(prob,   .lprob  , earg = .eprob  )

    dl.deta12 <- c(w) * cbind(dl.dpstr0 * dpstr0.deta,
                              dl.dprob  * dprob.deta)

    dl.deta12 <- dl.deta12[, interleave.VGAM(ncol(dl.deta12), M = M1)]
    dl.deta12
  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .eprob = eprob, .epstr0 = epstr0 ))),
  weight = eval(substitute(expression({
    if ( .expected ) {


      ned2l.dprob2 <- (1 - pstr0)^2 / pobs0 +
                      (1 - pstr0) * ((1 - prob) / prob) *
                                    (1 / prob + 1 / (1 - prob)^2)


      ned2l.dpstr0.prob <- 1 / pobs0
      ned2l.dpstr02 <- (1 - prob0) / ((1 - pstr0) * pobs0)
    } else {
      od2l.dprob2 <- ((1 - pstr0) / pobs0)^2
      od2l.dprob2[!index0] <- 1 / (prob[!index0])^2 +
                              y[!index0] / (1 - prob[!index0])^2
      od2l.dpstr0.prob <- (pobs0 + (1 - prob0) * (1 - pstr0)) / pobs0^2
      od2l.dpstr0.prob[!index0] <- 0

      od2l.dpstr02 <- ((1 - prob0) / pobs0)^2
      od2l.dpstr02[!index0] <- 1 / (1 - pstr0[!index0])^2
    }


    allvals <- if ( .expected )
                 c(c(w) * ned2l.dpstr02 * dpstr0.deta^2,
                   c(w) * ned2l.dprob2  *  dprob.deta^2,
                   c(w) * ned2l.dpstr0.prob * dprob.deta * dpstr0.deta) else
                 c(c(w) *  od2l.dpstr02 * dpstr0.deta^2,
                   c(w) *  od2l.dprob2  *  dprob.deta^2,
                   c(w) *  od2l.dpstr0.prob * dprob.deta * dpstr0.deta)
    wz <- array(allvals, dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)


    wz
  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .eprob = eprob, .epstr0 = epstr0,
            .expected = expected ))))
}




 zigeometricff <-
  function(lprob       = "logit",
           lonempstr0  = "logit",
           type.fitted = c("mean", "pobs0", "pstr0", "onempstr0"),
           iprob = NULL,   ionempstr0  = NULL,
           imethod = 1,
           bias.red = 0.5,
           zero = -2) {


  expected <- TRUE



  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lonempstr0 <- as.list(substitute(lonempstr0))
  eonempstr0 <- link2list(lonempstr0)
  lonempstr0 <- attr(eonempstr0, "function.name")


  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "pstr0", "onempstr0"))[1]


  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")

  if (length(ionempstr0))
    if (!is.Numeric(ionempstr0, positive = TRUE) ||
        ionempstr0 >= 1)
      stop("argument 'ionempstr0' is out of range")

  if (!is.Numeric(bias.red, length.arg = 1, positive = TRUE) ||
     bias.red > 1)
    stop("argument 'bias.red' must be between 0 and 1")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Zero-inflated geometric distribution,\n",
            "P[Y = 0] = 1 - onempstr0 + onempstr0 * prob,\n",
            "P[Y = y] = onempstr0 * prob * (1 - prob)^y, ",
            "y = 1, 2, ...\n\n",
            "Link:     ",
            namesof("prob",       lprob,       earg = eprob ), ", ",
            namesof("onempstr0",  lonempstr0,  earg = eonempstr0), "\n",
            "Mean:     onempstr0 * (1 - prob) / prob"),
  constraints = eval(substitute(expression({

    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted  = type.fitted ))),
  initialize = eval(substitute(expression({

    M1 <- 2
    if (any(y < 0))
      stop("the response must not have negative values")

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted      <- .type.fitted
    extra$dimnamesy <- dimnames(y)


    mynames1 <- if (ncoly == 1) "prob"      else
                paste("prob",      1:ncoly, sep = "")
    mynames2 <- if (ncoly == 1) "onempstr0" else
                paste("onempstr0", 1:ncoly, sep = "")

    predictors.names <-
      c(namesof(mynames1, .lprob      , earg = .eprob      , tag = FALSE),
        namesof(mynames2, .lonempstr0 , earg = .eonempstr0 , tag = FALSE))[
        interleave.VGAM(M1*NOS, M = M1)]


    if (!length(etastart)) {
      prob.init <- if ( .imethod == 3)
                       .bias.red / (1 + y + 1/8) else
                   if ( .imethod == 2)
                       .bias.red / (1 +
                   matrix(colMeans(y) + 1/8,
                          n, ncoly, byrow = TRUE)) else
                       .bias.red / (1 +
                   matrix(colSums(y * w) / colSums(w) + 1/8,
                          n, ncoly, byrow = TRUE))

      prob.init <- if (length( .iprob )) {
        matrix( .iprob , n, ncoly, byrow = TRUE)
      } else {
        prob.init  # Already a matrix
      }


      prob0.est <- psze.init <- matrix(0, n, NOS)
      for (jlocal in 1:NOS) {
        prob0.est[, jlocal] <-
          sum(w[y[, jlocal] == 0, jlocal]) / sum(w[, jlocal])
        psze.init[, jlocal] <- if ( .imethod == 3)
                         prob0.est[, jlocal] / 2 else
                     if ( .imethod == 1)
                         pmax(0.05, (prob0.est[, jlocal] -
                                     median(prob.init[, jlocal]))) else
                         prob0.est[, jlocal] / 5
      }
      psze.init <- if (length( .ionempstr0 )) {
        matrix( 1 - .ionempstr0 , n, ncoly, byrow = TRUE)
      } else {
        psze.init # Already a matrix
      }



      etastart <-
        cbind(theta2eta(    prob.init, .lprob      , earg = .eprob      ),
              theta2eta(1 - psze.init, .lonempstr0 , earg = .eonempstr0 ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = M1)]
    }
  }), list( .lprob = lprob, .lonempstr0 = lonempstr0,
            .eprob = eprob, .eonempstr0 = eonempstr0,
            .iprob = iprob, .ionempstr0 = ionempstr0,
            .type.fitted = type.fitted,
            .bias.red = bias.red,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob      <- eta2theta(eta[, c(TRUE, FALSE)], .lprob      ,
                           earg = .eprob  )
    onempstr0 <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                           earg = .eonempstr0 )

    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "pstr0", "onempstr0"))[1]

    ans <- switch(type.fitted,
                  "mean"      = onempstr0 * (1 - prob) / prob,
                  "pobs0"     = 1 - onempstr0 + onempstr0 * prob,  # P(Y=0)
                  "pstr0"     = 1 - onempstr0,
                  "onempstr0" =     onempstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lprob = lprob, .lonempstr0 = lonempstr0,
           .eprob = eprob, .eonempstr0 = eonempstr0,
           .type.fitted = type.fitted ))),
  last = eval(substitute(expression({
    temp.names <- c(rep( .lprob  , len = NOS),
                    rep( .lonempstr0 , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M = M1)]
    misc$link  <- temp.names


    misc$earg <- vector("list", M1 * NOS)
    names(misc$link) <-
    names(misc$earg) <-
        c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M = M1)]

    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .eprob
      misc$earg[[M1*ii  ]] <- .eonempstr0
    }


    misc$imethod  <- .imethod
    misc$zero     <- .zero
    misc$bias.red <- .bias.red
    misc$expected <- .expected
    misc$ionempstr0   <- .ionempstr0


    misc$pobs0 <- pobs0 
    if (length(dimnames(y)[[2]]) > 0)
      dimnames(misc$pobs0) <- dimnames(y)
    misc$onempstr0 <- onempstr0
    if (length(dimnames(y)[[2]]) > 0)
      dimnames(misc$onempstr0) <- dimnames(y)
  }), list( .lprob = lprob, .lonempstr0 = lonempstr0,
            .eprob = eprob, .eonempstr0 = eonempstr0,
                            .ionempstr0 = ionempstr0,
            .zero = zero,
            .expected = expected,
            .bias.red = bias.red,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    prob       <- eta2theta(eta[, c(TRUE, FALSE)], .lprob      ,
                            earg = .eprob )
    onempstr0  <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                            earg = .eonempstr0 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzigeom(x = y, prob = prob, pstr0 = 1 - onempstr0,
                       log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .lonempstr0 = lonempstr0,
           .eprob = eprob, .eonempstr0 = eonempstr0 ))),
  vfamily = c("zigeometricff"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    prob       <- eta2theta(eta[, c(TRUE, FALSE)], .lprob      ,
                            earg = .eprob )
    onempstr0  <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                            earg = .eonempstr0 )
    rzigeom(nsim * length(onempstr0), prob = prob, pstr0 = 1 - onempstr0)
  }, list( .lprob = lprob, .lonempstr0 = lonempstr0,
           .eprob = eprob, .eonempstr0 = eonempstr0 ))),





  deriv = eval(substitute(expression({
    M1 <- 2
    prob      <- eta2theta(eta[, c(TRUE, FALSE)], .lprob      ,
                           earg = .eprob  )
    onempstr0 <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                           earg = .eonempstr0 )


    prob0 <- prob  # P(Y == 0) from the parent distribution
    pobs0 <- 1 - onempstr0 + (onempstr0) * prob0  # P(Y == 0)
    index0 <- (y == 0)


    dl.donempstr0 <- -(1 - prob0) / pobs0  # zz
    dl.donempstr0[!index0] <-  1 / (onempstr0[!index0])  # zz

    dl.dprob <- (onempstr0) / pobs0
    dl.dprob[!index0]   <- 1 / prob[!index0] -
                           y[!index0] / (1 - prob[!index0])

    dprob.deta       <- dtheta.deta(prob      , .lprob      ,
                                    earg = .eprob )
    donempstr0.deta  <- dtheta.deta(onempstr0 , .lonempstr0 ,
                                    earg = .eonempstr0 )

    dl.deta12 <- c(w) * cbind(dl.dprob      * dprob.deta,
                              dl.donempstr0 *  donempstr0.deta)

    dl.deta12 <- dl.deta12[, interleave.VGAM(ncol(dl.deta12), M = M1)]
    dl.deta12
  }), list( .lprob = lprob, .lonempstr0 = lonempstr0,
            .eprob = eprob, .eonempstr0 = eonempstr0 ))),
  weight = eval(substitute(expression({
    if ( .expected ) {

      ned2l.dprob2 <- (onempstr0)^2 / pobs0 +
                      (onempstr0) * ((1 - prob) / prob) *
                                    (1 / prob + 1 / (1 - prob)^2)


      ned2l.donempstr0.prob <- -1 / pobs0
      ned2l.donempstr02 <- (1 - prob0) / ((    onempstr0) * pobs0)
    } else {
      od2l.dprob2 <- ((    onempstr0) / pobs0)^2
      od2l.dprob2[!index0] <- 1 / (prob[!index0])^2 +
                              y[!index0] / (1 - prob[!index0])^2
      od2l.donempstr0.prob <- -(pobs0 + (1 - prob0) * (onempstr0)) / pobs0^2
      od2l.donempstr0.prob[!index0] <- 0

      od2l.donempstr02 <- ((1 - prob0) / pobs0)^2
      od2l.donempstr02[!index0] <- 1 / (    onempstr0[!index0])^2
    }


    allvals <- if ( .expected )
                 c(c(w) * ned2l.dprob2  *  dprob.deta^2,
                   c(w) * ned2l.donempstr02 * donempstr0.deta^2,
                   c(w) * ned2l.donempstr0.prob * dprob.deta *
                                                  donempstr0.deta) else
                 c(c(w) *  od2l.dprob2  *  dprob.deta^2,
                   c(w) *  od2l.donempstr02 * donempstr0.deta^2,
                   c(w) *  od2l.donempstr0.prob * dprob.deta *
                                                  donempstr0.deta)
    wz <- array(allvals, dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)


    wz
  }), list( .lprob = lprob, .lonempstr0 = lonempstr0,
            .eprob = eprob, .eonempstr0 = eonempstr0,
            .expected = expected ))))
}




dzageom <- function(x, prob, pobs0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(prob), length(pobs0))
  if (length(x)      != LLL) x      <- rep(x,     len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,  len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0, len = LLL);
  ans <- rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  index0 <- (x == 0)

  if (log.arg) {
    ans[ index0] <- log(pobs0[index0])
    ans[!index0] <- log1p(-pobs0[!index0]) +
                   dposgeom(x[!index0],
                            prob = prob[!index0], log = TRUE)
  } else {
    ans[ index0] <- pobs0[index0]
    ans[!index0] <- (1-pobs0[!index0]) *
                   dposgeom(x[!index0],
                            prob = prob[!index0])
  }
  ans
}



pzageom <- function(q, prob, pobs0 = 0) {

  LLL <- max(length(q), length(prob), length(pobs0))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL);
  ans <- rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans[q >  0] <- pobs0[q > 0] +
                (1 - pobs0[q > 0]) *
                pposgeom(q[q > 0], prob = prob[q > 0])
  ans[q <  0] <- 0
  ans[q == 0] <- pobs0[q == 0]
  ans
}


qzageom <- function(p, prob, pobs0 = 0) {

  LLL <- max(length(p), length(prob), length(pobs0))
  if (length(p)      != LLL) p      <- rep(p,      len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL);

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans <- p
  ind4 <- (p > pobs0)
  ans[!ind4] <- 0.0
  ans[ ind4] <- qposgeom((p[ind4] - pobs0[ind4]) / (1 - pobs0[ind4]),
                         prob = prob[ind4])
  ans
}


rzageom <- function(n, prob, pobs0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
               stop("bad input for argument 'n'") else n

  ans <- rposgeom(use.n, prob)
  if (length(pobs0) != use.n)
    pobs0 <- rep(pobs0, len = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ifelse(runif(use.n) < pobs0, 0, ans)
}










dzabinom <- function(x, size, prob, pobs0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(size), length(prob), length(pobs0))
  if (length(x)      != LLL) x      <- rep(x,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL);
  ans <- rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  index0 <- (x == 0)

  if (log.arg) {
    ans[ index0] <- log(pobs0[index0])
    ans[!index0] <- log1p(-pobs0[!index0]) +
                   dposbinom(x[!index0], size = size[!index0],
                             prob = prob[!index0], log = TRUE)
  } else {
    ans[ index0] <- pobs0[index0]
    ans[!index0] <- (1-pobs0[!index0]) *
                   dposbinom(x[!index0], size = size[!index0],
                             prob = prob[!index0])
  }
  ans
}



pzabinom <- function(q, size, prob, pobs0 = 0) {

  LLL <- max(length(q), length(size), length(prob), length(pobs0))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL);
  ans <- rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) ||
      any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans[q >  0] <- pobs0[q > 0] +
                (1 - pobs0[q > 0]) *
                pposbinom(q[q > 0], size = size[q > 0], prob = prob[q > 0])
  ans[q <  0] <- 0
  ans[q == 0] <- pobs0[q == 0]
  ans
}


qzabinom <- function(p, size, prob, pobs0 = 0) {

  LLL <- max(length(p), length(size), length(prob), length(pobs0))
  if (length(p)      != LLL) p      <- rep(p,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL);

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans <- p
  ind4 <- (p > pobs0)
  ans[!ind4] <- 0.0
  ans[ ind4] <- qposbinom((p[ind4] - pobs0[ind4]) / (1 - pobs0[ind4]),
                         size = size[ind4],
                         prob = prob[ind4])
  ans
}


rzabinom <- function(n, size, prob, pobs0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
               stop("bad input for argument 'n'") else n

  ans <- rposbinom(use.n, size, prob)
  if (length(pobs0) != use.n)
    pobs0 <- rep(pobs0, len = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ifelse(runif(use.n) < pobs0, 0, ans)
}






 zabinomial <-
  function(lpobs0 = "logit",
           lprob  = "logit",
           type.fitted = c("mean", "pobs0"),
           ipobs0 = NULL, iprob = NULL,
           imethod = 1,
           zero = NULL  # Was zero = 2 prior to 20130917
          ) {



  lpobs0 <- as.list(substitute(lpobs0))
  epobs0 <- link2list(lpobs0)
  lpobs0 <- attr(epobs0, "function.name")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")


  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0"))[1]

  if (length(ipobs0))
    if (!is.Numeric(ipobs0, positive = TRUE) ||
        ipobs0 >= 1)
      stop("argument 'ipobs0' is out of range")

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")



  if (!is.Numeric(imethod, length.arg = 1,
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
            namesof("pobs0",   lpobs0, earg = epobs0), ", ",
            namesof("prob" ,   lprob,  earg = eprob),  "\n",
            "Mean:     (1 - pobs0) * prob / (1 - (1 - prob)^size)"),
  constraints = eval(substitute(expression({
      constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted ))),

  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w



    if (NCOL(y) == 1) {
      if (is.factor(y))
        y <- y != levels(y)[1]
      nn <- rep(1, n)
      if (!all(y >= 0 & y <= 1))
        stop("response values must be in [0, 1]")
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + w * y) / (1.0 + w)


      no.successes <- y
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(no.successes - round(no.successes)) > 1.0e-8))
        stop("Number of successes must be integer-valued")

    } else if (NCOL(y) == 2) {
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(y - round(y)) > 1.0e-8))
        stop("Count data must be integer-valued")
      y <- round(y)
      nvec <- y[, 1] + y[, 2]
      y <- ifelse(nvec > 0, y[, 1] / nvec, 0)
      w <- w * nvec
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + nvec * y) / (1 + nvec)
    } else {
      stop("for the binomialff family, response 'y' must be a ",
           "vector of 0 and 1's\n",
           "or a factor ",
           "(first level = fail, other levels = success),\n",
           "or a 2-column matrix where col 1 is the no. of ",
           "successes and col 2 is the no. of failures")
    }
    if (!all(w == 1))
      extra$new.w <- w


    y <- as.matrix(y)
    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)

    extra$dimnamesy <- dimnames(y)
    extra$type.fitted      <- .type.fitted


    predictors.names <-
        c(namesof("pobs0", .lpobs0 , earg = .epobs0 , tag = FALSE),
          namesof("prob" , .lprob  , earg = .eprob  , tag = FALSE))
          


    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    phi.init <- if (length( .ipobs0 )) .ipobs0 else {
        prob0.est <- sum(Size[y == 0]) / sum(Size)
        if ( .imethod == 1) {
          (prob0.est - (1 - mustart)^Size) / (1 - (1 - mustart)^Size)
        } else
        if ( .imethod == 2) {
          prob0.est
        } else {
          prob0.est * 0.5
        }
    }

    phi.init[phi.init <= -0.10] <- 0.50  # Lots of sample variation
    phi.init[phi.init <=  0.01] <- 0.05  # Last resort
    phi.init[phi.init >=  0.99] <- 0.95  # Last resort




    if (!length(etastart)) {
      etastart <-
        cbind(theta2eta(phi.init, .lpobs0, earg = .epobs0 ),
              theta2eta( mustart, .lprob,  earg = .eprob  ))
              

      mustart <- NULL
    }
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0,
            .iprob = iprob, .ipobs0 = ipobs0,
            .imethod = imethod,
            .type.fitted = type.fitted ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0"))[1]
    
    phi0  <- eta2theta(eta[, 1], .lpobs0, earg = .epobs0 )
    prob  <- eta2theta(eta[, 2], .lprob,  earg = .eprob  )
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    ans <- switch(type.fitted,
                  "mean"      = (1 - phi0) * prob / (1 - (1 - prob)^Size),
                  "pobs0"     = phi0)  # P(Y=0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lprob = lprob, .lpobs0 = lpobs0,
           .eprob = eprob, .epobs0 = epobs0 ))),

  last = eval(substitute(expression({
    misc$link <-    c(prob = .lprob, pobs0 = .lpobs0 )
    misc$earg <- list(prob = .eprob, pobs0 = .epobs0 )

    misc$imethod  <- .imethod
    misc$zero     <- .zero
    misc$expected <- TRUE
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0,
            .zero = zero,
            .imethod = imethod ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w
    pobs0 <- eta2theta(eta[, 1], .lpobs0 , earg = .epobs0 )
    prob  <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        orig.w * dzabinom(x = round(y * Size), size = Size,
                          prob = prob, pobs0 = pobs0,
                          log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .lpobs0 = lpobs0,
           .eprob = eprob, .epobs0 = epobs0 ))),
  vfamily = c("zabinomial"),

  deriv = eval(substitute(expression({
    NOS <- if (length(extra$NOS)) extra$NOS else 1
    M1 <- 2

    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    phi0 <- eta2theta(eta[, 1], .lpobs0 , earg = .epobs0 )
    prob <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )

    dphi0.deta <- dtheta.deta(phi0, .lpobs0, earg = .epobs0 )
    dprob.deta <- dtheta.deta(prob, .lprob , earg = .eprob  )

    df0.dprob   <- -Size *              (1 -  prob)^(Size - 1)
    df02.dprob2 <-  Size * (Size - 1) * (1 -  prob)^(Size - 2)
    prob0  <- (1 -  prob)^(Size)
    oneminusf0  <- 1 - prob0


    dl.dphi0 <- -1 / (1 - phi0)
    dl.dprob <-  c(w)      * (y / prob - (1 - y) / (1 - prob)) +
                 c(orig.w) * df0.dprob / oneminusf0


    dl.dphi0[y == 0] <- 1 / phi0[y == 0]  # Do it in one line
    skip <- extra$skip.these
    for (spp. in 1:NOS) {
      dl.dprob[skip[, spp.], spp.] <- 0
    }


    ans <- cbind(c(orig.w) * dl.dphi0 * dphi0.deta,
                             dl.dprob * dprob.deta)
                 
                 
    ans
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0 ))),


  weight = eval(substitute(expression({
    wz <- matrix(0.0, n, M1)

    usualmeanY <-  prob
    meanY <- (1 - phi0) * usualmeanY / oneminusf0


    term1 <-  c(Size) * (meanY /      prob^2 -
                         meanY / (1 - prob)^2) +
             c(Size) * (1 - phi0) / (1 - prob)^2

    term2 <-  -(1 - phi0) * df02.dprob2 / oneminusf0
    term3 <-  -(1 - phi0) * (df0.dprob  / oneminusf0)^2
    ned2l.dprob2 <- term1 + term2 + term3
    wz[, iam(2, 2, M)] <- ned2l.dprob2 * dprob.deta^2


    mu.phi0 <- phi0
    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if ( .lpobs0 == "logit" && is.empty.list( .epobs0 )) {
      tmp100
    } else {
      (dphi0.deta^2) / tmp100
    }
    wz[, iam(1, 1, M)] <- tmp200


    c(orig.w) * wz
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0 ))))
}





 zabinomialff <-
  function(lprob  = "logit",
           lonempobs0 = "logit",
           type.fitted = c("mean", "pobs0", "onempobs0"),
           iprob = NULL, ionempobs0 = NULL,
           imethod = 1,
           zero = 2) {


  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lonempobs0 <- as.list(substitute(lonempobs0))
  eonempobs0 <- link2list(lonempobs0)
  lonempobs0 <- attr(eonempobs0, "function.name")


  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "onempobs0"))[1]

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")
  if (length(ionempobs0))
    if (!is.Numeric(ionempobs0, positive = TRUE) ||
        ionempobs0 >= 1)
      stop("argument 'ionempobs0' is out of range")



  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Zero-altered binomial distribution ",
            "(Bernoulli and positive-binomial conditional model)\n\n",
            "P[Y = 0] = 1 - onempobs0,\n",
            "P[Y = y] = onempobs0 * dposbinom(x = y, size, prob), ",
            "y = 1, 2, ..., size,\n\n",
            "Link:     ",
            namesof("prob"     , lprob     , earg = eprob     ), ", ",
            namesof("onempobs0", lonempobs0, earg = eonempobs0), "\n",
            "Mean:     onempobs0 * prob / (1 - (1 - prob)^size)"),
  constraints = eval(substitute(expression({
      constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted ))),

  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w



    if (NCOL(y) == 1) {
      if (is.factor(y))
        y <- y != levels(y)[1]
      nn <- rep(1, n)
      if (!all(y >= 0 & y <= 1))
        stop("response values must be in [0, 1]")
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + w * y) / (1.0 + w)


      no.successes <- y
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(no.successes - round(no.successes)) > 1.0e-8))
        stop("Number of successes must be integer-valued")

    } else if (NCOL(y) == 2) {
      if (min(y) < 0)
        stop("Negative data not allowed!")
      if (any(abs(y - round(y)) > 1.0e-8))
        stop("Count data must be integer-valued")
      y <- round(y)
      nvec <- y[, 1] + y[, 2]
      y <- ifelse(nvec > 0, y[, 1] / nvec, 0)
      w <- w * nvec
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + nvec * y) / (1 + nvec)
    } else {
      stop("for the binomialff family, response 'y' must be a ",
           "vector of 0 and 1's\n",
           "or a factor ",
           "(first level = fail, other levels = success),\n",
           "or a 2-column matrix where col 1 is the no. of ",
           "successes and col 2 is the no. of failures")
    }
    if (!all(w == 1))
      extra$new.w <- w


    y <- as.matrix(y)
    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)

    extra$dimnamesy   <- dimnames(y)
    extra$type.fitted <- .type.fitted


    predictors.names <-
    c(namesof("prob"     , .lprob      , earg = .eprob      , tag = FALSE),
      namesof("onempobs0", .lonempobs0 , earg = .eonempobs0 , tag = FALSE))


    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    phi.init <- if (length( .ionempobs0 )) 1 - .ionempobs0 else {
        prob0.est <- sum(Size[y == 0]) / sum(Size)
        if ( .imethod == 1) {
          (prob0.est - (1 - mustart)^Size) / (1 - (1 - mustart)^Size)
        } else
        if ( .imethod == 2) {
          prob0.est
        } else {
          prob0.est * 0.5
        }
    }

    phi.init[phi.init <= -0.10] <- 0.50  # Lots of sample variation
    phi.init[phi.init <=  0.01] <- 0.05  # Last resort
    phi.init[phi.init >=  0.99] <- 0.95  # Last resort




    if (!length(etastart)) {
      etastart <-
        cbind(theta2eta(     mustart, .lprob      , earg = .eprob      ),
              theta2eta(1 - phi.init, .lonempobs0 , earg = .eonempobs0 ))

      mustart <- NULL
    }
  }), list( .lprob = lprob, .lonempobs0 = lonempobs0,
            .eprob = eprob, .eonempobs0 = eonempobs0,
            .iprob = iprob, .ionempobs0 = ionempobs0,
            .imethod = imethod,
            .type.fitted = type.fitted ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "onempobs0"))[1]
    
    prob      <- eta2theta(eta[, 1], .lprob      , earg = .eprob  )
    onempobs0 <- eta2theta(eta[, 2], .lonempobs0 , earg = .eonempobs0 )
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    ans <- switch(type.fitted,
                  "mean"      = onempobs0 * prob / (1 - (1 - prob)^Size),
                  "pobs0"     = 1 - onempobs0,  # P(Y=0)
                  "onempobs0" =     onempobs0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lprob = lprob, .lonempobs0 = lonempobs0,
           .eprob = eprob, .eonempobs0 = eonempobs0 ))),

  last = eval(substitute(expression({
    misc$link <-    c(prob = .lprob, onempobs0 = .lonempobs0 )
    misc$earg <- list(prob = .eprob, onempobs0 = .eonempobs0 )

    misc$imethod  <- .imethod
    misc$zero     <- .zero
    misc$expected <- TRUE
  }), list( .lprob = lprob, .lonempobs0 = lonempobs0,
            .eprob = eprob, .eonempobs0 = eonempobs0,
            .zero = zero,
            .imethod = imethod ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w
    prob      <- eta2theta(eta[, 1], .lprob      , earg = .eprob      )
    onempobs0 <- eta2theta(eta[, 2], .lonempobs0 , earg = .eonempobs0 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        orig.w * dzabinom(x = round(y * Size), size = Size,
                          prob = prob, pobs0 = 1 - onempobs0,
                          log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .lonempobs0 = lonempobs0,
           .eprob = eprob, .eonempobs0 = eonempobs0 ))),
  vfamily = c("zabinomialff"),

  deriv = eval(substitute(expression({
    NOS <- if (length(extra$NOS)) extra$NOS else 1
    M1 <- 2

    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    prob      <- eta2theta(eta[, 1], .lprob      , earg = .eprob      )
    onempobs0 <- eta2theta(eta[, 2], .lonempobs0 , earg = .eonempobs0 )
    phi0 <- 1 - onempobs0

    dprob.deta      <- dtheta.deta(prob     , .lprob      ,
                                   earg = .eprob      )
    donempobs0.deta <- dtheta.deta(onempobs0, .lonempobs0 ,
                                   earg = .eonempobs0 )

    df0.dprob   <- -Size *              (1 -  prob)^(Size - 1)
    df02.dprob2 <-  Size * (Size - 1) * (1 -  prob)^(Size - 2)
    prob0  <- (1 -  prob)^(Size)
    oneminusf0  <- 1 - prob0


    dl.dprob <-  c(w)      * (y / prob - (1 - y) / (1 - prob)) +
                 c(orig.w) * df0.dprob / oneminusf0
    dl.donempobs0 <- +1 / (onempobs0)


    dl.donempobs0[y == 0] <-
      -1 / (1 - onempobs0[y == 0])  # Do it in 1 line
    skip <- extra$skip.these
    for (spp. in 1:NOS) {
      dl.dprob[skip[, spp.], spp.] <- 0
    }


    ans <- cbind(            dl.dprob      * dprob.deta,
                 c(orig.w) * dl.donempobs0 * donempobs0.deta)
                 
    ans
  }), list( .lprob = lprob, .lonempobs0 = lonempobs0,
            .eprob = eprob, .eonempobs0 = eonempobs0 ))),


  weight = eval(substitute(expression({
    wz <- matrix(0.0, n, M1)

    usualmeanY <-  prob
    meanY <- (1 - phi0) * usualmeanY / oneminusf0


    term1 <-  c(Size) * (meanY /      prob^2 -
                         meanY / (1 - prob)^2) +
             c(Size) * (1 - phi0) / (1 - prob)^2

    term2 <-  -(1 - phi0) * df02.dprob2 / oneminusf0
    term3 <-  -(1 - phi0) * (df0.dprob  / oneminusf0)^2
    ned2l.dprob2 <- term1 + term2 + term3
    wz[, iam(1, 1, M)] <- ned2l.dprob2 * dprob.deta^2


    mu.phi0 <- phi0
    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if (FALSE &&
                  .lonempobs0 == "logit" &&
                  is.empty.list( .eonempobs0 )) {
      tmp100
    } else {
      (donempobs0.deta^2) / tmp100
    }
    wz[, iam(2, 2, M)] <- tmp200


    c(orig.w) * wz
  }), list( .lprob = lprob, .lonempobs0 = lonempobs0,
            .eprob = eprob, .eonempobs0 = eonempobs0 ))))
}






 zageometric <- function(lpobs0 = "logit", lprob = "logit",
                         type.fitted = c("mean", "pobs0", "onempobs0"),
                         imethod = 1,
                         ipobs0 = NULL, iprob = NULL,
                         zero = NULL) {



  lpobs0 <- as.list(substitute(lpobs0))
  epobs0 <- link2list(lpobs0)
  lpobs0 <- attr(epobs0, "function.name")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "onempobs0"))[1]


  if (!is.Numeric(imethod, length.arg = 1,
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
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 2
    if (any(y < 0))
      stop("the response must not have negative values")

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)

    extra$dimnamesy <- dimnames(y)
    extra$type.fitted      <- .type.fitted

    
    mynames1 <- if (ncoly == 1) "pobs0"  else
                paste("pobs0",  1:ncoly, sep = "")
    mynames2 <- if (ncoly == 1) "prob" else
                paste("prob", 1:ncoly, sep = "")
    predictors.names <-
        c(namesof(mynames1, .lpobs0 , earg = .epobs0 , tag = FALSE),
          namesof(mynames2, .lprob  , earg = .eprob  , tag = FALSE))[
          interleave.VGAM(M1*NOS, M = M1)]

    if (!length(etastart)) {

      foo <- function(x) mean(as.numeric(x == 0))
      phi0.init <- matrix(apply(y, 2, foo), n, ncoly, byrow = TRUE)
      if (length( .ipobs0 ))
        phi0.init <- matrix( .ipobs0 , n, ncoly, byrow = TRUE)


      prob.init <-
        if ( .imethod == 2)
          1 / (1 + y + 1/16) else
        if ( .imethod == 1)
          (1 - phi0.init) / (1 +
          matrix(colSums(y * w) / colSums(w) + 1/16,
                 n, ncoly, byrow = TRUE)) else
          (1 - phi0.init) / (1 +
          matrix(apply(y, 2, median), n, ncoly, byrow = TRUE) + 1/16)


      if (length( .iprob ))
        prob.init <- matrix( .iprob , n, ncoly, byrow = TRUE)



      etastart <- cbind(theta2eta(phi0.init, .lpobs0 , earg = .epobs0 ),
                       theta2eta(prob.init, .lprob ,  earg = .eprob ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = M1)]
    }
  }), list( .lpobs0 = lpobs0, .lprob = lprob,
            .epobs0 = epobs0, .eprob = eprob,
            .ipobs0 = ipobs0, .iprob = iprob,
            .imethod = imethod,
            .type.fitted = type.fitted ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "onempobs0"))[1]

    NOS <- extra$NOS
    M1 <- 2

    phi0 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                             .lpobs0 , earg = .epobs0 ))
    prob <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                             .lprob  , earg = .eprob ))


    ans <- switch(type.fitted,
                  "mean"      = (1 - phi0) / prob,
                  "pobs0"     =      phi0,  # P(Y=0)
                  "onempobs0" =  1 - phi0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpobs0 = lpobs0, .lprob = lprob,
           .epobs0 = epobs0, .eprob = eprob ))),
  last = eval(substitute(expression({
    temp.names <- c(rep( .lpobs0 , len = NOS),
                    rep( .lprob  , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M = M1)]
    misc$link  <- temp.names

    misc$earg <- vector("list", M1 * NOS)

    names(misc$link) <-
    names(misc$earg) <-
        c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M = M1)]

    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epobs0
      misc$earg[[M1*ii  ]] <- .eprob
    }


    misc$expected <- TRUE
    misc$imethod <- .imethod
    misc$ipobs0  <- .ipobs0
    misc$iprob   <- .iprob
    misc$multipleResponses <- TRUE
  }), list( .lpobs0 = lpobs0, .lprob = lprob,
            .epobs0 = epobs0, .eprob = eprob,
            .ipobs0 = ipobs0, .iprob = iprob,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    NOS <- extra$NOS
    M1 <- 2

    phi0 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                            .lpobs0 , earg = .epobs0 ))
    prob <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                            .lprob  , earg = .eprob  ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzageom(x = y, pobs0 = phi0, prob = prob, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs0 = lpobs0, .lprob = lprob,
           .epobs0 = epobs0, .eprob = eprob ))),
  vfamily = c("zageometric"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    phi0 <- cbind(eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                            .lpobs0 , earg = .epobs0 ))
    prob <- cbind(eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                            .lprob  , earg = .eprob  ))
    rzageom(nsim * length(prob), prob = prob, pobs0 = phi0)
  }, list( .lpobs0 = lpobs0, .lprob = lprob,
           .epobs0 = epobs0, .eprob = eprob ))),




  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- extra$NOS
    y0 <- extra$y0
    skip <- extra$skip.these

    phi0 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                            .lpobs0 , earg = .epobs0 ))
    prob <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                            .lprob  , earg = .eprob  ))


    dl.dprob <-  1 / prob - (y - 1) / (1 - prob)
    dl.dphi0 <- -1 / (1 - phi0)


    for (spp. in 1:NOS) {
      dl.dphi0[skip[, spp.], spp.] <- 1 / phi0[skip[, spp.], spp.]
      dl.dprob[skip[, spp.], spp.] <- 0
    }
    dphi0.deta <- dtheta.deta(phi0, .lpobs0 , earg = .epobs0 )
    dprob.deta <- dtheta.deta(prob, .lprob  , earg = .eprob  )


    ans <- c(w) * cbind(dl.dphi0 * dphi0.deta,
                        dl.dprob * dprob.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M = M1)]
    ans
  }), list( .lpobs0 = lpobs0, .lprob = lprob,
            .epobs0 = epobs0, .eprob = eprob ))),
  weight = eval(substitute(expression({

    wz <- matrix(0.0, n, M1*NOS)


    ned2l.dprob2 <- (1 - phi0) / (prob^2 * (1 - prob))

    wz[, NOS+(1:NOS)] <- c(w) * ned2l.dprob2 * dprob.deta^2


    mu.phi0 <- phi0
    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if ( .lpobs0 == "logit" && is.empty.list( .epobs0 )) {
      cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (dphi0.deta^2) / tmp100)
    }
    wz[, 1:NOS] <-  tmp200


    wz <- wz[, interleave.VGAM(ncol(wz), M = M1)]


    wz
  }), list( .lpobs0 = lpobs0,
            .epobs0 = epobs0 ))))
}  # End of zageometric




 zageometricff <- function(lprob = "logit", lonempobs0 = "logit",
                           type.fitted = c("mean", "pobs0", "onempobs0"),
                           imethod = 1,
                           iprob = NULL, ionempobs0 = NULL,
                           zero = -2) {


  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lonempobs0 <- as.list(substitute(lonempobs0))
  eonempobs0 <- link2list(lonempobs0)
  lonempobs0 <- attr(eonempobs0, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "pobs0", "onempobs0"))[1]


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
       max(iprob) >= 1)
    stop("argument 'iprob' out of range")

  if (length(ionempobs0))
    if (!is.Numeric(ionempobs0, positive = TRUE) ||
       max(ionempobs0) >= 1)
      stop("argument 'ionempobs0' out of range")


  new("vglmff",
  blurb = c("Zero-altered geometric ",
            "(Bernoulli and positive-geometric conditional model)\n\n",
            "Links:    ",
            namesof("prob"     , lprob     , earg = eprob     , tag = FALSE), ", ",
            namesof("onempobs0", lonempobs0, earg = eonempobs0, tag = FALSE), "\n",
            "Mean:     onempobs0 / prob"),

  constraints = eval(substitute(expression({

    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 2
    if (any(y < 0))
      stop("the response must not have negative values")

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)

    extra$dimnamesy   <- dimnames(y)
    extra$type.fitted <- .type.fitted

    
    mynames1 <- if (ncoly == 1) "prob"       else
                paste("prob",       1:ncoly, sep = "")
    mynames2 <- if (ncoly == 1) "onempobs0"  else
                paste("onempobs0",  1:ncoly, sep = "")
    predictors.names <-
        c(namesof(mynames1, .lprob      , earg = .eprob      , tag = FALSE),
          namesof(mynames2, .lonempobs0 , earg = .eonempobs0 , tag = FALSE))[
          interleave.VGAM(M1*NOS, M = M1)]

    if (!length(etastart)) {

      foo <- function(x) mean(as.numeric(x == 0))
      phi0.init <- matrix(apply(y, 2, foo), n, ncoly, byrow = TRUE)
      if (length( .ionempobs0 ))
        phi0.init <- matrix( 1 - .ionempobs0 , n, ncoly, byrow = TRUE)


      prob.init <-
        if ( .imethod == 2)
          1 / (1 + y + 1/16) else
        if ( .imethod == 1)
          (1 - phi0.init) / (1 +
          matrix(colSums(y * w) / colSums(w) + 1/16,
                 n, ncoly, byrow = TRUE)) else
          (1 - phi0.init) / (1 +
          matrix(apply(y, 2, median), n, ncoly, byrow = TRUE) + 1/16)


      if (length( .iprob ))
        prob.init <- matrix( .iprob , n, ncoly, byrow = TRUE)



      etastart <-
        cbind(theta2eta(    prob.init, .lprob      , earg = .eprob      ),
              theta2eta(1 - phi0.init, .lonempobs0 , earg = .eonempobs0 ))
                        
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = M1)]
    }
  }), list( .lonempobs0 = lonempobs0, .lprob = lprob,
            .eonempobs0 = eonempobs0, .eprob = eprob,
            .ionempobs0 = ionempobs0, .iprob = iprob,
            .imethod = imethod,
            .type.fitted = type.fitted ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "pobs0", "onempobs0"))[1]

    NOS <- extra$NOS
    M1 <- 2

    prob      <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                                 .lprob  , earg = .eprob ))
    onempobs0 <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                                 .lonempobs0 , earg = .eonempobs0 ))


    ans <- switch(type.fitted,
                  "mean"      =     (onempobs0) / prob,
                  "pobs0"     =  1 - onempobs0,  # P(Y=0)
                  "onempobs0" =      onempobs0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempobs0 = lonempobs0, .lprob = lprob,
           .eonempobs0 = eonempobs0, .eprob = eprob ))),
  last = eval(substitute(expression({
    temp.names <- c(rep( .lprob      , len = NOS),
                    rep( .lonempobs0 , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M = M1)]
    misc$link  <- temp.names

    misc$earg <- vector("list", M1 * NOS)

    names(misc$link) <-
    names(misc$earg) <-
        c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M = M1)]

    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .eprob
      misc$earg[[M1*ii  ]] <- .eonempobs0
    }


    misc$expected <- TRUE
    misc$imethod <- .imethod
    misc$ionempobs0  <- .ionempobs0
    misc$iprob   <- .iprob
    misc$multipleResponses <- TRUE
  }), list( .lonempobs0 = lonempobs0, .lprob = lprob,
            .eonempobs0 = eonempobs0, .eprob = eprob,
            .ionempobs0 = ionempobs0, .iprob = iprob,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    NOS <- extra$NOS
    M1 <- 2

    prob      <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                                 .lprob      , earg = .eprob      ))
    onempobs0 <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                                 .lonempobs0 , earg = .eonempobs0 ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzageom(x = y, pobs0 = 1 - onempobs0, prob = prob,
                       log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempobs0 = lonempobs0, .lprob = lprob,
           .eonempobs0 = eonempobs0, .eprob = eprob ))),
  vfamily = c("zageometricff"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    onempobs0 <- cbind(eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                                 .lonempobs0 , earg = .eonempobs0 ))
    prob      <- cbind(eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                                 .lprob      , earg = .eprob      ))
    rzageom(nsim * length(prob), pobs0 = 1 - onempobs0, prob = prob)
  }, list( .lonempobs0 = lonempobs0, .lprob = lprob,
           .eonempobs0 = eonempobs0, .eprob = eprob ))),




  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- extra$NOS
    y0 <- extra$y0
    skip <- extra$skip.these

    prob      <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                       .lprob      , earg = .eprob      ))
    onempobs0 <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                       .lonempobs0 , earg = .eonempobs0 ))
    pobs0 <- 1 - onempobs0


    dl.dprob      <-  1 / prob - (y - 1) / (1 - prob)
    dl.donempobs0 <- +1 / (onempobs0)


    for (spp. in 1:NOS) {
      dl.donempobs0[skip[, spp.], spp.] <- -1 / pobs0[skip[, spp.], spp.]
      dl.dprob[skip[, spp.], spp.] <- 0
    }
    dprob.deta      <- dtheta.deta(prob,      .lprob  , earg = .eprob  )
    donempobs0.deta <- dtheta.deta(onempobs0, .lonempobs0 ,
                                   earg = .eonempobs0 )


    ans <- c(w) * cbind(dl.dprob      * dprob.deta,
                        dl.donempobs0 * donempobs0.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M = M1)]
    ans
  }), list( .lonempobs0 = lonempobs0, .lprob = lprob,
            .eonempobs0 = eonempobs0, .eprob = eprob ))),
  weight = eval(substitute(expression({

    wz <- matrix(0.0, n, M1*NOS)


    ned2l.dprob2 <- (1 - pobs0) / (prob^2 * (1 - prob))

    wz[, (1:NOS)] <- c(w) * ned2l.dprob2 * dprob.deta^2


    mu.phi0 <- pobs0  # phi0
    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if ( FALSE &&
                  .lonempobs0 == "logit" &&
                  is.empty.list( .eonempobs0 )) {

      cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (donempobs0.deta^2) / tmp100)
    }
    wz[, NOS+(1:NOS)] <- tmp200


    wz <- wz[, interleave.VGAM(ncol(wz), M = M1)]


    wz
  }), list( .lonempobs0 = lonempobs0,
            .eonempobs0 = eonempobs0 ))))
}  # End of zageometricff






